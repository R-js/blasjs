/*>  
    -- Jacob Bogers, 03/2018,
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/


import { Complex, errMissingIm, errWrongArg, FortranArr } from '../../f_func';



export function chpr2(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    ap: FortranArr
): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('ap.i'));
    }

    // faster then String.toLowerCase()
    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0X20);

    let info = 0;
    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (incy === 0) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chpr2', info));
    }
    const { re: AlphaRe, im: AlphaIm } = alpha;
    const alphaIsZero = AlphaRe === 0 && AlphaIm === 0;

    if (n === 0 || alphaIsZero) return;

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;
    let ky = incy < 0 ? 1 - (n - 1) * incy : 1;

    let jx = kx - x.base;
    let jy = ky - y.base;

    let kk = 1;

    if (ul === 'u') {
        for (let j = 1; j <= n; j++) {
            const apkk = kk - ap.base; //note kk is updated within this loop
            const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
            const yIsZero = y.r[jy] === 0 && y.i[jy] === 0;
            if (!xIsZero || !yIsZero) {
                //(a+ib)(c-id)=(ac+bd)+i(-ad+ bc)
                let temp1Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp1Im = -AlphaRe * y.i[jy] + AlphaIm * y.r[jy];

                //(a-ib)(c+id)=(ac+bd)+i(ad-bc)

                let temp2Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp2Im = AlphaRe * y.i[jy] - AlphaIm * y.i[jy];


                ap.i[apkk] = 0;
                ap.r[apkk] += (x.r[jx] * temp1Re - x.i[jx] * temp1Im)
                    + (y.r[jy] * temp2Re - y.i[jy] * temp2Im);

                let ix = jx;
                let iy = jy;
                for (let k = kk; k <= kk + j - 2; k++) {
                    ap.r[k - ap.base] +=
                        (x.r[ix] * temp1Re - x.i[ix] * temp1Im) +
                        (y.r[iy] * temp2Re - y.i[iy] * temp2Im);
                    ix += incx;
                    iy += incy;
                }
                ap.i[apkk + j - 1] = 0;
                ap.r[apkk + j - 1] =
                    // REAL(x[jx]*temp1)    
                    (x.r[jx] * temp1Re - x.i[jx] * temp1Im)
                    + //REAL(y[jy]*temp2)
                    (y.r[jy] * temp2Re - y.i[jy] * temp2Im);

            } //if  !xIsZero || !yIsZero
            else {
                // nuke imaginary part
                ap.i[apkk] = 0;
            }// if
            jx += incx;
            jy += incy;
            kk += j;
        }//for
    }
    else {
        // Form  A  when lower triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            const apkk = kk - ap.base; //note kk is updated within this loop
            const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
            const yIsZero = y.r[jy] === 0 && y.i[jy] === 0;

            if (!(xIsZero && yIsZero)) {

                //(a+ib)(c-id)=(ac+bd)+i(-ad+ bc)
                let temp1Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp1Im = -AlphaRe * y.i[jy] + AlphaIm * y.r[jy];

                //(a-ib)(c+id)=(ac+bd)+i(ad-bc)

                let temp2Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp2Im = AlphaRe * y.i[jy] - AlphaIm * y.i[jy];

                ap.i[apkk] = 0;
                ap.r[apkk] =
                    // REAL(x[jx]*temp1)    
                    (x.r[jx] * temp1Re - x.i[jx] * temp1Im)
                    + //REAL(y[jy]*temp2)
                    (y.r[jy] * temp2Re - y.i[jy] * temp2Im);

                let ix = jx;
                let iy = jy;

                for (let k = kk + 1; k <= kk + n - j; k++) {

                    ix += incx;
                    iy += incy;
                    ap.r[k - ap.base] +=
                        (x.r[ix] * temp1Re - x.i[ix] * temp1Im)
                        +
                        (y.r[iy] * temp2Re - y.i[iy] * temp2Re);
                }
            } else {
                ap.i[apkk] = 0;
            }
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}

