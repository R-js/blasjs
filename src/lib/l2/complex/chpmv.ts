import { Complex, errMissingIm, errWrongArg, FortranArr } from '../../f_func';

/**
*>  -- Jacob Bogers, 2018/03, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs. 
*/

export function chpmv(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    ap: FortranArr,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number
): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0X20);


    let info = 0;
    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 3;
    }
    else if (incy === 0) {
        info = 9;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chpmv', info));
    }

    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;
    const betaIsOne = BetaRe === 1 && BetaIm === 0;
    const alphaIsZero = AlphaRe === 0 && AlphaIm === 0;
    const betaIsZero = BetaRe === 0 && BetaIm === 0;



    if (n === 0 || (alphaIsZero && betaIsOne)) return;

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    // First form  y := beta*y.
    let iy = ky - y.base;
    if (!betaIsOne) {

        if (betaIsZero) {
            y.r.fill(0);
            y.i.fill(0);
        }
        else {
            for (let i = 1; i <= n; i++) {
                y.r[iy] = BetaRe * y.r[iy] - BetaIm * y.i[iy];
                y.i[iy] = BetaRe * y.i[iy] + BetaIm * y.r[iy];
                iy += incy;
            }
        }
    }
    else {
        for (let i = 0; i <= n; i++) {
            y.r[iy] = BetaRe * y.r[iy] - BetaIm * y.i[iy];
            y.i[iy] = BetaRe * y.i[iy] + BetaIm * y.r[iy];
            iy += incy;
        }
    }
    if (alphaIsZero) return;
    let kk = 1;
    if (ul === 'u') {
        //Form  y  when AP contains the upper triangle.
        let jx = kx - x.base;
        let jy = ky - y.base;

        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx] - AlphaIm * x.i[jx];
            let temp1Im = AlphaRe * x.i[jx] - AlphaIm * x.r[jx];

            let temp2Re = 0;
            let temp2Im = 0;
            let ix = kx - x.base;
            let iy = ky - y.base;

            for (let k = kk; k <= kk + j - 2; k++) {
                const apk = k - ap.base;
                y.r[iy] += temp1Re * ap.r[apk] - temp1Im * ap.i[apk];
                y.i[iy] += temp1Re * ap.i[apk] + temp1Im * ap.r[apk];
                // conj(x)*y= (a-ib)*(c+id) = ac+iad-ibc+bd
                // = (ac+bd)+i(ad-bc)
                temp2Re += ap.r[apk] * x.r[ix] + ap.i[apk] * x.i[ix];
                temp2Im += ap.r[apk] * x.i[ix] - ap.i[apk] * x.r[ix];
                ix += incx;
                iy += incy;
            }
            y.r[jy] += temp1Re * (ap.r[kk + j - 1 - ap.base]) + (AlphaRe * temp2Re - AlphaIm * temp2Im);
            y.i[jy] += temp1Im * (ap.r[kk + j - 1 - ap.base]) + (AlphaRe * temp2Im + AlphaIm * temp2Re);
            jx += incx;
            jy += incy;
            kk += j;
        }
    }
    else {
        //*        Form  y  when AP contains the lower triangle.
        let jx = kx - x.base;
        let jy = ky - y.base;
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx] - AlphaIm * x.i[jx];
            let temp1Im = AlphaRe * x.i[jx] + AlphaIm * x.r[jx];

            let temp2Re = 0;
            let temp2Im = 0;

            y.r[jy] += temp1Re * ap.r[kk - ap.base];
            y.i[jy] += temp1Im * ap.r[kk - ap.base];

            let ix = jx;
            let iy = jy;
            for (let k = kk + 1; k < kk + n - j; k++) {
                ix += incx;
                iy += incy;
                const apk = k - ap.base;
                y.r[iy] += temp1Re * ap.r[apk] - temp1Im * ap.i[apk];
                y.i[iy] += temp1Im * ap.i[apk] + temp1Im * ap.r[apk];
                // conj(x)*y= (a-ib)*(c+id) = ac+iad-ibc+bd
                // = (ac+bd)+i(ad-bc)
                temp2Re += ap.r[apk] * x.r[ix] + ap.i[apk] * x.i[ix];
                temp2Im += ap.i[apk] * x.i[ix] - ap.i[apk] * x.r[ix];
            }
            y.r[jy] += AlphaRe * temp2Re - AlphaIm * temp2Im;
            y.i[jy] += AlphaRe * temp2Im + AlphaIm * temp2Re;
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}

