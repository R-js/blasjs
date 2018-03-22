/*
*>  -- Jacob Boggers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

import { Complex, errMissingIm, errWrongArg, FortranArr, Matrix } from '../../f_func';

const { max } = Math;
/*
*> CHER2  performs the hermitian rank 2 operation
*>
*>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
*>
*> where alpha is a scalar, x and y are n element vectors and A is an n
*> by n hermitian matrix.
*/

export function cher2(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    a: Matrix,
    lda: number): void {


    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (a.i === undefined) {
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
        info = 5;
    }
    else if (incy === 0) {
        info = 7;
    }
    else if (lda < max(1, n)) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('cher2', info));
    }

    const { re: AlphaRe, im: AlphaIm } = alpha;

    const alphaIsZero = AlphaRe === 0 && AlphaIm === 0;

    if (n === 0 || alphaIsZero) return; //nothing to do

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;
    let jx = kx - x.base;
    let jy = ky - y.base;

    if (ul === 'u') {
        //        Form  A  when A is stored in the upper triangle.
        for (let j = 1; j <= n; j++) {
            const coords = a.colOfEx(j);
            if (!(x.r[jx] === 0 && x.i[jx] === 0)
                || !(y.r[jy] === 0 && y.i[jy] === 0)
            ) {
                let temp1Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp1Im = -AlphaRe * y.i[jy] + AlphaIm * y.r[jy];

                let temp2Re = AlphaRe * x.r[jx] - AlphaIm * x.i[jx];
                let temp2Im = AlphaRe * x.i[jx] + AlphaIm * x.r[jx];

                let ix = kx - x.base;
                let iy = ky - y.base;

                for (let i = 1; i <= j - 1; i++) {
                    a.r[coords + 1] +=
                        (x.r[ix] * temp1Re - x.i[ix] * temp1Im) +
                        (y.r[iy] * temp2Re - y.i[iy] * temp2Im);
                    a.i[coords + 1] +=
                        (x.r[ix] * temp1Im + x.i[ix] * temp1Re) +
                        (y.r[iy] * temp2Im + y.i[iy] * temp2Re);

                    ix += incx;
                    iy += incy;
                }//for
                a.i[coords + j] = 0;
                a.r[coords + j] +=
                    (x.r[jx] * temp1Re - x.i[jx] * temp1Im)
                    +
                    (y.r[jy] * temp2Re - y.i[jy] * temp2Im)
            } else {
                a.i[coords + j] = 0;
            }
            jx += incx;
            jy += incy;
        }//for
    }
    else {
        //  Form  A  when A is stored in the lower triangle.

        for (let j = 1; j <= n; j++) {

            const coords = a.colOfEx(j);
            if (
                !(x.r[jx] === 0 && x.r[jx] === 0)
                ||
                !(y.r[jy] === 0 && x.i[jy] === 0)) {
                //TEMP1 = ALPHA*CONJG(Y(JY))
                // (a+ib)*(c-di)= ac -iad + ibc+ bd = (ac+bd) + i(-ad+bc)
                // (a-ib)*(c+di)= ac+adi -ibc+bd = (ac+bd)+i(ad-bc)

                let temp1Re = AlphaRe * y.r[jy] + AlphaIm * y.i[jy];
                let temp1Im = -AlphaRe * y.i[jy] + AlphaIm * y.r[jy];

                let temp2Re = AlphaRe * y.r[jy] - AlphaIm * y.i[jy];
                let temp2Im = -(AlphaRe * y.i[jy] + AlphaIm * y.i[jy]);

                a.i[coords + j] = 0;
                a.r[coords + j] +=
                    (x.r[jx] * temp1Re - x.i[jx] * temp1Im) +
                    (y.r[jy] * temp2Re - y.i[jy * temp2Im]);
                let ix = jx;
                let iy = jy;
                for (let i = j + 1; i <= n; i++) {
                    ix += incx;
                    iy += incy;
                    a.r[coords + i] +=
                        (x.r[ix] * temp1Re - x.i[ix] * temp1Im) +
                        (y.r[iy] * temp2Re - y.i[iy] * temp2Im);
                    a.i[coords + i] +=
                        (x.r[ix] * temp1Im + x.i[ix] * temp1Re) +
                        (y.r[iy] * temp2Im + y.i[iy] * temp2Re);
                }
            }
            else {
                a.i[coords + j] = 0;
            }
            jx += incx;
            jy += incy;
        }

    }
}
