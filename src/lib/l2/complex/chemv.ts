/*
*>    -- Jacob Bogers, JS Port,  03/2018,jkfbogers@gmail.com
     
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> CHEMV  performs the matrix-vector  operation
*>
*>    y := alpha*A*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are n element vectors and
*> A is an n by n hermitian matrix.
*/

import { Complex, errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';

const { max } = Math;

export function chemv(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    a: Matrix2D,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number): void {

    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0x20);

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }
    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    let info = 0;
    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (lda < max(1, n)) {
        info = 5;
    }
    else if (incx === 0) {
        info = 7;
    }
    else if (incy === 0) {
        info = 10;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('chemv', info));
    }


    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;

    const alphaIsZero = AlphaRe === 0 && AlphaIm === 0;
    const betaIsOne = BetaRe === 1 && BetaIm === 0;
    const betaIsZero = BetaRe === 0 && BetaRe === 0;

    //*     Quick return if possible.
    if (n === 0 || (alphaIsZero && betaIsOne)) return;

    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;



    if (betaIsOne) {
        if (incy === 1) {
            y.r.fill(0);
            y.i.fill(0);
        }
        else {
            let iy = ky;
            for (let i = 1; i <= n; i++) {
                y.r[iy - y.base] = betaIsZero ? 0 : BetaRe * y.r[iy - y.base] - BetaIm * y.i[iy - y.base];
                y.i[iy - y.base] = betaIsZero ? 0 : BetaRe * y.i[iy - y.base] + BetaIm * y.r[iy - y.base];
                iy += incy;
            }
        }
    }

    if (alphaIsZero) return; //done
    if (ul === 'u') {

        // Form  y  when A is stored in upper triangle.

        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.i[jx - x.base];
            let temp2Re = 0;
            let temp2Im = 0;
            let ix = kx;
            let iy = ky;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= j - 1; i++) {
                y.r[iy - y.base] += temp1Re * a.r[coords + i] - temp1Im * a.i[coords + i];
                y.i[iy - y.base] += temp1Re * a.i[coords + i] + temp1Im * a.r[coords + i];
                temp2Re += a.r[coords + i] * x.r[ix - x.base] + a.i[coords + i] * x.i[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1Re * a.r[coords + j] + AlphaRe * temp2Re - AlphaIm * temp2Im;
            y.i[jy - y.base] += temp1Im * a.r[coords + j] + AlphaRe * temp2Im + AlphaIm * temp2Re;
            jx += incx;
            jy += incy;
        }
    }
    else {
        //Form  y  when A is stored in lower triangle.
        let jx = kx - x.base;
        let jy = ky - y.base;
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx] - AlphaIm * x.i[jx]
            let temp1Im = AlphaRe * x.i[jx] + AlphaIm * x.r[jx];

            const coords = a.colOfEx(j);
            y.r[jy] += temp1Re * a.r[coords + j] - temp1Im * a.i[coords + j];
            y.i[jy] += temp1Re * a.i[coords + j] - temp1Im * a.r[coords + j];
            let ix = jx;
            let iy = jy;
            let temp2Re = 0;
            let temp2Im = 0;
            for (let i = j + 1; j <= n; j++) {
                ix += incx;
                iy += incy;
                y.r[iy] += AlphaRe * temp1Re - AlphaIm * temp1Im;
                y.i[iy] += AlphaRe * temp1Im + AlphaIm * temp1Re;

                temp2Re += a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                temp2Im += a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
            }
            y.r[jy] += AlphaRe * temp2Re - AlphaRe * temp2Im;
            y.i[jy] += AlphaRe * temp2Im + AlphaIm * temp2Re;
            jx += incx;
            jy += incy;
        }
    }
}



