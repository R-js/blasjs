/*

*  -- Jacob Bogers, JS port, 03/2018, jkfbogers@gmail.com
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*/

/*
*>
*> CGERC  performs the rank 1 operation
*>
*>    A := alpha*x*y**H + A,
*>
*> where alpha is a scalar, x is an m element vector, y is an n element
*> vector and A is an m by n matrix.
*/


import {
    Complex,
    errMissingIm,
    errWrongArg,
    FortranArr,
    isOne,
    isZero,
    lowerChar,
    Matrix
} from '../../f_func';

const { max, min } = Math;

export function chbmv(
    uplo: 'u' | 'l',
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number
): void {


    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    const ul = lowerChar(uplo);

    let info = 0;

    if (!'ul'.includes(uplo)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (k < 0) {
        info = 3;
    }
    else if (lda < (k + 1)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    else if (incy === 0) {
        info = 11;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('chbmv', info));
    }

    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;

    const alphaIsZero = isZero(alpha);
    const betaIsOne = isOne(beta);
    const betaIsZero = isZero(beta);

    if (n === 0 || (alphaIsZero && betaIsOne)) return;

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    if (!betaIsOne) {
        let iy = ky;
        for (let i = 1; i <= n; i++) {
            const re = betaIsZero ? 0 : BetaRe * y.r[iy - y.base] - BetaIm * y.i[iy - y.base];
            const im = betaIsZero ? 0 : BetaRe * y.i[iy - y.base] + BetaIm * y.r[iy - y.base];
            y.r[iy - y.base] = re;
            y.i[iy - y.base] = im;
            iy += incy;
        }

    }
    if (alphaIsZero) return;

    let jx = kx;
    let jy = ky;

    if (ul === 'u') {
        // Form  y  when upper triangle of A is stored.
        let kplus1 = k + 1;
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];

            let temp2Re = 0;
            let temp2Im = 0;

            let ix = kx;
            let iy = ky;

            let L = kplus1 - j;
            //
            let coords = a.colOfEx(j);
            //
            for (let i = max(1, j - k); i <= j - 1; i++) {
                y.r[iy - y.base] += temp1Re * a.r[coords + L + i] - temp1Im * a.i[coords + L + i];
                y.i[iy - y.base] += temp1Re * a.i[coords + L + i] + temp1Im * a.r[coords + L + i];
                // CONJG(A(L+I,J))*X(IX)
                //(a - ib) * (c +id)=(ac+bd)+i(ad-bc)
                temp2Re += a.r[coords + L + i] * x.r[ix - x.base] + a.i[coords + L + i] * x.i[ix - x.base];
                temp2Im += a.r[coords + L + i] * x.i[ix - x.base] - a.i[coords + L + i] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[iy - y.base] += temp1Re * a.r[coords + kplus1] + (AlphaRe * temp2Re - AlphaIm * temp2Im);
            y.i[iy - y.base] += temp1Im * a.r[coords + kplus1] + (AlphaRe * temp2Im + AlphaIm * temp2Re);

            jx += incx;
            jy += incy;
            if (j > k) {
                kx += incx;
                ky += incy;
            }
        }
    }
    else {
        //Form  y  when lower triangle of A is stored.
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];

            let temp2Re = 0;
            let temp2Im = 0;

            const coords = a.colOfEx(j);
            y.r[jy - y.base] += temp1Re * a.r[coords + 1];
            y.i[jy - y.base] += temp1Im * a.r[coords + 1];

            let L = 1 - j;
            let ix = jx;
            let iy = jy;

            for (let i = j + 1; i <= min(n, j + k); i++) {
                ix += incx;
                iy += incy;

                y.r[iy - y.base] += temp1Re * a.r[coords + L + i] - temp1Im * a.i[coords + L + i];
                y.i[iy - y.base] += temp1Re * a.i[coords + L + i] + temp1Im * a.r[coords + L + i];

                // conjugate!!!
                temp2Re += a.r[coords + L + i] * x.r[ix - x.base] + a.i[coords + L + i] * x.i[ix - x.base];
                temp2Im += a.r[coords + L + i] * x.i[ix - x.base] - a.i[coords + L + i] * x.r[ix - x.base];
            }
            y.r[jy - y.base] += AlphaRe * temp2Re - AlphaIm * temp2Im;
            y.i[jy - y.base] += AlphaRe * temp2Im + AlphaIm * temp2Re;

            jx += incx;
            jy += incy;
        }
    }
}
