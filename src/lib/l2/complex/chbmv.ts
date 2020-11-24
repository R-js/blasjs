/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import { Complex, errMissingIm, errWrongArg, FortranArr, isOne, isZero, lowerChar, Matrix } from '../../f_func';

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
    incy: number,
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
    } else if (n < 0) {
        info = 2;
    } else if (k < 0) {
        info = 3;
    } else if (lda < k + 1) {
        info = 6;
    } else if (incx === 0) {
        info = 8;
    } else if (incy === 0) {
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
        const kplus1 = k + 1;
        for (let j = 1; j <= n; j++) {
            const temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            const temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];

            let temp2Re = 0;
            let temp2Im = 0;

            let ix = kx;
            let iy = ky;

            const L = kplus1 - j;
            //
            const coords = a.colOfEx(j);
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
    } else {
        //Form  y  when lower triangle of A is stored.
        for (let j = 1; j <= n; j++) {
            const temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            const temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];

            let temp2Re = 0;
            let temp2Im = 0;

            const coords = a.colOfEx(j);
            y.r[jy - y.base] += temp1Re * a.r[coords + 1];
            y.i[jy - y.base] += temp1Im * a.r[coords + 1];

            const L = 1 - j;
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
