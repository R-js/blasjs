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

const { max } = Math;

export function chemv(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number): void {

    const ul = lowerChar(uplo);

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
    if (!'ul'.includes(ul)) {
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

    const alphaIsZero = isZero(alpha);
    const betaIsOne = isOne(beta);
    const betaIsZero = isZero(beta);

    //*     Quick return if possible.
    if (n === 0 || (alphaIsZero && betaIsOne)) return;

    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;



    if (!betaIsOne) {

        let iy = ky;
        for (let i = 1; i <= n; i++) {
            const re = betaIsZero ? 0 : BetaRe * y.r[iy - y.base] - BetaIm * y.i[iy - y.base];
            const im = betaIsZero ? 0 : BetaRe * y.i[iy - y.base] + BetaIm * y.r[iy - y.base];
            y.r[iy - y.base] = re;
            y.i[iy - y.base] = im;
            iy += incy;
        } y.r[iy - y.base]
    }

    if (alphaIsZero) return; //done
    if (ul === 'u') {

        // Form  y  when A is stored in upper triangle.

        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];
            let temp2Re = 0;
            let temp2Im = 0;
            let ix = kx;
            let iy = ky;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= j - 1; i++) {
                y.r[iy - y.base] += temp1Re * a.r[coords + i] - temp1Im * a.i[coords + i];
                y.i[iy - y.base] += temp1Re * a.i[coords + i] + temp1Im * a.r[coords + i];
                // TEMP2 = TEMP2 + CONJG(A(I,J))*X(IX)
                // (a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                temp2Re += a.r[coords + i] * x.r[ix - x.base] + a.i[coords + i] * x.i[ix - x.base];
                temp2Im += a.r[coords + i] * x.i[ix - x.base] - a.i[coords + i] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            //  Y(JY) = Y(JY) + TEMP1*REAL(A(J,J)) + ALPHA*TEMP2
            y.r[jy - y.base] += temp1Re * a.r[coords + j] + (AlphaRe * temp2Re - AlphaIm * temp2Im);
            y.i[jy - y.base] += temp1Im * a.r[coords + j] + (AlphaRe * temp2Im + AlphaIm * temp2Re);
            jx += incx;
            jy += incy;
        }
    }
    else {
        //Form  y  when A is stored in lower triangle.
        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);

            let temp1Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let temp1Im = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];

            let temp2Re = 0;
            let temp2Im = 0;

            //Y(JY) = Y(JY) + TEMP1*REAL(A(J,J))
            y.r[jy - y.base] += temp1Re * a.r[coorAJ + j];
            y.i[jy - y.base] += temp1Im * a.r[coorAJ + j];

            let ix = jx;
            let iy = jy;

            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                iy += incy;
                //Y(IY) = Y(IY) + TEMP1*A(I,J)
                y.r[iy - y.base] += temp1Re * a.r[coorAJ + i] - temp1Im * a.i[coorAJ + i];
                y.i[iy - y.base] += temp1Re * a.i[coorAJ + i] + temp1Im * a.r[coorAJ + i];
                //TEMP2 = TEMP2 + CONJG(A(I,J))*X(IX)
                temp2Re += a.r[coorAJ + i] * x.r[ix - x.base] + a.i[coorAJ + i] * x.i[ix - x.base];
                temp2Im += a.r[coorAJ + i] * x.i[ix - x.base] - a.i[coorAJ + i] * x.r[ix - x.base];
            }
            // Y(JY) = Y(JY) + ALPHA*TEMP2
            y.r[jy - y.base] += AlphaRe * temp2Re - AlphaIm * temp2Im;
            y.i[jy - y.base] += AlphaRe * temp2Im + AlphaIm * temp2Re;

            jx += incx;
            jy += incy;
        }
    }
}



