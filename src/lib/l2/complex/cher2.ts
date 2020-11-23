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

import { Complex, errMissingIm, errWrongArg, FortranArr, isZero, lowerChar, Matrix } from '../../f_func';

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
    lda: number,
): void {
    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (incx === 0) {
        info = 5;
    } else if (incy === 0) {
        info = 7;
    } else if (lda < max(1, n)) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('cher2', info));
    }

    const { re: AlphaRe, im: AlphaIm } = alpha;

    const alphaIsZero = isZero(alpha);

    if (n === 0 || alphaIsZero) return; //nothing to do

    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;
    let jx = kx;
    let jy = ky;

    if (ul === 'u') {
        // Form  A  when A is stored in the upper triangle.
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
            const yIsZero = y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0;

            if (!xIsZero || !yIsZero) {
                //console.log(`${jx},${jy}  : ${xIsZero},${yIsZero}`);
                // TEMP1 = ALPHA*CONJG(Y(JY))
                const temp1Re = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
                const temp1Im = -AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];

                // TEMP2 = CONJG(ALPHA*X(JX))
                const temp2Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
                const temp2Im = -(AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base]);

                //console.log(`${jx},${jy}, (${temp1Re},${temp1Im}),(${temp2Re},${temp2Im})`);

                let ix = kx;
                let iy = ky;

                for (let i = 1; i <= j - 1; i++) {
                    // X(IX)*TEMP1
                    const re1 = x.r[ix - x.base] * temp1Re - x.i[ix - x.base] * temp1Im;
                    const im1 = x.r[ix - x.base] * temp1Im + x.i[ix - x.base] * temp1Re;
                    //Y(IY)*TEMP2
                    const re2 = y.r[iy - y.base] * temp2Re - y.i[iy - y.base] * temp2Im;
                    const im2 = y.r[iy - y.base] * temp2Im + y.i[iy - y.base] * temp2Re;

                    //console.log(`${jx},${jy}, (${re1},${im1}),(${re2},${im2})`);
                    a.r[coorAJ + i] += re1 + re2;
                    a.i[coorAJ + i] += im1 + im2;
                    //console.log(`i:${i},j:${j}, (${re1 + re2},${im1 + im2})`);
                    ix += incx;
                    iy += incy;
                } //for
                //X(JX)*TEMP1
                const re1 = x.r[jx - x.base] * temp1Re - x.i[jx - x.base] * temp1Im;
                //Y(JY)*TEMP2
                const re2 = y.r[jy - y.base] * temp2Re - y.i[jy - y.base] * temp2Im;

                a.r[coorAJ + j] += re1 + re2;
                a.i[coorAJ + j] = 0;
            } else {
                a.i[coorAJ + j] = 0;
            }
            jx += incx;
            jy += incy;
        } //for
    } else {
        //  Form  A  when A is stored in the lower triangle.

        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
            const yIsZero = y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0;

            if (!xIsZero || !yIsZero) {
                //TEMP1 = ALPHA*CONJG(Y(JY))
                // (a+ib)*(c-di)=  (ac+bd) + i(-ad+bc)
                const temp1Re = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
                const temp1Im = -AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];

                //TEMP2 = DCONJG(ALPHA * X(JX))
                const temp2Re = AlphaRe * x.r[jx - y.base] - AlphaIm * x.i[jx - x.base];
                const temp2Im = -(AlphaRe * x.i[jx - y.base] + AlphaIm * x.r[jx - x.base]);

                a.i[coorAJ + j] = 0;
                a.r[coorAJ + j] +=
                    // REAL( X(JX)*TEMP1 )
                    x.r[jx - x.base] * temp1Re -
                    x.i[jx - x.base] * temp1Im +
                    // REAL(  Y(JY)*TEMP2  )
                    (y.r[jy - y.base] * temp2Re - y.i[jy - y.base] * temp2Im);
                let ix = jx;
                let iy = jy;
                for (let i = j + 1; i <= n; i++) {
                    ix += incx;
                    iy += incy;
                    //X(IX)*TEMP1
                    const re1 = x.r[ix - x.base] * temp1Re - x.i[ix - x.base] * temp1Im;
                    const im1 = x.r[ix - x.base] * temp1Im + x.i[ix - x.base] * temp1Re;

                    //Y(IY)*TEMP2
                    const re2 = y.r[iy - y.base] * temp2Re - y.i[iy - y.base] * temp2Im;
                    const im2 = y.r[iy - y.base] * temp2Im + y.i[iy - y.base] * temp2Re;

                    a.r[coorAJ + i] += re1 + re2;
                    a.i[coorAJ + i] += im1 + im2;
                }
            } else {
                a.i[coorAJ + j] = 0;
            }
            jx += incx;
            jy += incy;
        }
    }
}
