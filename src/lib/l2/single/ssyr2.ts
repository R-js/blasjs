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

const { max } = Math;

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

export function ssyr2(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    A: Matrix,
    lda: number): void {

    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(uplo)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 3;
    }
    else if (incy === 0) {
        info = 5;
    }
    else if (lda < max(1, n)) { //n can be 0?
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('ssyr2', info));
    }

    // Quick return if possible.

    if (n === 0 || alpha === 0) return;

    //     Set up the start points in X and Y if the increments are not both
    //     unity
    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;
    //
    let jx = kx;
    let jy = ky;

    // Start the operations. In this version the elements of A are
    //    accessed sequentially with one pass through the triangular part
    //    of A.

    if (ul === 'u') {
        //Form  A  when A is stored in the upper triangle.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0 || y.r[jy - y.base] !== 0) {
                let temp1 = alpha * y.r[jy - y.base];
                let temp2 = alpha * x.r[jx - x.base];
                let ix = kx;
                let iy = ky;
                const coords = A.colOfEx(j);
                for (let i = 1; i <= j; i++) {
                    A.r[coords + i] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
        }//for
    }
    else {
        //Form  A  when A is stored in the lower triangle.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0 || y.r[jy - y.base] !== 0) {
                let temp1 = alpha * y.r[jy - y.base];
                let temp2 = alpha * x.r[jx - x.base];
                let ix = jx;
                let iy = jy;
                const coords = A.colOfEx(j);
                for (let i = j; i <= n; i++) {
                    A.r[coords + i] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
        }
    }
}
