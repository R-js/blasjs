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

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function ssymv(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number): void {

    const lu = lowerChar(uplo);

    //input check

    let info = 0;

    if (!'ul'.includes(lu)) {
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
        info = 9
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssymv', info));
    }

    // Quick return if possible.
    if (n === 0 || (alpha === 0 && beta === 1)) return;

    // Set up the start points in  X  and  Y.

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.

    //     First form  y := beta*y.


    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.
    //     First form  y := beta*y.

    if (beta !== 1) {

        let iy = ky;
        for (let i = 1; i <= n; i++) {
            y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
            iy += incy;
        }

    }

    if (alpha === 0) return;

    let jx = kx;
    let jy = ky;

    if (lu === 'u') {
        // Form  y  when A is stored in upper triangle.
        // First form  y := beta*y.  
        for (let j = 1; j <= n; j++) {
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            let ix = kx;
            let iy = ky;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= j - 1; i++) {
                y.r[iy - y.base] += temp1 * a.r[coords + i];
                temp2 += a.r[coords + i] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1 * a.r[coords + j] + alpha * temp2;
            jx += incx;
            jy += incy;
        }
    }
    else {
        // Form  y,  when A is stored in lower triangle.
        for (let j = 1; j <= n; j++) {
            //
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            const coords = a.colOfEx(j);
            y.r[jy - y.base] += temp1 * a.r[coords + j];
            let ix = jx;
            let iy = jy;
            //
            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                iy += incy;
                y.r[iy - y.base] += temp1 * a.r[coords + i];
                temp2 += a.r[coords + i] * x.r[ix - x.base];
            }
            //
            y.r[jy - y.base] += alpha * temp2;
            jx += incx;
            jy += incy;
            //
        }
    }
}
