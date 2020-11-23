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

export function sgemv(
    trans: string,
    m: number,
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number,
): void {
    const tr = lowerChar(trans);

    let info = 0;
    if (!'ntc'.includes(tr)) {
        info = 1;
    } else if (m < 0) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (lda < max(1, m)) {
        info = 6;
    } else if (incx === 0) {
        info = 8;
    } else if (incy === 0) {
        info = 11;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('sgemv', info));
    }

    //Quick return if possible.

    if (m === 0 || n === 0 || (alpha === 0 && beta === 1)) return;

    const lenx = tr === 'n' ? n : m;
    const leny = tr === 'n' ? m : n;

    const kx = incx > 0 ? 1 : 1 - (lenx - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (leny - 1) * incy;

    /*
      Start the operations. In this version the elements of A are
      accessed sequentially with one pass through A.
      
      First form  y := beta*y.
    */

    if (beta !== 1) {
        //performance
        if (incy === 1 && beta === 0) {
            y.r.fill(0, 1 - y.base, 1 - y.base + leny);
        } else {
            let iy = ky;
            for (let i = 1; i <= leny; i++) {
                y.r[iy - y.base] = beta === 0 ? 0 : y.r[iy - y.base] * beta;
                iy += incy;
            }
        }
    }
    if (alpha === 0) return;
    if (tr === 'n') {
        //Form  y := alpha*A*x + y.
        let jx = kx;
        for (let j = 1; j <= n; j++) {
            const temp = alpha * x.r[jx - x.base];
            let iy = ky;
            const coorAJ = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                y.r[iy - y.base] += temp * a.r[coorAJ + i];
                iy += incy;
            }
            jx += incx;
        }
    } else {
        //Form  y:= alpha * A ** T * x + y.
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp = 0;
            let ix = kx;
            const coorAJ = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                temp += a.r[coorAJ + i] * x.r[ix - x.base];
                ix += incx;
            }
            y.r[jy - y.base] += alpha * temp;
            jy += incy;
        }
    }
}
