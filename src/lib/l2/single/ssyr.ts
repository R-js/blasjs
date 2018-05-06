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

export function ssyr(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    a: Matrix,
    lda: number
): void {

    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(uplo)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (lda < max(1, n)) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssyr', info));
    }

    // Quick return if possible

    if (n === 0 || alpha === 0) return;

    //Original had .LE.0 but the 0 case is already handled
    const kx = (incx < 0) ? 1 - (n - 1) * incx : 1;
    //let ky = (incy < 0) ? 1 - (n - 1) * incy : 1;

    let jx = kx;

    if (ul === 'u') {
        //Form  A  when A is stored in upper triangle. 
        for (let j = 1; j <= n; j++) {

            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = kx;
                const coords = a.colOfEx(j);
                for (let i = 1; i <= j; i++) {
                    a.r[coords + i] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
        }
    }
    else {
        //  Form  A  when A is stored in lower triangle.
        jx = kx;
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = jx;
                const coorAJ = a.colOfEx(j);
                let i = j
                for (; i <= n; i++) {
                    const delta = x.r[ix - x.base] * temp;
                    //   console.log(`i:${i},j:${j}, A_${i}.${j}=${a.r[coorAJ + i]}, x_${i}.x_${j}*alpha=${delta},A_N=${a.r[coorAJ + i] + delta}`);

                    a.r[coorAJ + i] = a.r[coorAJ + i] + delta;
                    // console.log(a.r[coorAJ + i]);
                    ix += incx;
                }
            }
            jx += incx;
        }
    }
}









