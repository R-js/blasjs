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

import { errWrongArg, FortranArr, Matrix } from '../../f_func';

const { max } = Math;

export function sger(
    m: number,
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    a: Matrix,
    lda: number,
): void {
    let err = 0;
    switch (true) {
        case m < 0:
            err = 1;
            break;
        case n < 0:
            err = 2;
            break;
        case incx === 0:
            err = 5;
            break;
        case incy === 0:
            err = 7;
            break;
        case lda < max(1, m):
            err = 9;
            break;
        default:
            err = 0;
    }

    if (err) {
        throw new Error(errWrongArg('sger', err));
    }

    if (m === 0 || n === 0 || alpha === 0) return;

    let jy = incy < 0 ? 1 - (n - 1) * incy : 1;
    const kx = incx < 0 ? 1 - (m - 1) * incx : 1;
    for (let j = 1; j <= n; j++) {
        if (y.r[jy - y.base] !== 0) {
            const temp = alpha * y.r[jy - y.base];
            let ix = kx;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                a.r[coords + i] += x.r[ix - x.base] * temp;
                ix += incx;
            }
        }
        jy += incy;
    }
}
