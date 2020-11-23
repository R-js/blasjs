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

import type {
    FortranArrEComplex,
    MatrixEComplex
} from '../../../f_func';

export function normalUpper(
    kx: number,
    _noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number
): void {

    let jx = kx - x.base;
    for (let j = 1; j <= n; j++) {
        const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
        if (!xIsZero) {
            let tempRe = x.r[jx];
            let tempIm = x.i[jx];
            let ix = kx - x.base;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= j - 1; i++) {
                x.r[ix] += tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                x.i[ix] += tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                ix += incx;
            }
            if (nounit) {
                const tr = x.r[jx] * a.r[coords + j] - x.i[jx] * a.i[coords + j];
                const ti = x.r[jx] * a.i[coords + j] + x.i[jx] * a.r[coords + j];
                x.r[jx] = tr;
                x.i[jx] = ti;
            }
        }
        jx += incx;
    }
}
