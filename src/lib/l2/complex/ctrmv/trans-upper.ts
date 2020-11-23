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

import type { FortranArrEComplex, MatrixEComplex } from '../../../f_func';

export function transUpper(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number,
): void {
    let jx = kx + (n - 1) * incx - x.base;
    for (let j = n; j >= 1; j--) {
        let tempRe = x.r[jx];
        let tempIm = x.i[jx];
        let ix = jx;
        const coords = a.colOfEx(j);
        if (noconj) {
            if (nounit) {
                const tr = tempRe * a.r[coords + j] - tempIm * a.i[coords + j];
                const ti = tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                tempRe = tr;
                tempIm = ti;
            }
            for (let i = j - 1; i >= 1; i--) {
                ix -= incx;
                tempRe += a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                tempIm += a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
            }
        } else {
            if (nounit) {
                //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                const tr = tempRe * a.r[coords + j] + tempIm * a.i[coords + j];
                const ti = -tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                tempRe = tr;
                tempIm = ti;
            }
            for (let i = j - 1; i >= 1; i--) {
                ix -= incx;
                //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                tempRe += a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                tempIm += a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
            }
        }
        x.r[jx] = tempRe;
        x.i[jx] = tempIm;
        jx -= incx;
    }
}
