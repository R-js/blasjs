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
    div_rxr,
    FortranArrEComplex,
    MatrixEComplex,
    mul_rxr
} from '../../../f_func';

export function normalUpper(
    kx: number,
    _noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number
) {

    let jx = kx + (n - 1) * incx;
    for (let j = n; j >= 1; j--) {
        const coorAJ = a.colOfEx(j);
        const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!xIsZero) {
            if (nounit) {
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[coorAJ + j],
                    a.i[coorAJ + j]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = jx;
            for (let i = j - 1; i >= 1; i--) {
                ix -= incx;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    a.r[coorAJ + i],
                    a.i[coorAJ + i]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
            }
        }
        jx -= incx;
    }
}
