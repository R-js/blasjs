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

import { div_rxr, FortranArrEComplex, MatrixEComplex, mul_rxr } from '../../../f_func';
const { max } = Math;

export function normUpper(
    kx: number,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    _noconj: boolean,
    nounit: boolean,
    n: number,
    k: number,
): void {
    const kplus1 = k + 1;
    kx += (n - 1) * incx;
    let jx = kx;

    for (let j = n; j >= 1; j--) {
        const coorAJ = a.colOfEx(j);
        kx -= incx;
        const isXZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        const extrI = max(1, j - k);
        if (!isXZero) {
            let ix = kx;
            const L = kplus1 - j;
            if (nounit) {
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[kplus1 + coorAJ],
                    a.i[kplus1 + coorAJ],
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            const tempRe = x.r[jx - x.base];
            const tempIm = x.i[jx - x.base];

            for (let i = j - 1; i >= extrI; i--) {
                const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAJ + i + L], a.i[coorAJ + i + L]);
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
                ix -= incx;
            }
        }
        jx -= incx;
    }
}
