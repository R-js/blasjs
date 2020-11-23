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

import { div_rxr, FortranArrEComplex, isZeroE, MatrixEComplex, mul_rxr } from '../../../f_func';

const { min } = Math;

export function normLower(
    kx: number,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    _noconj: boolean,
    nounit: boolean,
    n: number,
    k: number) {

    /*if (!x.i || !a.i) {
        return;
    }*/

    let jx = kx;
    for (let j = 1; j <= n; j++) {
        kx += incx;
        const xIsZero = isZeroE(x.r[jx - x.base], x.i[jx - x.base]);
        const extrI = min(n, j + k);
        if (!xIsZero) {
            let ix = kx;
            let L = 1 - j;
            const coorAJ = a.colOfEx(j);
            if (nounit) {
                // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[coorAJ + 1],
                    a.i[coorAJ + 1]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            for (let i = j + 1; i <= extrI; i++) {
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    a.r[coorAJ + L + i],
                    a.i[coorAJ + L + i]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
                ix += incx;
            }
        }
        jx += incx;
    }
}
