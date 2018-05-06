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

import { FortranArrEComplex, mul_rxr } from '../../../f_func';

export function normUpper(
    kx: number,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {


    let kk = 1;
    let jx = kx;
    for (let j = 1; j <= n; j++) {
        const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!xIsZero) {
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = kx;
            for (let k = kk; k <= kk + j - 2; k++) {
                const apk = k - ap.base;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    ap.r[apk],
                    ap.i[apk]
                );
                x.r[ix - x.base] += re;
                x.i[ix - x.base] += im;
                ix += incx;
            }
            if (nounit) {
                const apk2 = kk + j - 1 - ap.base;
                const { re, im } = mul_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    ap.r[apk2],
                    ap.i[apk2]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
        }
        jx += incx;
        kk += j;
    }
}
