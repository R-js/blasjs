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
    mul_rxr
} from '../../../f_func';

export function normUpper(
    kx: number,
    _noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {
    // console.log({ kx, noconj, nounit, incx, n });

    let kk = n * (n + 1) / 2;
    let jx = kx + (n - 1) * incx;
    //console.log({ kk, jx });
    for (let j = n; j >= 1; j--) {
        const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!xIsZero) {
            //console.log({ j, jx, kk })
            if (nounit) {
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    ap.r[kk - ap.base],
                    ap.i[kk - ap.base]
                );
                //  console.log(`j:${j}\tx=(${x.r[jx - x.base]},${x.i[jx - x.base]})\tx/ap=(${re},${im}`);
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
                // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = jx;
            for (let k = kk - 1; k >= kk - j + 1; k--) {
                ix -= incx;
                const apk = k - ap.base;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    ap.r[apk],
                    ap.i[apk]
                );
                //console.log(`ix:${ix}\tx(ix):(${x.r[ix - x.base]},${x.i[ix - x.base]}\tte*ap(k):(${re},${im})`);
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
            }
        }
        jx -= incx;
        kk -= j;
    }
}
