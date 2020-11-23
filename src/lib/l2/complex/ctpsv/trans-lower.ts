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

import { div_rxr, FortranArrEComplex, mul_rxr } from '../../../f_func';

export function transLower(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number,
): void {
    let kk = (n * (n + 1)) / 2;
    kx += (n - 1) * incx;
    let jx = kx;
    for (let j = n; j >= 1; j--) {
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        let ix = kx;

        for (let k = kk; k >= kk - (n - (j + 1)); k--) {
            const apk = k - ap.base;
            const { re, im } = mul_rxr(ap.r[apk], noconj ? ap.i[apk] : -ap.i[apk], x.r[ix - x.base], x.i[ix - x.base]);
            //console.log(`k:${k},apk=(${ap.r[apk]},${(noconj ? ap.i[apk] : -ap.i[apk])})`);
            tempRe -= re;
            tempIm -= im;
            ix -= incx;
        }
        if (nounit) {
            // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
            const apkk = kk - n + j - ap.base;
            const { re, im } = div_rxr(tempRe, tempIm, ap.r[apkk], noconj ? ap.i[apkk] : -ap.i[apkk]);
            tempRe = re;
            tempIm = im;
            //console.log(`${kk - n + j},\t(${re},${im})`)
        }
        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx -= incx;
        kk -= n - j + 1;
    }
}
