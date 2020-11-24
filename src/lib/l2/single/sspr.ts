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

import { errWrongArg, FortranArr, lowerChar } from '../../f_func';

export function sspr(uplo: 'u' | 'l', n: number, alpha: number, x: FortranArr, incx: number, ap: FortranArr): void {
    // test parameters

    let info = 0;
    const ul = lowerChar(uplo);

    if (!'ul'.includes(uplo)) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (incx === 0) {
        info = 5;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('sspr', info));
    }

    if (n === 0 || alpha === 0) return;

    // corrected the comparison was incx <= 0
    const kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    let kk = 1;
    let jx = kx;

    if (ul === 'u') {
        //  Form  A  when upper triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                const temp = alpha * x.r[jx - x.base];
                let ix = kx;
                for (let k = kk; k <= kk + j - 1; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
            kk += j;
        }
    } else {
        //  Form  A  when lower triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                const temp = alpha * x.r[jx - x.base];
                let ix = jx;
                for (let k = kk; k <= kk + n - j; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
            kk += n - j + 1;
        }
    }
}
