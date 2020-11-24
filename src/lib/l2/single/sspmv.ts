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

export function sspmv(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    ap: FortranArr, // a symmetric matrix in packed form
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number,
): void {
    const ul = lowerChar(uplo);

    let info = 0;

    if (!'ul'.includes(uplo)) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (incx === 0) {
        info = 6;
    } else if (incy === 0) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('SSPMV', info));
    }

    if (n === 0 || (alpha === 0 && beta === 1)) return;

    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    //    First form  y := beta*y.

    if (beta !== 1) {
        let iy = ky;
        for (let i = 1; i <= n; i++) {
            y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
            iy += incy;
        }
    }
    //
    if (alpha === 0) return;

    let kk = 1;

    if (ul === 'u') {
        // Form  y  when AP contains the upper triangle.

        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            const temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            let ix = kx;
            let iy = ky;
            for (let k = kk; k <= kk + j - 2; k++) {
                y.r[iy - y.base] += temp1 * ap.r[k - ap.base];
                temp2 += ap.r[k - ap.base] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1 * ap.r[kk + j - 1 - ap.base] + alpha * temp2;
            jx += incx;
            jy += incy;
            kk += j;
        }
    } else {
        // Form  y  when AP contains the lower triangle.

        let jx = kx;
        let jy = ky;

        for (let j = 1; j <= n; j++) {
            const temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            y.r[jy - y.base] += temp1 * ap.r[kk - ap.base];
            let ix = jx;
            let iy = jy;

            for (let k = kk + 1; k <= kk + n - j; k++) {
                ix += incx;
                iy += incy;
                y.r[iy - y.base] += temp1 * ap.r[k - ap.base];
                temp2 += ap.r[k - ap.base] * x.r[ix - x.base];
            }
            y.r[jy - y.base] += alpha * temp2;
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}
