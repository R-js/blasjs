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

export function sspr2(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    ap: FortranArr,
): void {
    // validate input parameters

    let info = 0;
    const ul = lowerChar(uplo);
    //console.log('start', ul);

    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (n < 0) {
        info = 2;
    } else if (incx === 0) {
        info = 5;
    } else if (incy === 0) {
        info = 7;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('sspr2', info));
    }
    //console.log('startw2');
    //     Quick return if possible.

    if (n === 0 || alpha === 0) return;

    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    let jx = kx;
    let jy = ky;
    //console.log('startw3');
    //    Start the operations. In this version the elements of the array AP
    //    are accessed sequentially with one pass through AP.

    let kk = 1;
    if (ul === 'u') {
        //Form  A  when upper triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0 || y.r[jy - y.base] !== 0) {
                const temp1 = alpha * y.r[jy - y.base];
                //console.log({ temp1 });
                const temp2 = alpha * x.r[jx - x.base];
                //console.log({ temp2 });
                let ix = kx;
                let iy = ky;
                for (let k = kk; k <= kk + j - 1; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
            kk += j;
        }
    } else {
        // Form  A  when lower triangle is stored in AP.

        for (let j = 1; j <= n; j++) {
            //console.log('startw5');
            const xAndyIsZero = x.r[jx - x.base] === 0 && y.r[jy - y.base] === 0;
            //console.log(`xyIsZero=${xAndyIsZero}`);
            if (!xAndyIsZero) {
                //console.log('startw6', alpha, y.r);
                const temp1 = alpha * y.r[jy - y.base];
                const temp2 = alpha * x.r[jx - x.base];
                //console.log({ temp1 });
                let ix = jx;
                let iy = jy;
                for (let k = kk; k <= kk + n - j; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    //console.log({ temp1 });
                    ix += incx;
                    iy += incy;
                }
            }
            //console.log('startw7');
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}
