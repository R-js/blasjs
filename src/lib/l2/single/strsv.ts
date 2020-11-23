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

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function strsv(
    uplo: 'u' | 'l',
    trans: 't' | 'c' | 'n',
    diag: 'u' | 'n',
    n: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
): void {
    // lowerCase it all in a fast way
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0;

    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (!'ntc'.includes(tr)) {
        info = 2;
    } else if (!'un'.includes(dg)) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (lda < max(1, n)) {
        info = 6;
    } else if (incx === 0) {
        info = 8;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('strsv', info));
    }

    if (n === 0) return;

    const nounit = dg === 'n';

    /*
    Set up the start point in X if the increment is not unity. This
         will be  ( N - 1 )*INCX  too small for descending loops.
    */

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        //Form  x := inv( A )*x.
        if (ul === 'u') {
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx] !== 0) {
                    const coords = a.colOfEx(j);

                    if (nounit) x.r[jx] /= a.r[coords + j];
                    const temp = x.r[jx];
                    let ix = jx;
                    for (let i = j - 1; i >= 1; i--) {
                        ix -= incx;
                        x.r[ix] -= temp * a.r[coords + i];
                    }
                }
                jx -= incx;
            }
        } else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const coords = a.colOfEx(j);
                if (x.r[jx] !== 0) {
                    if (nounit) x.r[jx] /= a.r[coords + j];
                    const temp = x.r[jx];
                    let ix = jx;
                    for (let i = j + 1; i <= n; i++) {
                        ix += incx;
                        x.r[ix] -= temp * a.r[coords + i];
                    }
                }
                jx += incx;
            }
        }
    } else {
        if (ul === 'u') {
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                let temp = x.r[jx - x.base];
                let ix = kx;
                const coords = a.colOfEx(j);
                for (let i = 1; i <= j - 1; i++) {
                    temp -= a.r[coords + i] * x.r[ix - x.base];
                    ix += incx;
                }
                if (nounit) temp /= a.r[coords + j];
                x.r[jx - x.base] = temp;
                jx += incx;
            }
        } else {
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx - x.base];
                let ix = kx;
                const coords = a.colOfEx(j);
                for (let i = n; i >= j + 1; i--) {
                    temp -= a.r[coords + i] * x.r[ix - x.base];
                    ix -= incx;
                }
                if (nounit) temp /= a.r[coords + j];
                x.r[jx - x.base] = temp;
                jx -= incx;
            }
        }
    }
}
