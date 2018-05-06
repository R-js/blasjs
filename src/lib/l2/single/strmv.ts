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

export function strmv(
    uplo: 'u' | 'l',
    trans: 't' | 'c' | 'n',
    diag: 'u' | 'n',
    n: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0

    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (!'un'.includes(dg)) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (lda < max(1, n)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('strmv', info));
    }

    if (n === 0) return;

    let nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        // Form  x := A*x.
        if (ul === 'u') {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                //console.log(`jx:${jx}, x:${x.r[jx]}`);
                if (x.r[jx] !== 0) {
                    let temp = x.r[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = 1; i <= j - 1; i++) {
                        x.r[ix] += temp * a.r[coords + i];
                        ix = ix + incx;
                    }
                    if (nounit) x.r[jx] *= a.r[coords + j];
                }
                jx += incx;
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx] !== 0) {
                    let temp = x.r[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = n; i >= j + 1; i--) {
                        x.r[ix] += temp * a.r[coords + i];
                        ix -= incx;
                    }
                    if (nounit) x.r[jx] *= a.r[coords + j];
                }
                jx -= incx;
            }
        }
    }
    else {
        //  Form  x := A**T*x.
        if (ul === 'u') {
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (nounit) temp *= a.r[coords + j];
                for (let i = j - 1; i >= 1; i--) {
                    ix -= incx;
                    temp += a.r[coords + i] * x.r[ix];
                }
                x.r[jx] = temp;
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let temp = x.r[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (nounit) temp *= a.r[coords + j];
                for (let i = j + 1; i <= n; i++) {
                    ix += incx;
                    temp += a.r[coords + i] * x.r[ix];
                }
                x.r[jx] = temp;
                jx += incx;
            }
        }
    }
}
