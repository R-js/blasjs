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

const { max, min } = Math;

export function stbmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    A: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0;

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
    else if (k < 0) {
        info = 5;
    }
    else if (lda < (k + 1)) {
        info = 7;
    }
    else if (incx === 0) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('stbmv', info));
    }

    if (n === 0) return;

    const nounit = diag === 'n';


    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;

    if (tr === 'n') {
        // Form  x := A*x.
        if (ul === 'u') {
            const kplus1 = k + 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    let L = kplus1 - j;
                    let coords = A.colOfEx(j);
                    for (let i = max(1, j - k); i <= j - 1; i++) {
                        x.r[ix - x.base] += temp * A.r[coords + L + i];
                        ix += incx;
                    }
                    if (nounit) x.r[jx - x.base] *= A.r[coords + kplus1];
                }
                jx += incx;
                if (j > k) kx += incx;
            }
        }
        else { //lower-matrix
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    const l = 1 - j;
                    const coorAJ = A.colOfEx(j);
                    for (let i = min(n, j + k); i >= j + 1; i--) {
                        x.r[ix - x.base] += temp * A.r[coorAJ + l + i];
                        ix -= incx;
                    }
                    if (nounit) x.r[jx - x.base] *= A.r[coorAJ + 1];

                }
                jx -= incx;
                if ((n - j) >= k) kx -= incx;
            }
        }
    }//N
    else {
        // Form  x := A**T*x.
        if (uplo === 'u') {
            const kplus1 = k + 1;
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                const coorAJ = A.colOfEx(j);
                let temp = x.r[jx - x.base];
                kx -= incx;
                let ix = kx;
                let l = kplus1 - j;
                if (nounit) temp *= A.r[kplus1 + coorAJ]
                for (let i = j - 1; i >= max(1, j - k); i--) {
                    temp += A.r[l + i + coorAJ] * x.r[ix - x.base];
                    ix -= incx;
                }
                x.r[jx - x.base] = temp;
                jx -= incx;

            }

        }
        //lower matrix
        else {
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                const coorAJ = A.colOfEx(j);
                let temp = x.r[jx - x.base];
                kx += incx;
                let ix = kx;
                const l = 1 - j;
                if (nounit) temp *= A.r[coorAJ + 1];
                for (let i = j + 1; i <= min(n, j + k); i++) {
                    temp += A.r[coorAJ + l + i] * x.r[ix - x.base];
                    ix += incx;
                }
                x.r[jx - x.base] = temp;
                jx += incx;
            }
        }
    }
}
