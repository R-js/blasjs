import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

/*
 This is a conversion from BLAS to Typescript/Javascript
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

const { max, min } = Math;

export function stbsv(
    uplo: 'u' | 'l',
    trans: 't' | 'n' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    A: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

    const ul = lowerChar(uplo);
    const dg = lowerChar(diag);
    const tr = lowerChar(trans);

    let info = 0;

    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (tr !== 'n' && tr !== 't' && tr !== 'c') {
        info = 2;
    }
    else if (dg !== 'u' && dg !== 'n') {
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
        throw new Error(errWrongArg('stbsv', info));
    }

    //     Quick return if possible.

    if (n === 0) return;

    //     Set up the start point in X if the increment is not unity. This
    //     will be  ( N - 1 )*INCX  too small for descending loops.

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;

    //     Start the operations. In this version the elements of A are
    //    accessed by sequentially with one pass through A.


    if (tr === 'n') {
        //  Form  x := inv( A )*x.

        if (ul === 'u') {
            let kplus1 = k + 1;
            kx += (n - 1) * incx;
            let jx = kx;

            for (let j = n; j >= 1; j--) {
                kx = kx - incx;
                if (x.r[jx - x.base] !== 0) {
                    let ix = kx;
                    let L = kplus1 - j;
                    const coords = A.colOfEx(j);
                    if (dg === 'n') x.r[jx - x.base] /= A.r[kplus1 + coords];
                    let temp = x.r[jx - x.base];
                    for (let i = j - 1; i >= max(1, j - k); i--) {
                        x.r[ix - x.base] -= temp * A.r[coords + L + i];
                        ix -= incx;
                    }
                }
                jx -= incx;
            }
        }
        // ul === 'L'
        else {
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                kx += incx;
                if (x.r[jx - x.base] !== 0) {
                    let ix = kx;
                    let L = 1 - j;
                    const coords = A.colOfEx(j);
                    if (dg === 'n') x.r[jx - x.base] /= A.r[1 + coords];
                    const temp = x.r[jx - x.base];
                    for (let i = j + 1; i <= min(n, j + k); i++) {
                        x.r[ix - x.base] -= temp * A.r[coords + L + i];
                        ix += incx;
                    }
                }
                jx += incx;
            }
        }
    }
    // tr !== 'N', aka tr in[']
    else {
        // Form  x := inv( A**T)*x.
        if (ul === 'u') {
            let kplus1 = k + 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {

                let temp = x.r[jx - x.base];
                let ix = kx;
                let L = kplus1 - j;
                const coords = A.colOfEx(j);
                for (let i = max(1, j - k); i <= j - 1; i++) {
                    temp -= A.r[coords + L + i] * x.r[ix - x.base];
                    ix += incx;
                }
                if (dg === 'n') temp /= A.r[coords + kplus1];
                x.r[jx - x.base] = temp;
                jx += incx;
                if (j > k) kx += incx;
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx - x.base];
                let ix = kx;
                let L = 1 - j;
                const coords = A.colOfEx(j);
                for (let i = min(n, j + k); i >= j + 1; i--) {
                    temp -= A.r[L + i + coords] * x.r[ix - x.base];
                    ix -= incx;
                }
                if (dg === 'n') temp /= A.r[coords + 1];
                x.r[jx - x.base] = temp;
                jx -= incx;
                if ((n - j) >= k) kx -= incx;
            }
        }
    }
}
