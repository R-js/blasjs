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


import { errWrongArg, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function ssyrk(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    beta: number,
    c: Matrix,
    ldc: number): void {

    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);

    const nrowA = (tr === 'n') ? n : k;

    //Test the input parameters.

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (k < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (ldc < max(1, n)) {
        info = 10;
    }

    if (info) {
        throw new Error(errWrongArg('ssyrk', info));
    }

    //*     Quick return if possible.

    if (n === 0 || ((alpha === 0 || k === 0) && beta === 1)) return;

    //*      And when  alpha.eq.zero.
    if (alpha === 0) {
        // if (beta !== 1) {
        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = beta === 0 ? 0 : beta * c.r[coorCJ + i];
            }
        }
        // }
        return;
    }

    //*     Start the operations.

    if (tr === 'n') {
        //Form  C := alpha*A*A**T + beta*C.

        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            if (beta === 0) {
                //zap it
                c.r.fill(0, coorCJ + start, coorCJ + stop + 1);
            }
            else if (beta !== 1) {
                for (let i = start; i <= stop; i++) {
                    c.r[coorCJ + i] *= beta;
                }
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                if (a.r[coorAL + j] !== 0) {
                    let temp = alpha * a.r[coorAL + j];
                    for (let i = start; i <= stop; i++) {
                        c.r[coorCJ + i] += temp * a.r[coorAL + i];
                    }
                }
            }
        }
    }
    else {
        // Form  C := alpha*A**T*A + beta*C.
        for (let j = 1; j <= n; j++) {
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            // console.log(`${j},${start},${stop}`);
            for (let i = start; i <= stop; i++) {
                let temp = 0;
                const coorAI = a.colOfEx(i);
                for (let l = 1; l <= k; l++) {
                    temp += a.r[coorAI + l] * a.r[coorAJ + l];
                }
                //console.log(`${j},${i}\t${temp}`);

                let re = alpha * temp;
                if (beta !== 0) {
                    re += beta * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
            }
        }
    }
}
