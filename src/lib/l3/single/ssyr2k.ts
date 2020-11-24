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

export function ssyr2k(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: number,
    c: Matrix,
    ldc: number,
): void {
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);

    const nrowA = tr === 'n' ? n : k;

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (!'ntc'.includes(tr)) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (k < 0) {
        info = 4;
    } else if (lda < max(1, nrowA)) {
        info = 7;
    } else if (ldb < max(1, nrowA)) {
        info = 9;
    } else if (ldc < max(1, n)) {
        info = 12;
    }

    if (info) {
        throw new Error(errWrongArg('ssyr2k', info));
    }

    //  Quick return if possible.

    if (n === 0 || ((alpha === 0 || k === 0) && beta === 1)) return;

    //*     And when  alpha.eq.zero.

    if (alpha === 0) {
        for (let j = 1; j <= n; j++) {
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            if (beta === 0) {
                c.setCol(j, start, stop, 0);
            } else {
                const coor = c.colOfEx(j);
                for (let i = start; i <= stop; i++) {
                    c.r[coor + i] *= beta;
                }
            }
        }
        return;
    }

    // Start the operations.
    if (tr === 'n') {
        // Form  C := alpha*A*B**T + alpha*B*A**T + C.
        if (ul === 'u') {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, 1, j, 0);
                } else if (beta !== 1) {
                    for (let i = 1; i <= j; i++) {
                        c.r[coorCJ + i] *= beta;
                    }
                }
                for (let l = 1; l <= k; l++) {
                    const coorAL = a.colOfEx(l);
                    const coorBL = b.colOfEx(l);
                    if (a.r[coorAL + j] !== 0 || b.r[coorBL + j] !== 0) {
                        const temp1 = alpha * b.r[coorBL + j];
                        const temp2 = alpha * a.r[coorAL + j];
                        for (let i = 1; i <= j; i++) {
                            c.r[coorCJ + i] += a.r[coorAL + i] * temp1 + b.r[coorBL + i] * temp2;
                        }
                    }
                }
            }
        }
        //uplu="l"
        else {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);

                if (beta === 0) {
                    c.setCol(j, j, n, 0);
                } else if (beta !== 1) {
                    for (let i = j; i <= n; i++) {
                        c.r[coorCJ + i] *= beta;
                    }
                }
                for (let L = 1; L <= k; L++) {
                    const coorBL = b.colOfEx(L);
                    const coorAL = a.colOfEx(L);
                    if (a.r[coorAL + j] !== 0 || b.r[coorBL + j] !== 0) {
                        const temp1 = alpha * b.r[coorBL + j];
                        const temp2 = alpha * a.r[coorAL + j];
                        for (let i = j; i <= n; i++) {
                            c.r[coorCJ + i] += a.r[coorAL + i] * temp1 + b.r[coorBL + i] * temp2;
                        }
                    }
                }
            }
        }
    } else {
        //Form  C := alpha*A**T*B + alpha*B**T*A + C.
        if (ul === 'u') {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= j; i++) {
                    let temp1 = 0;
                    let temp2 = 0;
                    const coorAI = a.colOfEx(i);
                    const coorBI = b.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        temp1 += a.r[coorAI + L] * b.r[coorBJ + L];
                        temp2 += b.r[coorBI + L] * b.r[coorBJ + L];
                    }
                    if (beta === 0) {
                        c.r[coorCJ + i] = alpha * temp1 + alpha * temp2;
                    } else {
                        c.r[coorCJ + i] = beta * c.r[coorCJ + i] + alpha * temp1 + alpha * temp2;
                    }
                }
            }
        } else {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorAJ = a.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = j; i <= n; i++) {
                    let temp1 = 0;
                    let temp2 = 0;
                    const coorAI = a.colOfEx(i);
                    const coorBI = b.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        temp1 += a.r[coorAI + L] * b.r[coorBJ + L];
                        temp2 += b.r[coorBI + L] * a.r[coorAJ + L];
                    }
                    c.r[coorCJ + i] = alpha * (temp1 + temp2) + (beta !== 0 ? beta * c.r[coorCJ + i] : 0);
                }
            }
        }
    }
}
