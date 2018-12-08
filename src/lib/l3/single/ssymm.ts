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

export function ssymm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    m: number,
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: number,
    c: Matrix,
    ldc: number): void {


    const si = lowerChar(side);
    const ul = lowerChar(uplo);

    const nrowA = si === 'l' ? m : n;

    //Test the input parameters.

    let info = 0;
    if (si !== 'l' && si !== 'r') {
        info = 1;
    }
    else if (ul !== 'u' && ul !== 'l') {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (ldb < max(1, m)) {
        info = 9;
    }
    else if (ldc < max(1, m)) {
        info = 12;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssymm', info));
    }

    //*     Quick return if possible.

    if (m === 0 || n === 0 ||
        (alpha === 0 && beta === 1)
    ) return;

    //when alpha is zero
    if (alpha === 0) {
        if (beta === 0) {
            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                const coords = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    c.r[coords + i] *= beta;
                }
            }
        }
        return;
    }

    //   Start the operations.
    if (si === 'l') {
        //Form  C := alpha*A*B + beta*C.
        if (ul === 'u') {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    const coorAI = a.colOfEx(i);
                    let temp1 = alpha * b.r[coorBJ + i];
                    let temp2 = 0;
                    for (let k = 1; k <= i - 1; k++) {
                        c.r[coorCJ + k] += temp1 * a.r[coorAI + k];
                        temp2 += b.r[coorBJ + k] * a.r[coorAI + k];
                    }
                    let re = temp1 * a.r[coorAI + i] + alpha * temp2;
                    if (beta !== 0) {
                        re += beta * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                }
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = b.colOfEx(j);
                for (let i = m; i >= 1; i--) {
                    const coorAI = a.colOfEx(i);
                    const temp1 = alpha * b.r[coorBJ + i];
                    let temp2 = 0;
                    for (let k = i + 1; k <= m; k++) {
                        c.r[coorCJ + k] += temp1 * a.r[coorAI + k];
                        temp2 += b.r[coorBJ + k] * a.r[coorAI + k];
                    }
                    let re = temp1 * a.r[coorAI + i] + alpha * temp2;
                    if (beta !== 0) {
                        re += beta * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                }
            }
        }
    }
    else {
        //  Form  C := alpha*B*A + beta*C.
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            const coorBJ = b.colOfEx(j);

            const temp1 = alpha * a.r[coorAJ + j];
            if (beta === 0) {
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] = temp1 * b.r[coorBJ + i];
                }
            }
            else {
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] = beta * c.r[coorCJ + i] + temp1 * b.r[coorBJ + i];
                }
            }
            for (let k = 1; k <= j - 1; k++) {
                const coorAK = a.colOfEx(k);
                const coorBK = b.colOfEx(k);
                let temp1 = alpha * (ul === 'u' ? a.r[coorAJ + k] : a.r[coorAK + j]);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1 * b.r[coorBK + i];
                }
            }
            for (let k = j + 1; k <= n; k++) {
                const coorAK = a.colOfEx(k);
                const coorBK = b.colOfEx(k);

                let temp1 = alpha * (ul === 'u' ? a.r[coorAK + j] : a.r[coorAJ + k]);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1 * b.r[coorBK + i];
                }
            }
        }
    }
}
