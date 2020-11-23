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

import type { Matrix } from '../../f_func';
import { errWrongArg, lowerChar } from '../../f_func';

const { max } = Math;

export function strsm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    m: number,
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number): void {

    const si = lowerChar(side);
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const di = lowerChar(diag);

    const lside = si === 'l';
    const nrowA = lside ? m : n;
    const nounit = di === 'n';
    const upper = ul === 'u';


    let info = 0;
    if (!'lr'.includes(si)) {
        info = 1;
    }
    else if (!'ul'.includes(ul)) {
        info = 2;
    }
    else if (!'ntc'.includes(tr)) {
        info = 3;
    }
    else if (!'un'.includes(di)) {
        info = 4;
    }
    else if (m < 0) {
        info = 5;
    }
    else if (n < 0) {
        info = 6;
    }
    else if (lda < max(1, nrowA)) {
        info = 9;
    }
    else if (ldb < max(1, m)) {
        info = 11;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('strsm', info));
    }

    // Quick return if possible.

    if (m === 0 || n === 0) return;

    if (alpha === 0) {
        for (let j = 1; j <= n; j++) {
            b.setCol(j, 1, m, 0);
        }
        return;
    }

    //*     Start the operations.

    if (lside) {
        if (tr === 'n') {
            // Form  B := alpha*inv( A )*B.
            if (upper) {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] *= alpha;
                        }
                    }
                    for (let k = m; k >= 1; k--) {
                        const coorAK = a.colOfEx(k);
                        if (b.r[coorBJ + k] !== 0) {
                            if (nounit) b.r[coorBJ + k] /= a.r[coorAK + k];
                            for (let i = 1; i <= k - 1; i++) {
                                b.r[coorBJ + i] -= b.r[coorBJ + k] * a.r[coorAK + i];
                            }
                        }
                    }
                }
                60
            }
            // not upper
            else {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] *= alpha;
                            //       console.log(`${j},${i}\t${b.r[coorBJ + i]}`)
                        }
                    }
                    for (let k = 1; k <= m; k++) {
                        //const coorBK = b.colOfEx(k);
                        const coorAK = a.colOfEx(k);
                        if (b.r[coorBJ + k] !== 0) {
                            if (nounit) {
                                b.r[coorBJ + k] /= a.r[coorAK + k];
                            }
                            //console.log(`${j},${k}\t${b.r[coorBJ + k]}`)
                            for (let i = k + 1; i <= m; i++) {
                                b.r[coorBJ + i] -= b.r[coorBJ + k] * a.r[coorAK + i];
                            }
                        }
                    }

                }
            }
        }
        else {
            //Form  B := alpha*inv( A**T )*B.
            if (upper) {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let i = 1; i <= m; i++) {
                        const coorAI = a.colOfEx(i);
                        let temp = alpha * b.r[coorBJ + i];
                        for (let k = 1; k <= i - 1; k++) {
                            temp -= a.r[coorAI + k] * b.r[coorBJ + k];
                        }
                        if (nounit) temp /= a.r[coorAI + i];
                        b.r[coorBJ + i] = temp;
                    }
                }
            }
            else {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let i = m; i >= 1; i--) {
                        const coorAI = a.colOfEx(i);
                        let temp = alpha * b.r[coorBJ + i];
                        for (let k = i + 1; k <= m; k++) {
                            temp -= a.r[coorAI + k] * b.r[coorBJ + k];
                        }
                        if (nounit) temp /= a.r[coorAI + i];
                        b.r[coorBJ + i] = temp;
                    }
                }
            }
        }
    }
    // if (lside)
    else {
        if (tr === 'n') {
            // Form  B := alpha*B*inv( A ).
            if (upper) {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    const coorAJ = a.colOfEx(j);
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] = alpha * b.r[coorBJ + i];
                        }
                    }
                    for (let k = 1; k <= j - 1; k++) {
                        const coorBK = b.colOfEx(k);
                        if (a.r[coorAJ + k] !== 0) {
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] -= a.r[coorAJ + k] * b.r[coorBK + i];
                            }
                        }
                    }
                    if (nounit) {
                        let temp = 1 / a.r[coorAJ + j];
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] *= temp;
                        }
                    }
                }//for
            }//if upper
            else {
                for (let j = n; j >= 1; j--) {
                    const coorBJ = b.colOfEx(j);
                    const coorAJ = a.colOfEx(j);
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] *= alpha;
                        }
                    }
                    for (let k = j + 1; k <= n; k++) {
                        const coorBK = b.colOfEx(k);
                        if (a.r[coorAJ + k] !== 0) {
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] -= a.r[coorAJ + k] * b.r[coorBK + i];
                            }
                        }
                    }
                    if (nounit) {
                        let temp = 1 / a.r[coorAJ + j];
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBJ + i] *= temp;
                        }
                    }
                }//for
            }//upper
        }
        else {
            // Form  B := alpha*B*inv( A**T ).
            if (upper) {
                for (let k = n; k >= 1; k--) {
                    const coorAK = a.colOfEx(k);
                    const coorBK = b.colOfEx(k);
                    if (nounit) {
                        let temp = 1 / a.r[coorAK + k];
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] *= temp;
                        }
                    }
                    for (let j = 1; j <= k - 1; j++) {
                        const coorBJ = b.colOfEx(j);
                        if (a.r[coorAK + j] !== 0) {
                            let temp = a.r[coorAK + j];
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] -= temp * b.r[coorBK + i];
                            }
                        }
                    }
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] *= alpha;
                        }
                    }
                }
            }
            else {
                for (let k = 1; k <= n; k++) {
                    const coorAK = a.colOfEx(k);
                    const coorBK = b.colOfEx(k);
                    if (nounit) {
                        let temp = 1 / a.r[coorAK + k];

                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] *= temp;
                        }
                    }
                    for (let j = k + 1; j <= n; j++) {
                        const coorBJ = b.colOfEx(j);
                        if (a.r[coorAK + j] !== 0) {
                            let temp = a.r[coorAK + j];
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] -= temp * b.r[coorBK + i];
                            }
                        }
                    }
                    if (alpha !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] *= alpha;
                        }
                    }
                }
            }
        }
    }
}
