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

export function strmm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    transA: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    m: number,
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number
): void {

    const si = lowerChar(side); //String.fromCharCode(side.charCodeAt(0) | 0X20);
    const ul = lowerChar(uplo); //String.fromCharCode(uplo.charCodeAt(0) | 0X20);
    const tr = lowerChar(transA); //(String.fromCharCode(transA.charCodeAt(0) | 0X20);
    const di = lowerChar(diag); //String.fromCharCode(diag.charCodeAt(0) | 0X20);

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
        throw new Error(errWrongArg('strmm', info));
    }

    // Quick return if possible.

    if (m === 0 || n === 0) return;

    // And when  alpha.eq.zero.

    if (alpha === 0) {
        for (let j = 1; j <= n; j++) {
            b.setCol(j, 1, m, 0);
        }
        return;
    }

    //Start the operations.

    if (lside) {
        if (tr === 'n') {
            //Form  B := alpha*A*B.
            if (upper) {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let k = 1; k <= m; k++) {
                        const coorAK = a.colOfEx(k);
                        if (b.r[coorBJ + k] !== 0) {
                            //let temp = alpha * b.r[coorBJ + k]; JKF
                            let temp = b.r[coorBJ + k];
                            for (let i = 1; i <= k - 1; i++) {
                                b.r[coorBJ + i] += temp * a.r[coorAK + i]
                            }
                            if (nounit) temp *= a.r[coorAK + k];
                            b.r[coorBJ + k] = temp;
                        }
                    }
                }
            }
            else {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let k = m; k >= 1; k--) {
                        const coorAK = a.colOfEx(k);
                        if (b.r[coorBJ + k] !== 0) {
                            let temp = alpha * b.r[coorBJ + k];
                            b.r[coorBJ + k] = temp;
                            if (nounit) {
                                b.r[coorBJ + k] *= a.r[coorAK + k]
                            }
                            for (let i = k + 1; i <= m; i++) {
                                b.r[coorBJ + i] += temp * a.r[coorAK + i]
                            }
                        }
                    }
                }
            }
        }
        else {
            //  Form  B := alpha*A**T*B.
            if (upper) {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let i = m; i >= 1; i--) {
                        let temp = b.r[coorBJ + i];
                        const coorAI = a.colOfEx(i);
                        if (nounit) temp *= a.r[coorAI + i];
                        for (let k = 1; k <= i - 1; k++) {
                            temp += a.r[coorAI + k] * b.r[coorBJ + k];
                        }
                        b.r[coorBJ + i] = alpha * temp;
                    }
                }
            }
            else {
                for (let j = 1; j <= n; j++) {
                    const coorBJ = b.colOfEx(j);
                    for (let i = 1; i <= m; i++) {
                        const coorAI = a.colOfEx(i);
                        let temp = b.r[coorBJ + i];
                        if (nounit) temp *= a.r[coorAI + i];
                        for (let k = i + 1; k <= m; k++) {
                            temp += a.r[coorAI + k] * b.r[coorBJ + k];
                        }
                        b.r[coorBJ + i] = alpha * temp;
                    }
                }
            }
        }
    }
    else {
        if (tr === 'n') {
            // Form  B := alpha*B*A.
            if (upper) {
                for (let j = n; j >= 1; j--) {
                    let temp = alpha;
                    const coorAJ = a.colOfEx(j);
                    const coorBJ = b.colOfEx(j);
                    if (nounit) {
                        temp *= a.r[coorAJ + j];
                    }
                    for (let i = 1; i <= m; i++) {
                        b.r[coorBJ + i] *= temp;
                    }
                    for (let k = 1; k <= j - 1; k++) {
                        const coorBK = b.colOfEx(k);
                        if (a.r[coorAJ + k] !== 0) {
                            temp = alpha * a.r[coorAJ + k];
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] += temp * b.r[coorBK + i];
                            }
                        }
                    }
                }
            }
            //lower
            else {
                for (let j = 1; j <= n; j++) {
                    const coorAJ = a.colOfEx(j);
                    const coorBJ = b.colOfEx(j);
                    let temp = alpha;
                    if (nounit) temp *= a.r[coorAJ + j];
                    //console.log(`${j}\t ${temp}`)
                    for (let i = 1; i <= m; i++) {
                        b.r[coorBJ + i] *= temp;
                    }
                    for (let k = j + 1; k <= n; k++) {
                        const coorBK = b.colOfEx(k);
                        if (a.r[coorAJ + k] !== 0) {
                            temp = alpha * a.r[coorAJ + k];

                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] += temp * b.r[coorBK + i];
                            }
                        }
                    }
                }
            }
        }
        else {
            // Form  B := alpha*B*A**T.
            if (upper) {
                for (let k = 1; k <= n; k++) {
                    const coorAK = a.colOfEx(k);
                    const coorBK = b.colOfEx(k);
                    for (let j = 1; j <= k - 1; j++) {
                        const coorBJ = b.colOfEx(j);
                        if (a.r[coorAK + j] !== 0) {
                            let temp = alpha * a.r[coorAK + j];
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] += temp * b.r[coorBK + i];
                            }
                        }
                    }
                    let temp = alpha
                    if (nounit) temp *= a.r[coorAK + k];
                    if (temp !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] = temp * b.r[coorBK + i];
                        }
                    }
                }
            }
            else {
                for (let k = n; k >= 1; k--) {
                    const coorAK = a.colOfEx(k);
                    const coorBK = b.colOfEx(k);
                    for (let j = k + 1; j <= n; j++) {
                        const coorBJ = b.colOfEx(j);
                        if (a.r[coorAK + j] !== 0) {
                            let temp = alpha * a.r[coorAK + j];
                            for (let i = 1; i <= m; i++) {
                                b.r[coorBJ + i] += temp * b.r[coorBK + i];
                            }
                        }
                    }
                    let temp = alpha;
                    if (nounit) temp *= a.r[coorAK + k];
                    if (temp !== 1) {
                        for (let i = 1; i <= m; i++) {
                            b.r[coorBK + i] *= temp;
                        }
                    }
                }
            }
        }
    }
}


