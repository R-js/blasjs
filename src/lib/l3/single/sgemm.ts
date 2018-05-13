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

export function sgemm(
    transA: 'n' | 't' | 'c',
    transB: 'n' | 't' | 'c',
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: number,
    c: Matrix,
    ldc: number): void {

    // faster then String.toLowerCase()
    const trA = lowerChar(transA);
    const trB = lowerChar(transB);

    const notA = trA === 'n';
    const notB = trB === 'n';

    const nrowA = notA ? m : k;
    //ncolA is never used, I checked in the original F-code
    // also checked online http://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html#gafe51bacb54592ff5de056acabd83c260
    // ncolA is not used
    // const ncolA = notA ? k : n;

    const nrowB = notB ? k : n;

    let info = 0;
    if (!'ntc'.includes(trA)) {
        info = 1;
    }
    else if (!'ntc'.includes(trB)) {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (k < 0) {
        info = 5;
    }
    else if (lda < max(1, nrowA)) {
        info = 8;
    }
    else if (ldb < max(1, nrowB)) {
        info = 10;
    }
    else if (ldc < max(1, m)) {
        info = 13;
    }
    // ok?
    if (info !== 0) {
        throw new Error(errWrongArg('sgemm', info));
    }

    //*     Quick return if possible.
    if (m === 0 || n === 0 ||
        (
            (alpha === 0 || k === 0) && beta === 1
        )
    ) return;


    if (alpha === 0) {
        if (beta === 0) {
            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        }
        else {
            // I have to typecast it this way for TS compiler not to nag!!
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] *= beta;
                }
            }
        }
        return;
    }

    //*     Start the operations.
    if (notB) {
        if (notA) {
            //Form  C := alpha*A*B + beta*C.
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                const coorBJ = b.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, 1, m, 0);
                }
                else if (beta !== 1) {
                    for (let i = 1; i <= m; i++) {
                        c.r[coorCJ + i] *= beta;
                    }
                }
                for (let l = 1; l <= k; l++) {
                    let temp = alpha * b.r[coorBJ + l];
                    const coorAL = a.colOfEx(l);
                    for (let i = 1; i <= m; i++) {
                        c.r[coorCJ + i] += temp * a.r[coorAL + i]
                    }
                }
            }

        }
        else {
            //transA !== 'n'
            //     Form  C := alpha*A**T*B + beta*C
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    let temp = 0;
                    const coorA = a.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        temp += a.r[L + coorA] * b.r[coorBJ + L];
                    }
                    if (beta === 0) {
                        c.r[coorCJ + i] = alpha * temp;
                    }
                    else {
                        c.r[coorCJ + i] = alpha * temp + beta * c.r[coorCJ + i];
                    }
                }
            }
        }
    }
    else {
        if (notA) {
            // Form  C := alpha*A*B**T + beta*C
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, 1, m, 0);
                } else if (beta !== 1) {
                    for (let i = 1; i <= m; i++) {
                        c.r[coorCJ + i] *= beta;
                    }
                }
                for (let l = 1; l <= k; l++) {
                    let temp = alpha * b.r[(l - b.colBase) * b.nrRows - b.rowBase + j];
                    const coorAL = a.colOfEx(l);
                    for (let i = 1; i <= m; i++) {
                        c.r[coorCJ + i] += temp * a.r[coorAL + i];
                    }
                }
            }
        }
        else {
            //  Form  C := alpha*A**T*B**T + beta*C
            for (let j = 1; j <= n; j++) {
                for (let i = 1; i <= m; i++) {
                    let temp = 0;
                    const coorA = a.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        //TODO: write out colOfEx
                        temp += a.r[coorA + L] * b.r[b.colOfEx(L) + j];
                    }
                    const coorC = c.colOfEx(j);
                    if (beta === 0) {
                        c.r[coorC + i] = alpha * temp;
                    }
                    else {
                        c.r[coorC + i] = alpha * temp + beta * c.r[coorC + i];
                    }
                }
            }
        }
    }
}
