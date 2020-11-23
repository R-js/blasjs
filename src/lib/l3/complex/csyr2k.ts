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

import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function csyr2k(
    uplo: 'u' | 'l',
    trans: 'n' | 't',
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: Complex,
    c: Matrix,
    ldc: number,
): void {
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    const tr = lowerChar(trans);
    const ul = lowerChar(uplo);

    const nrowA = tr === 'n' ? n : k;
    const upper = ul === 'u';

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (!'nt'.includes(tr)) {
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
    if (info !== 0) {
        throw new Error(errWrongArg('csyr2k', info));
    }

    //Quick return if possible.

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;
    const betaIsZero = beta.re === 0 && beta.im === 0;

    if (n === 0 || ((alphaIsZero || k === 0) && betaIsOne)) return;

    //And when  alpha.eq.zero.

    if (alphaIsZero) {
        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                const re = betaIsZero ? 0 : beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = betaIsZero ? 0 : beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];

                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }
        return;
    }

    //     Start the operations.

    if (tr === 'n') {
        //  Form  C := alpha*A*B**T + alpha*B*A**T + C.

        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            if (betaIsZero) {
                c.setCol(j, start, stop, 0);
            } else if (!betaIsOne) {
                for (let i = start; i <= stop; i++) {
                    c.r[coorCJ + i] = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    c.i[coorCJ + i] = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                const coorBL = b.colOfEx(l);
                const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j] === 0;
                const bIsZero = b.r[coorBL + j] === 0 && b.i[coorBL + j] === 0;
                if (!aIsZero || !bIsZero) {
                    // TEMP1 = ALPHA * B(J, L)
                    const temp1Re = alpha.re * b.r[coorBL + j] - alpha.im * b.i[coorBL + j];
                    const temp1Im = alpha.re * b.i[coorBL + j] + alpha.im * b.r[coorBL + j];
                    // TEMP2 = ALPHA * A(J, L)
                    const temp2Re = alpha.re * a.r[coorBL + j] - alpha.im * a.i[coorBL + j];
                    const temp2Im = alpha.re * a.i[coorBL + j] + alpha.im * a.r[coorBL + j];

                    for (let i = start; i <= stop; i++) {
                        //A(I, L) * TEMP1
                        let at1Re = a.r[coorAL + i] * temp1Re - a.i[coorAL + i] * temp1Im;
                        let at1Im = a.r[coorAL + i] * temp1Im + a.i[coorAL + i] * temp1Re;
                        //B(I, L) * TEMP2
                        at1Re += b.r[coorBL + i] * temp2Re - b.i[coorBL + i] * temp2Im;
                        at1Im += b.r[coorBL + i] * temp2Im + b.i[coorBL + i] * temp2Re;
                        // C(I, J) = C(I, J) + ....
                        c.r[coorCJ + i] += at1Re;
                        c.i[coorCJ + i] += at1Im;
                    }
                }
            }
        } //for(j)
    } //tr === 't'
    else {
        //Form  C := alpha*A**T*B + alpha*B**T*A + C.

        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorBJ = b.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                const coorAI = a.colOfEx(i);
                const coorBI = b.colOfEx(i);
                let temp1Re = 0;
                let temp1Im = 0;
                let temp2Re = 0;
                let temp2Im = 0;
                for (let l = 1; l <= k; l++) {
                    //TEMP1 = TEMP1 + A(L,I)*B(L,J)
                    temp1Re += a.r[coorAI + l] * b.r[coorBJ + l] - a.i[coorAI + l] * b.i[coorBJ + l];
                    temp1Im += a.r[coorAI + l] * b.i[coorBJ + l] + a.i[coorAI + l] * b.r[coorBJ + l];
                    //TEMP2 = TEMP2 + B(L,I)*A(L,J)
                    temp2Re += b.r[coorBI + l] * a.r[coorAJ + l] - b.i[coorBI + l] * a.i[coorAJ + l];
                    temp2Im += b.r[coorBI + l] * a.i[coorAJ + l] + b.i[coorBI + l] * a.r[coorAJ + l];
                }
                //console.log(`j:${j},i:${i},\t(${temp2Re},${temp2Im})`);
                //ALPHA*TEMP1 + ALPHA*TEMP2
                let re = alpha.re * (temp1Re + temp2Re) - alpha.im * (temp1Im + temp2Im);
                let im = alpha.re * (temp1Im + temp2Im) + alpha.im * (temp1Re + temp2Re);
                //console.log(`j:${j},i:${i},\t(${re},${im})`);

                if (!betaIsZero) {
                    re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            } //for(i)
        } //for(j)
    } //tr==='t'
}
