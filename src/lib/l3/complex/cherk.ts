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

import {
    errMissingIm,
    errWrongArg,
    lowerChar,
    Matrix
} from '../../f_func';
const { max } = Math;

export function cherk(
    uplo: 'u' | 'l',
    trans: 'n' | 'c',
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    beta: number,
    c: Matrix,
    ldc: number): void {



    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
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
    }
    else if (!'nc'.includes(tr)) {
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
    if (info !== 0) {
        throw new Error(errWrongArg('cherk', info));
    }

    //Quick return if possible.

    if (n === 0 || ((alpha === 0 || k === 0) && beta === 1)) return;


    //And when  alpha.eq.zero.

    if (alpha === 0) {

        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = beta === 0 ? 0 : c.r[coorCJ + i] * beta;
                c.i[coorCJ + i] = (beta === 0 || j === i) ? 0 : c.i[coorCJ + i] * beta;
            }
        }
        return;
    }

    //*     Start the operations.

    if (tr === 'n') {
        // Form  C := alpha*A*A**H + beta*C.
        // upper is weaved into the process
        for (let j = 1; j <= n; j++) {
            const fs = upper ? 1 : j;
            const fe = upper ? j : n;
            const start = upper ? 1 : j + 1;
            const stop = upper ? j - 1 : n;
            const coorCJ = c.colOfEx(j);
            if (beta === 0) {
                c.setCol(j, fs, fe, 0);
            }
            else if (beta !== 1) {
                for (let i = start; i <= stop; i++) {
                    c.r[coorCJ + i] = beta * c.r[coorCJ + i];
                    c.i[coorCJ + i] = beta * c.i[coorCJ + i];
                }
                c.r[coorCJ + j] = beta * c.r[coorCJ + j];
                c.i[coorCJ + j] = 0;
            }
            else {
                c.i[coorCJ + j] = 0;
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j] === 0;
                if (!aIsZero) {
                    // TEMP = ALPHA*DCONJG(A(J,L))
                    let tempRe = alpha * a.r[coorAL + j];
                    let tempIm = -alpha * a.i[coorAL + j];
                    for (let i = start; i <= stop; i++) {
                        //  C(I,J) = C(I,J) + TEMP*A(I,L)
                        c.r[coorCJ + i] += (tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i])
                        c.i[coorCJ + i] += (tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i]);
                    }
                    c.r[coorCJ + j] += (tempRe * a.r[coorAL + j] - tempIm * a.i[coorAL + j])
                    c.i[coorCJ + j] = 0;
                }//if
            }//forl
        }//forj
    }
    //tr === "c"
    else {
        // Form  C := alpha*A**H*A + beta*C.
        if (upper) {
            for (let j = 1; j <= n; j++) {
                const coorAJ = a.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= j - 1; i++) {
                    const coorAI = a.colOfEx(i);
                    let tempRe = 0;
                    let tempIm = 0;
                    for (let l = 1; l <= k; l++) {
                        //(a-ib)*(c+id)=(ac+bd)+i(ad-bc)
                        // CONJG(A(L,I))*A(L,J)
                        tempRe += a.r[coorAI + l] * a.r[coorAJ + l] + a.i[coorAI + l] * a.i[coorAJ + l];
                        tempIm += a.r[coorAI + l] * a.i[coorAJ + l] - a.i[coorAI + l] * a.r[coorAJ + l];
                    }
                    //ALPHA*TEMP
                    let re = alpha * tempRe;
                    let im = alpha * tempIm;

                    if (beta !== 0) {
                        //[ALPHA*TEMP] + BETA*C(I,J)
                        re += beta * c.r[coorCJ + i];
                        im += beta * c.i[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
                let rTempRe = 0;
                let rTempIm = 0;
                for (let l = 1; l <= k; l++) {
                    //RTEMP = RTEMP + CONJG(A(L,J))*A(L,J)
                    //(a-ib)(c+id)=(ac+bd)+i(ad-bc)
                    const alj = coorAJ + l;
                    rTempRe += a.r[alj] * a.r[alj] + a.i[alj] * a.i[alj];
                    rTempIm += a.r[alj] * a.i[alj] - a.i[alj] * a.r[alj];
                }
                let re = alpha * rTempRe;
                let im = alpha * rTempIm;
                if (beta !== 0) {
                    //[ALPHA*RTEMP] + BETA*REAL(C(J,J))
                    re += beta * c.r[coorCJ + j];
                    //im += beta* c.i[coorCJ + j];
                }
                c.r[coorCJ + j] = re;
                c.i[coorCJ + j] = im;
            }
        }
        // not upper
        else {
            for (let j = 1; j <= n; j++) {
                let rTempRe = 0;

                const coorAJ = a.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let l = 1; l <= k; l++) {
                    //RTEMP = RTEMP + CONJG(A(L,J))*A(L,J)
                    // (a-ib)*(c+id)=(ac+ad)+i(ad-bc)
                    const ajl = coorAJ + l;
                    const ar = a.r[ajl];
                    const ai = a.i[ajl];
                    rTempRe += ar * ar + ai * ai;
                    //   rTempIm += ar * ai - ai * ar;
                }
                //console.log(`${j}\t${rTempRe}`)
                //ALPHA*RTEMP
                let re = alpha * rTempRe;
                //let im = alpha * rTempIm;
                if (beta !== 0) {
                    re += beta * c.r[coorAJ + j];
                    //im += beta * c.i[coorAJ + j];
                }
                c.r[coorAJ + j] = re;
                c.i[coorAJ + j] = 0;

                for (let i = j + 1; i <= n; i++) {
                    const coorAI = a.colOfEx(i);
                    let tempRe = 0;
                    let tempIm = 0;
                    for (let l = 1; l <= k; l++) {
                        //TEMP = TEMP + CONJG(A(L,I))*A(L,J)
                        //(a-ib)*(c+id)=(ac+bd)+i(ad-bd)
                        const alj = coorAJ + l;
                        const ali = coorAI + l;
                        tempRe += a.r[ali] * a.r[alj] + a.i[ali] * a.i[alj];
                        tempIm += a.r[ali] * a.i[alj] - a.i[ali] * a.r[alj];
                    }
                    //ALPHA*TEMP
                    let re = alpha * tempRe;
                    let im = alpha * tempIm;
                    if (beta !== 0) {
                        re += beta * c.r[coorCJ + i];
                        im += beta * c.i[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }//for(i)
            }// for(j)
        }// not upper
    }//
}
