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

export function csyrk(
    uplo: 'u' | 'l',
    trans: 't' | 'n',
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    beta: Complex,
    c: Matrix,
    ldc: number,
): void {
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    const tr = lowerChar(trans);
    const ul = lowerChar(uplo);
    const alphaIsZero = alpha.re === 0 && alpha.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;
    const betaIsZero = beta.re === 0 && beta.im === 0;

    const upper = ul === 'u';
    const nrowA = tr === 'n' ? n : k;

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (!'tn'.includes(tr)) {
        info = 2;
    } else if (n < 0) {
        info = 3;
    } else if (k < 0) {
        info = 4;
    } else if (lda < max(1, nrowA)) {
        info = 7;
    } else if (ldc < max(1, n)) {
        info = 10;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('csyrk', info));
    }

    //  Quick return if possible.
    if (n === 0 || ((alphaIsZero || k === 0) && betaIsOne)) return;

    //  And when  alpha.eq.zero.

    if (alphaIsZero) {
        for (let j = 1; j <= n; j++) {
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            if (betaIsZero) {
                c.setCol(j, start, stop, 0);
            } else {
                const coorCJ = c.colOfEx(j);
                for (let i = start; i <= stop; i++) {
                    const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
        return;
    }

    // Start the operations.
    if (tr === 'n') {
        // Form  C := alpha*A*A**T + beta*C.
        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            if (betaIsZero) {
                c.setCol(j, start, stop, 0);
            } else if (!betaIsOne) {
                for (let i = start; i <= stop; i++) {
                    const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j] === 0;
                if (!aIsZero) {
                    const tempRe = alpha.re * a.r[coorAL + j] - alpha.im * a.i[coorAL + j];
                    const tempIm = alpha.re * a.i[coorAL + j] + alpha.im * a.r[coorAL + j];
                    for (let i = start; i <= stop; i++) {
                        c.r[coorCJ + i] += tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                        c.i[coorCJ + i] += tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
                    }
                }
            }
        }
    } else {
        //  Form  C := alpha*A**T*A + beta*C.
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                const coorAI = a.colOfEx(i);
                let tempRe = 0;
                let tempIm = 0;
                for (let l = 1; l <= k; l++) {
                    tempRe += a.r[coorAI + l] * a.r[coorAJ + l] - a.i[coorAI + l] * a.i[coorAJ + l];
                    tempIm += a.r[coorAI + l] * a.i[coorAJ + l] + a.i[coorAI + l] * a.r[coorAJ + l];
                }
                //ALPHA*TEMP
                let re = alpha.re * tempRe - alpha.im * tempIm;
                let im = alpha.re * tempIm + alpha.im * tempRe;
                if (!betaIsZero) {
                    re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }
    }
}
