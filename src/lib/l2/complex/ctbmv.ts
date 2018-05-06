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
    FortranArr,
    lowerChar,
    Matrix,
    mul_rxr
} from '../../f_func';
/*
*>
*> CTBMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.
*/

const { max, min } = Math;

export function ctbmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    // faster then String.toLowerCase()
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (!'un'.includes(dg)) {
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
        throw new Error(errWrongArg('ctbmv', info));
    }

    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        if (ul === 'u') {
            let kplus1 = k + 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                const xIsZero = x.r[jx - x.base] === 0
                    && x.i[jx - x.base] === 0;
                if (!xIsZero) {
                    //  TEMP = X(JX)
                    let tempRe = x.r[jx - x.base];
                    let tempIm = x.i[jx - x.base];
                    let ix = kx;
                    let L = kplus1 - j;
                    const coorAJ = a.colOfEx(j);
                    for (let i = max(1, j - k); i <= j - 1; i++) {
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAJ + L + i], a.i[coorAJ + L + i]);
                        x.r[ix - x.base] += re;
                        x.i[ix - x.base] += im;
                        ix += incx;
                    }
                    if (nounit) {
                        const { re, im } = mul_rxr(x.r[jx - x.base], x.i[jx - x.base], a.r[coorAJ + kplus1], a.i[coorAJ + kplus1]); // * a.r[coorAJ + kplus1] - x.i[jx - x.base] * a.i[coorAJ + kplus1];
                        x.r[jx - x.base] = re;
                        x.i[jx - x.base] = im;
                    }
                }
                jx += incx;
                if (j > k) kx += incx;

            } //for
        }
        // ul !== ['u']
        else {
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                const coorAJ = a.colOfEx(j);
                const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx - x.base];
                    let tempIm = x.i[jx - x.base];
                    let ix = kx;
                    let L = 1 - j;
                    for (let i = min(n, j + k); i >= j + 1; i--) {
                        // X(IX) = X(IX) + TEMP*A(L+I,J)
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAJ + L + i], a.i[coorAJ + L + i]);
                        x.r[ix - x.base] += re; //tempRe * a.r[coorAJ + L + i] - tempIm * a.i[coorAJ + L + i];
                        x.i[ix - x.base] += im; //tempIm * a.i[coorAJ + L + i] + tempRe * a.r[coorAJ + L + i];
                        ix -= incx;
                    }
                    if (nounit) {
                        const { re, im } = mul_rxr(x.r[jx - x.base], x.i[jx - x.base], a.r[coorAJ + 1], a.i[coorAJ + 1]);
                        x.r[jx - x.base] = re; //x.r[jx - x.base] * a.r[coorAJ + 1] - x.i[jx - x.base] * a.i[coorAJ + 1];
                        x.i[jx - x.base] = im; //x.r[jx - x.base] * a.i[coorAJ + 1] + x.i[jx - x.base] * a.r[coorAJ + 1];
                    }
                }
                jx -= incx;
                if (n - j >= k) kx -= incx;
            }
        }
    } else {
        // trans != 'n'    
        if (ul === 'u') {
            let kplus1 = k + 1;
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx - x.base];
                let tempIm = x.i[jx - x.base];

                kx -= incx;
                let ix = kx;
                let L = kplus1 - j;
                const coorAJ = a.colOfEx(j);
                const extrI = max(1, j - k); // evaluate once!
                const sign: 1 | -1 = noconj ? 1 : -1;
                if (nounit) {
                    // TEMP*A(1,J)
                    const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAJ + kplus1], sign * a.i[coorAJ + kplus1]);
                    tempRe = re;
                    tempIm = im;
                }
                for (let i = j - 1; i >= extrI; i--) {
                    //A(L+I,J)*X(IX)
                    const { re, im } = mul_rxr(a.r[coorAJ + L + i], sign * a.i[coorAJ + L + i], x.r[ix - x.base], x.i[ix - x.base]);
                    tempRe += re;
                    tempIm += im;
                    ix -= incx;
                }
                x.r[jx - x.base] = tempRe;
                x.i[jx - x.base] = tempIm;
                jx -= incx;
            }
        }
        else {
            //upl=lower
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                let tempRe = x.r[jx - x.base];
                let tempIm = x.i[jx - x.base];
                kx += incx;
                let ix = kx;
                const L = 1 - j;
                const coords = a.colOfEx(j);
                const extrI = min(n, j + k); // evaluate once!
                const sign = noconj ? 1 : -1;
                //   IF (NOUNIT) TEMP = TEMP*(N|C)A(1,J)
                if (nounit) {
                    const { re, im } = mul_rxr(tempRe, tempIm, a.r[coords + 1], sign * a.i[coords + 1]);
                    tempRe = re;
                    tempIm = im;
                }
                for (let i = j + 1; i <= extrI; i++) {
                    const { re, im } = mul_rxr(a.r[coords + L + i], sign * a.i[coords + L + i], x.r[ix - x.base], x.i[ix - x.base]);
                    tempRe += re;
                    tempIm += im;
                    ix += incx;
                }
                x.r[jx - x.base] = tempRe;
                x.i[jx - x.base] = tempIm;
                jx += incx;
            }
        }
    }
}
