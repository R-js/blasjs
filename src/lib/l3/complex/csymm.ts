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
import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix, mul_cxr, mul_rxr } from '../../f_func';

const { max } = Math;

export function csymm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    m: number,
    n: number,
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

    const si = lowerChar(side);
    const ul = lowerChar(uplo);

    //Set NROWA as the number of rows of A.

    const nrowA = si === 'l' ? m : n;
    const upper = ul === 'u';

    //Test the input parameters.

    let info = 0;

    if (!'lr'.includes(si)) {
        info = 1;
    } else if (!'ul'.includes(ul)) {
        info = 2;
    } else if (m < 0) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (lda < max(1, nrowA)) {
        info = 7;
    } else if (ldb < max(1, m)) {
        info = 9;
    } else if (ldc < max(1, m)) {
        info = 12;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('csymm', info));
    }

    //Quick return if possible.
    const alphaIsZero = alpha.re === 0 && alpha.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;
    const betaIsZero = beta.re === 0 && beta.im === 0;

    if (m === 0 || n === 0 || (alphaIsZero && betaIsOne)) return;

    //And when  alpha.eq.zero.

    if (alphaIsZero) {
        if (betaIsZero) {
            // console.log('Setting all to zero!!');
            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        } else {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
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

    if (si === 'l') {
        // Form  C := alpha*A*B + beta*C
        if (upper) {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                const coorBJ = b.colOfEx(j);

                for (let i = 1; i <= m; i++) {
                    const coorAI = a.colOfEx(i);
                    //TEMP1 = ALPHA*B(I,J)
                    const temp1Re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const temp1Im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];

                    let temp2Re = 0;
                    let temp2Im = 0;

                    for (let k = 1; k <= i - 1; k++) {
                        // C(K,J) = C(K,J) + TEMP1*A(K,I)
                        c.r[coorCJ + k] += temp1Re * a.r[coorAI + k] - temp1Im * a.i[coorAI + k];
                        c.i[coorCJ + k] += temp1Re * a.i[coorAI + k] + temp1Im * a.r[coorAI + k];
                        //TEMP2 = TEMP2 + B(K,J)*A(K,I)
                        temp2Re += b.r[coorBJ + k] * a.r[coorAI + k] - b.i[coorBJ + k] * a.i[coorAI + k];
                        temp2Im += b.r[coorBJ + k] * a.i[coorAI + k] + b.i[coorBJ + k] * a.r[coorAI + k];
                    }
                    //C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2
                    let re =
                        temp1Re * a.r[coorAI + i] -
                        temp1Im * a.i[coorAI + i] +
                        (alpha.re * temp2Re - alpha.im * temp2Im);
                    let im =
                        temp1Re * a.i[coorAI + i] +
                        temp1Im * a.r[coorAI + i] +
                        (alpha.re * temp2Im + alpha.im * temp2Re);
                    if (!betaIsZero) {
                        // BETA*C(I,J)
                        re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                        im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
        //not upper
        else {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = m; i >= 1; i--) {
                    const coorAI = a.colOfEx(i);
                    //TEMP1 = ALPHA * B(I, J)
                    //let temp1Re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    //let temp1Im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    const { re, im } = mul_cxr(alpha, b.r[coorBJ + i], b.i[coorBJ + i]);
                    const temp1Re = re;
                    const temp1Im = im;
                    let temp2Re = 0;
                    let temp2Im = 0;
                    for (let k = i + 1; k <= m; k++) {
                        // C(K,J) = C(K,J) + TEMP1*A(K,I)
                        const { re, im } = mul_rxr(temp1Re, temp1Im, a.r[coorAI + k], a.i[coorAI + k]);

                        c.r[coorCJ + k] += re; //temp1Re * a.r[coorAI + k] - temp1Im * a.i[coorAI + k];
                        c.i[coorCJ + k] += im; //temp1Re * a.i[coorAI + k] + temp1Im * a.r[coorAI + k];
                        // TEMP2 = TEMP2 + B(K,J)*A(K,I)
                        const { re: re1, im: im1 } = mul_rxr(
                            b.r[coorBJ + k],
                            b.i[coorBJ + k],
                            a.r[coorAI + k],
                            a.i[coorAI + k],
                        );
                        temp2Re += re1; //b.r[coorBJ + k] * a.r[coorAI + k] - b.i[coorBJ + k] * a.i[coorAI + k];
                        temp2Im += im1; //b.r[coorBJ + k] * a.i[coorAI + k] + b.i[coorBJ + k] * a.r[coorAI + k];
                    }
                    //TEMP1*A(I,I) + ALPHA*TEMP2
                    const { re: re2, im: im2 } = mul_rxr(temp1Re, temp1Im, a.r[coorAI + i], a.i[coorAI + i]);
                    const { re: re3, im: im3 } = mul_cxr(alpha, temp2Re, temp2Im);
                    let re0 = re2 + re3;
                    let im0 = im2 + im3;
                    //let re = (temp1Re * a.r[coorAI + i] - temp1Im * a.i[coorAI + i]) + (alpha.re * temp2Re - alpha.im * temp2Im);
                    //let im = (temp1Re * a.i[coorAI + i] + temp1Im * a.r[coorAI + i]) + (alpha.re * temp2Im + alpha.im * temp2Re);
                    //BETA * C(I, J)
                    if (!betaIsZero) {
                        const { re, im } = mul_cxr(beta, c.r[coorCJ + i], c.i[coorCJ + i]);
                        re0 += re;
                        im0 += im;
                        //  re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                        //  im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re0;
                    c.i[coorCJ + i] = im0;
                } //for(i)
            } //for(j)
        } //upper
    }
    //si==='r'
    else {
        // Form  C := alpha*B*A + beta*C.
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorBJ = b.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            //TEMP1 = ALPHA*A(J,J)
            let temp1Re = alpha.re * a.r[coorCJ + j] - alpha.im * a.i[coorCJ + j];
            let temp1Im = alpha.re * a.i[coorCJ + j] + alpha.im * a.r[coorCJ + j];
            for (let i = 1; i <= m; i++) {
                //TEMP1*B(I,J)
                let re = temp1Re * b.r[coorBJ + i] - temp1Im * b.i[coorBJ + i];
                let im = temp1Re * b.i[coorBJ + i] + temp1Im * b.r[coorBJ + i];
                //console.log(`${i},${j}\t(${re},${im})`);
                if (!betaIsZero) {
                    //BETA*C(I,J)
                    re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
            for (let k = 1; k <= j - 1; k++) {
                const coorAK = a.colOfEx(k);
                const coorBK = b.colOfEx(k);
                //A(K,J) : A(J,K)
                const aRe = upper ? a.r[coorAJ + k] : a.r[coorAK + j];
                const aIm = upper ? a.i[coorAJ + k] : a.i[coorAK + j];

                temp1Re = alpha.re * aRe - alpha.im * aIm;
                temp1Im = alpha.re * aIm + alpha.im * aRe;

                for (let i = 1; i <= m; i++) {
                    //C(I,J) = C(I,J) + TEMP1*B(I,K)
                    c.r[coorCJ + i] += temp1Re * b.r[coorBK + i] - temp1Im * b.i[coorBK + i];
                    c.i[coorCJ + i] += temp1Re * b.i[coorBK + i] + temp1Im * b.r[coorBK + i];
                }
            }
            for (let k = j + 1; k <= n; k++) {
                const coorAK = a.colOfEx(k);
                const coorBK = b.colOfEx(k);
                //A(K,J) : A(J,K)
                const aRe = upper ? a.r[coorAK + j] : a.r[coorAJ + k];
                const aIm = upper ? a.i[coorAK + j] : a.i[coorAJ + k];

                temp1Re = alpha.re * aRe - alpha.im * aIm;
                temp1Im = alpha.re * aIm + alpha.im * aRe;

                for (let i = 1; i <= m; i++) {
                    //C(I,J) = C(I,J) + TEMP1*B(I,K)
                    c.r[coorCJ + i] += temp1Re * b.r[coorBK + i] - temp1Im * b.i[coorBK + i];
                    c.i[coorCJ + i] += temp1Re * b.i[coorBK + i] + temp1Im * b.r[coorBK + i];
                }
            }
        } //for(j)
    } //si==='r'
}
