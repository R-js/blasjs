
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


import { Complex, MatrixEComplex } from '../../../f_func';

export function tranAB(
    nounit: boolean,
    upper: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex): void {



    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);

            for (let i = m; i >= 1; i--) {
                const coorAI = a.colOfEx(i);
                let tempRe = b.r[coorBJ + i];
                let tempIm = b.i[coorBJ + i];
                /*if (noconj) {
                    // IF (NOUNIT) TEMP = TEMP*A(I,I)
                    if (nounit) {
                        let re = tempRe * a.r[coorAI + i] - tempIm * a.i[coorAI + i];
                        let im = tempRe * a.i[coorAI + i] + tempIm * a.r[coorAI + i];
                        tempRe = re;
                        tempIm = im;
                    }
                    for (let k = 1; k <= i - 1; k++) {
                        //A(K,I)*B(K,J)
                        tempRe += a.r[coorAI + k] * b.r[coorBJ + k] - a.i[coorAI + k] * b.i[coorBJ + k];
                        tempIm += a.r[coorAI + k] * b.i[coorBJ + k] + a.i[coorAI + k] * b.r[coorBJ + k];
                    }
                }
                else {*/
                /*
                IF (NOUNIT) TEMP = TEMP*CONJG(A(I,I))
                DO 100 K = 1,I - 1
                    TEMP = TEMP + CONJG(A(K,I))*B(K,J)
100                 CONTINUE
                */
                if (nounit) {
                    // TEMP = TEMP*CONJG(A(I,I))
                    //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                    const re = tempRe * a.r[coorAI + i] - tempIm * (noconj ? a.i[coorAI + i] : -a.i[coorAI + i]);
                    const im = tempRe * (noconj ? a.i[coorAI + i] : -a.i[coorAI + i]) + tempIm * a.r[coorAI + i];
                    tempRe = re;
                    tempIm = im;
                }
                for (let k = 1; k <= i - 1; k++) {
                    // CONJG(A(K,I))*B(K,J)
                    // (a-ib)*(c+id) = (ab+bd)+i(ad-bc)
                    tempRe += a.r[coorAI + k] * b.r[coorBJ + k] - (noconj ? a.i[coorAI + k] : -a.i[coorAI + k]) * b.i[coorBJ + k];
                    tempIm += a.r[coorAI + k] * b.i[coorBJ + k] + (noconj ? a.i[coorAI + k] : -a.i[coorAI + k]) * b.r[coorBJ + k];
                }
                //}
                // B(I,J) = ALPHA*TEMP
                b.r[coorBJ + i] = alpha.re * tempRe - alpha.im * tempIm;
                b.i[coorBJ + i] = alpha.re * tempIm + alpha.im * tempRe;

            }//for(i)
        }//for(j)
    }
    //upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                const coorAI = a.colOfEx(i);
                // TEMP = B(I,J)
                let tempRe = b.r[coorBJ + i];
                let tempIm = b.i[coorBJ + i];
                /* if (noconj) {
                     // IF (NOUNIT) TEMP = TEMP*A(I,I)
                     if (nounit) {
                         let re = tempRe * a.r[coorAI + i] - tempIm * a.i[coorAI + i];
                         let im = tempRe * a.i[coorAI + i] + tempIm * a.r[coorAI + i];
                         tempRe = re;
                         tempIm = im;
                     }
                     for (let k = i + 1; k <= m; k++) {
                         tempRe += a.r[coorAI + k] * b.r[coorBJ + k] - a.i[coorAI + k] * b.i[coorBJ + k];
                         tempIm += a.r[coorAI + k] * b.i[coorBJ + k] + a.i[coorAI + k] * b.r[coorBJ + k];
                     }
                 }
                 else {*/
                //  IF (NOUNIT) TEMP = TEMP*CONJG(A(I,I))
                if (nounit) {
                    //(a-ib)*(c+id) =(ac+bd)+i(ad-bc)
                    let re = tempRe * a.r[coorAI + i] - tempIm * (noconj ? a.i[coorAI + i] : -a.i[coorAI + i]);
                    let im = tempRe * (noconj ? a.i[coorAI + i] : -a.i[coorAI + i]) + tempIm * a.r[coorAI + i];
                    tempRe = re;
                    tempIm = im;
                }
                for (let k = i + 1; k <= m; k++) {
                    // TEMP = TEMP + CONJG(A(K,I))*B(K,J)
                    //(a-ib)*(c+id) =(ac+bd)+i(ad-bc)
                    tempRe += a.r[coorAI + k] * b.r[coorBJ + k] - (noconj ? a.i[coorAI + k] : -a.i[coorAI + k]) * b.i[coorBJ + k];
                    tempIm += a.r[coorAI + k] * b.i[coorBJ + k] + (noconj ? a.i[coorAI + k] : -a.i[coorAI + k]) * b.r[coorBJ + k];
                }

                //}
                //B(I,J) = ALPHA*TEMP
                b.r[coorBJ + i] = alpha.re * tempRe - alpha.im * tempIm;
                b.i[coorBJ + i] = alpha.re * tempIm + alpha.im * tempRe;


            }
        }
    }
}
