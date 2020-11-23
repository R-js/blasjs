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

import { Complex, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';

export function conjAtransB(
    betaIsZero: boolean,
    _betaIsOne: boolean,
    beta: Complex,
    alpha: Complex,
    a: MatrixEComplex,
    b: MatrixEComplex,
    c: MatrixEComplex,
    n: number,
    m: number,
    k: number,
): void {
    // DO 310 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        // DO 300 I = 1,M
        for (let i = 1; i <= m; i++) {
            const coorAI = a.colOfEx(i);
            // TEMP = ZERO
            let tempRe = 0;
            let tempIm = 0;
            //DO 290 L = 1,K
            for (let l = 1; l <= k; l++) {
                const coorBL = b.colOfEx(l);
                // TEMP = TEMP + CONJG(A(L,I)) * B(J,L)
                // (a-ib)*(c+id) = (ac+bd) + i(ad-bc)
                const { re, im } = mul_rxr(a.r[coorAI + l], -a.i[coorAI + l], b.r[coorBL + j], b.i[coorBL + j]);
                tempRe += re;
                tempIm += im;
            }
            /*
            IF (BETA.EQ.ZERO) THEN
              C(I,J) = ALPHA*TEMP
            ELSE
              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
            END IF
            */
            // C(I,J) = ALPHA*TEMP
            let { re, im } = mul_cxr(alpha, tempRe, tempIm);
            if (!betaIsZero) {
                const { re: re1, im: im1 } = mul_cxr(beta, c.r[coorCJ + i], c.i[coorCJ + i]);
                re += re1;
                im += im1;
            }
            c.r[coorCJ + i] = re;
            c.i[coorCJ + i] = im;
        }
    }
}
