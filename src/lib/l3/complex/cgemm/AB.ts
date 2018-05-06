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

export function AB(
    betaIsZero: boolean,
    betaIsOne: boolean,
    beta: Complex,
    alpha: Complex,
    a: MatrixEComplex,
    b: MatrixEComplex,
    c: MatrixEComplex,
    n: number,
    m: number,
    k: number): void {

    // DO 90 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        const coorBJ = b.colOfEx(j);
        // IF (BETA.EQ.ZERO) THEN
        if (betaIsZero) {
            c.setCol(j, 1, m, 0);
        }
        //IF (BETA.NE.ONE) THEN
        else if (!betaIsOne) {
            /// DO 60 I = 1,M
            for (let i = 1; i <= m; i++) {
                const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }
        //DO 80 L = 1,K
        for (let l = 1; l <= k; l++) {
            const coorAL = a.colOfEx(l);
            let tempRe = alpha.re * b.r[coorBJ + l] - alpha.im * b.i[coorBJ + l];
            let tempIm = alpha.re * b.i[coorBJ + l] + alpha.im * b.r[coorBJ + l];

            for (let i = 1; i <= m; i++) {
                const re3 = tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                const im3 = tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
                c.r[coorCJ + i] += re3;
                c.i[coorCJ + i] += im3;
            }
        }
    }
}
