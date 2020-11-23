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

import type { Complex, MatrixEComplex } from '../../../f_func';

export function BtranA(
    nounit: boolean,
    upper: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex,
): void {
    if (upper) {
        for (let k = 1; k <= n; k++) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            for (let j = 1; j <= k - 1; j++) {
                const coorBJ = b.colOfEx(j);
                const aIsZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!aIsZero) {
                    //console.log(`${k},${j}`);
                    const ajkRe = a.r[coorAK + j];
                    const ajkIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];
                    //TEMP = ALPHA*A(J,K)
                    //TEMP = ALPHA*CONJG(A(J,K))
                    const tempRe = alpha.re * ajkRe - alpha.im * ajkIm;
                    const tempIm = alpha.re * ajkIm + alpha.im * ajkRe;
                    //console.log(`${k},${j},\t(${tempRe},${tempIm})`);
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) + TEMP*B(I,K)
                        const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                        // console.log(`${k},${j},${i}\t${re},${im}`);
                    }
                }
            } //j
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            if (nounit) {
                const akkRe = a.r[coorAK + k];
                const akkIm = noconj ? a.i[coorAK + k] : -a.i[coorAK + k];

                const re = tempRe * akkRe - tempIm * akkIm;
                const im = tempRe * akkIm + tempIm * akkRe;

                tempRe = re;
                tempIm = im;
            }
            if (!(tempRe === 1 && tempIm === 0)) {
                for (let i = 1; i <= m; i++) {
                    const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        } //k
    } //upper
    else {
        for (let k = n; k >= 1; k--) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            for (let j = k + 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const aIsZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!aIsZero) {
                    const ajkRe = a.r[coorAK + j];
                    const ajkIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];

                    const tempRe = alpha.re * ajkRe - alpha.im * ajkIm;
                    const tempIm = alpha.re * ajkIm + alpha.im * ajkRe;
                    for (let i = 1; i <= m; i++) {
                        //B(I,J) = B(I,J) + TEMP*B(I,K)
                        const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            } //for(j)
            let tempRe = alpha.re;
            let tempIm = alpha.im;

            if (nounit) {
                const akkRe = a.r[coorAK + k];
                const akkIm = noconj ? a.i[coorAK + k] : -a.i[coorAK + k];

                const re = tempRe * akkRe - tempIm * akkIm;
                const im = tempRe * akkIm + tempIm * akkRe;

                tempRe = re;
                tempIm = im;
            }
            if (!(tempRe === 1 && tempIm === 0)) {
                for (let i = 1; i <= m; i++) {
                    //  B(I,K) = TEMP*B(I,K)
                    const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        } //k
    } //upper
}
