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

import { Complex, div_rxr, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';
//Form  B := alpha*inv( A )*B.

export function invAB(
    nounit: boolean,
    upper: boolean,
    alphaIsOne: boolean,
    //alphaIsZero: boolean,
    //noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex,
): void {
    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    //const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    //const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    const { re, im } = mul_cxr(alpha, b.r[coorBJ + i], b.i[coorBJ + i]);
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                    //console.log(`${j},${i}\t(${re},${im})`)
                }
            }
            for (let k = m; k >= 1; k--) {
                const coorAK = a.colOfEx(k);
                const isBZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!isBZero) {
                    if (nounit) {
                        const { re, im } = div_rxr(b.r[coorBJ + k], b.i[coorBJ + k], a.r[coorAK + k], a.i[coorAK + k]);
                        b.r[coorBJ + k] = re;
                        b.i[coorBJ + k] = im;
                        //       console.log(`*${j},${k}\t(${b.r[coorBJ + k]},${b.i[coorBJ + k]})`);
                    }
                    for (let i = 1; i <= k - 1; i++) {
                        const { re, im } = mul_rxr(b.r[coorBJ + k], b.i[coorBJ + k], a.r[coorAK + i], a.i[coorAK + i]);
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    } //for(i)
                } //if b(k,j)!=0
            } //k
        } //j
    } //upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    /*const { re, im } = mul_cxr(
                         alpha,
                         b.r[coorBJ + i],
                         b.i[coorBJ + i]
                     );*/
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = 1; k <= m; k++) {
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                const coorAK = a.colOfEx(k);
                if (!bIsZero) {
                    if (nounit) {
                        const { re, im } = div_rxr(b.r[coorBJ + k], b.i[coorBJ + k], a.r[coorAK + k], a.i[coorAK + k]);
                        b.r[coorBJ + k] = re;
                        b.i[coorBJ + k] = im;
                    }
                    for (let i = k + 1; i <= m; i++) {
                        /*const { re, im } = mul_rxr(
                            b.r[coorBJ + k],
                            b.i[coorBJ + k],
                            a.r[coorAK + i],
                            a.i[coorAK + i]
                        );*/
                        const re = b.r[coorBJ + k] * a.r[coorAK + i] - b.i[coorBJ + k] * a.i[coorAK + i];
                        const im = b.r[coorBJ + k] * a.i[coorAK + i] + b.i[coorBJ + k] * a.r[coorAK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                } //bNoZero
            } //k
        } //j
    } //upper
}
