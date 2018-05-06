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
    Complex,
    MatrixEComplex,
    mul_cxr,
    mul_rxr
} from '../../../f_func';


export function AB(
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
            for (let k = 1; k <= m; k++) {
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!bIsZero) {
                    // TEMP = ALPHA*B(K,J)
                    const coorAK = a.colOfEx(k);
                    const { re, im } = mul_cxr(alpha, b.r[coorBJ + k], b.i[coorBJ + k]);
                    let tempRe = re;
                    let tempIm = im;

                    for (let i = 1; i <= k - 1; i++) {
                        //TEMP*A(I,K)
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAK + i], a.i[coorAK + i]);
                        //B(I,J) = B(I,J) + ...
                        // console.log(`${j},${k},${i},\t (${re},${im})`);
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                        //                console.log(`${j},${k},${i},\t (${b.r[coorBJ + i]},${b.i[coorBJ + i]})`);

                    }
                    if (nounit) {
                        //TEMP = TEMP*A(K,K)
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAK + k], a.i[coorAK + k]);
                        tempRe = re;
                        tempIm = im;
                    }
                    //B(K,J) = TEMP
                    b.r[coorBJ + k] = tempRe;
                    b.i[coorBJ + k] = tempIm;
                }//if b[*]!=zero
            }//for k
        }//for j
    }
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let k = m; k >= 1; k--) {
                const coorAK = a.colOfEx(k);
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!bIsZero) {

                    let tempRe = alpha.re * b.r[coorBJ + k] - alpha.im * b.i[coorBJ + k];
                    let tempIm = alpha.re * b.i[coorBJ + k] + alpha.im * b.r[coorBJ + k];
                    b.r[coorBJ + k] = tempRe;
                    b.i[coorBJ + k] = tempIm;
                    //                    console.log(`${j},${k},\t${tempRe},${tempIm}`);
                    //IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                    if (nounit) {
                        /*const { re, im } = mul_rxr(
                            b.r[coorBJ + k],
                            b.i[coorBJ + k],
                            a.r[coorAK + k],
                            a.i[coorAK + k]
                        );*/
                        const re1 = b.r[coorBJ + k] * a.r[coorAK + k] - b.i[coorBJ + k] * a.i[coorAK + k];
                        const im1 = b.r[coorBJ + k] * a.i[coorAK + k] + b.i[coorBJ + k] * a.r[coorAK + k];
                        //     console.log(`${j},${k},\t${a.r[coorAK + k]},${a.i[coorAK + k]}`);
                        b.r[coorBJ + k] = re1;
                        b.i[coorBJ + k] = im1;
                    }
                    //console.log(`${j},${k},\t${b.r[coorBJ + k]},${b.i[coorBJ + k]}`);
                    for (let i = k + 1; i <= m; i++) {
                        //TEMP * A(I,K)
                        let re = tempRe * a.r[coorAK + i] - tempIm * a.i[coorAK + i];
                        let im = tempRe * a.i[coorAK + i] + tempIm * a.r[coorAK + i];
                        //B(I,J) = B(I,J) + TEMP*A(I,K)
                        //console.log(`${j},${k},${i}\t${re},${im}`);
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }
        }
    }
}
