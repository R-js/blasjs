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


export function BA(
    nounit: boolean,
    upper: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex): void {



    if (upper) {
        for (let j = n; j >= 1; j--) {
            const coorAJ = a.colOfEx(j);
            const coorBJ = b.colOfEx(j);
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            //IF (NOUNIT) TEMP = TEMP*A(J,J)
            if (nounit) {
                let re = tempRe * a.r[coorAJ + j] - tempIm * a.i[coorAJ + j];
                let im = tempRe * a.i[coorAJ + j] + tempIm * a.r[coorAJ + j];
                tempRe = re;
                tempIm = im;
            }
            for (let i = 1; i <= m; i++) {
                // B(I,J) = TEMP*B(I,J)
                // (a-ib)*(c+id) = (ab+bd)+i(ad-bc)
                let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                b.r[coorBJ + i] = re;
                b.i[coorBJ + i] = im;
            }
            //DO 190 K = 1,J - 1
            for (let k = 1; k <= j - 1; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {

                    //TEMP = ALPHA * A(K, J)
                    let tempRe = alpha.re * a.r[coorAJ + k] - alpha.im * a.i[coorAJ + k];
                    let tempIm = alpha.re * a.i[coorAJ + k] + alpha.im * a.r[coorAJ + k];
                    // console.log(`${k},${j}\t${tempRe},${tempIm}`);
                    for (let i = 1; i <= m; i++) {
                        // B(I, J) = B(I, J) + TEMP * B(I, K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        // console.log(`${k},${i}, ${j}\t${re},${im}`);
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                        //     console.log(`${k},${i}, ${j}\t${b.r[coorBJ + i]},${b.i[coorBJ + i]}`);
                    }
                }
            }//k
        }//j
    }//upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            if (nounit) {
                //  TEMP = TEMP * A(J, J)
                let re = tempRe * a.r[coorAJ + j] - tempIm * a.i[coorAJ + j];
                let im = tempRe * a.i[coorAJ + j] + tempIm * a.r[coorAJ + j];
                tempRe = re;
                tempIm = im;
            }
            for (let i = 1; i <= m; i++) {
                //   B(I,J) = TEMP*B(I,J)
                let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                b.r[coorBJ + i] = re;
                b.i[coorBJ + i] = im;
            }
            for (let k = j + 1; k <= n; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    tempRe = alpha.re * a.r[coorAJ + k] - alpha.im * a.i[coorAJ + k];
                    tempIm = alpha.re * a.i[coorAJ + k] + alpha.im * a.r[coorAJ + k];
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) + TEMP*B(I,K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }//k
        }
    }
}
