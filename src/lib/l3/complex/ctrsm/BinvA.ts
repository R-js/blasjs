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

import { Complex, div_rxr, MatrixEComplex } from '../../../f_func';
//Form  B := alpha*inv( A )*B.

export function BinvA(
    nounit: boolean,
    upper: boolean,
    alphaIsOne: boolean,
    _alphaIsZero: boolean,
    _noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex,
): void {
    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = 1; k <= j - 1; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) - A(K,J)*B(I,K)
                        const re = a.r[coorAJ + k] * b.r[coorBK + i] - a.i[coorAJ + k] * b.i[coorBK + i];
                        const im = a.r[coorAJ + k] * b.i[coorBK + i] + a.i[coorAJ + k] * b.r[coorBK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                }
            }
            if (nounit) {
                // TEMP = ONE/A(J,J)
                //(a+ib)/(c+id)
                // re= (ac+bd)/(c*c+d*d)
                // im =(bc-ad)/(c*c+d*d)
                // (1+i0)/(c+id), a=1,b=0
                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                const { re: tempRe, im: tempIm } = div_rxr(1, 0, a.r[coorAJ + j], a.i[coorAJ + j]);
                /*
                                let _c = a.r[coorAJ + j];
                                let _d = a.i[coorAJ + j];
                                n = _c * _c + _d * _d;
                                let tempRe = _c / n;
                                let tempIm = _d / n;*/
                for (let i = 1; i <= m; i++) {
                    const re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                    const im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            } // nounit
        } //k
    } //upper
    else {
        for (let j = n; j >= 1; j--) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    /*const { re, im } = mul_cxr(
                        alpha,
                        b.r[coorBJ + i],
                        b.i[coorBJ + i]
                    );*/
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            //console.log(`j+1:${j + 1},n:${n}`);
            for (let k = j + 1; k <= n; k++) {
                const coorBK = b.colOfEx(k);
                // console.log(`j:${j},k:${k}`);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) - A(K,J)*B(I,K)
                        const re = a.r[coorAJ + k] * b.r[coorBK + i] - a.i[coorAJ + k] * b.i[coorBK + i];
                        const im = a.r[coorAJ + k] * b.i[coorBK + i] + a.i[coorAJ + k] * b.r[coorBK + i];
                        /*const { re, im } = mul_rxr(
                            a.r[coorAJ + k],
                            a.i[coorAJ + k],
                            b.r[coorBK + i],
                            b.i[coorBK + i]
                        );*/
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                }
            }
            if (nounit) {
                const { re: tempRe, im: tempIm } = div_rxr(1, 0, a.r[coorAJ + j], a.i[coorAJ + j]);
                // TEMP = ONE/A(J,J)
                //(a+ib)/(c+id)
                // re= (ac+bd)/(c*c+d*d)
                // im =(bc-ad)/(c*c+d*d)
                // (1+i0)/(c+id), a=1,b=0
                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                //let _c = a.r[coorAJ + j];
                //let _d = a.i[coorAJ + j];
                //n = _c * _c + _d * _d;
                //let tempRe = _c / n;
                //let tempIm = _d / n;
                for (let i = 1; i <= m; i++) {
                    const re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                    const im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            } // nounit
        }
    }
}
