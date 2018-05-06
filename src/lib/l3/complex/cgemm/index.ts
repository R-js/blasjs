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
    errMissingIm,
    errWrongArg,
    lowerChar,
    Matrix,
    MatrixEComplex,
    mul_cxr
} from '../../../f_func';

import { AB } from './AB';
import { AconjB } from './AconjB';
import { AtransB } from './AtransB';
import { conjAB } from './conjAB';
import { conjAconjB } from './conjAconjB';
import { conjAtransB } from './conjAtransB';
import { transAB } from './transAB';
import { transAconjB } from './transAconjB';
import { transAtransB } from './transAtransB';


const { max } = Math;

export function cgemm(
    transA: 'n' | 't' | 'c',
    transB: 'n' | 't' | 'c',
    m: number,
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: Complex,
    c: Matrix,
    ldc: number
): void {

    const trA = lowerChar(transA);
    const trB = lowerChar(transB);

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    //const notA = trA === 'n';
    //const notB = trB === 'n';

    const nrowA = trA === 'n' ? m : k;
    const nrowB = trB === 'n' ? k : m;

    //* Test the input parameters.



    let info = 0;
    if (!'ntc'.includes(trA)) {
        info = 1;
    }
    else if (!'ntc'.includes(trB)) {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (k < 0) {
        info = 5;
    }
    else if (lda < max(1, nrowA)) {
        info = 8;
    }
    else if (ldb < max(1, nrowB)) {
        info = 10;
    }
    else if (ldc < max(1, m)) {
        info = 13;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('cgemm', info));
    }

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;
    const betaIsZero = beta.re === 0 && beta.im === 0;

    // Quick return if possible.

    if (m === 0 || n === 0 || (
        (alphaIsZero || k === 0)
        && betaIsOne)) {
        return;
    }

    // And when  alpha.eq.zero.
    if (alphaIsZero) {
        //console.log('alpha is zero')
        if (betaIsZero) {//fast shortcut

            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    const { re, im } = mul_cxr(beta, c.r[coorCJ + i], c.i[coorCJ + i]);
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
        return;
    }

    let proc: (
        bIsZero: boolean,
        bIsOne: boolean,
        beta: Complex,
        alpha: Complex,
        a: MatrixEComplex,
        b: MatrixEComplex,
        c: MatrixEComplex,
        n: number,
        m: number,
        k: number) => void;

    //    Start the operations.
    if (trA === 'n' && trB === 'n') {
        //Form  C := alpha*A*B + beta*C.
        //console.log('AB');
        proc = AB;
    } else if (trA === 'n' && trB === 'c') {
        //Form  C := alpha*A*B**H + beta*C.
        proc = AconjB;
    } else if (trA === 'n' && trB === 't') {
        //Form  C := alpha*A*B**T + beta*C
        proc = AtransB;
    } else if (trA === 'c' && trB === 'n') {
        // Form  C := alpha*A**H*B + beta*C.   
        proc = conjAB;
    } else if (trA === 'c' && trB === 'c') {
        //  Form  C := alpha*A**H*B**H + beta*C. 
        proc = conjAconjB;
    } else if (trA === 'c' && trB === 't') {
        //  Form  C := alpha*A**H*B**T + beta*C 
        proc = conjAtransB;
    } else if (trA === 't' && trB === 'n') {
        //Form  C := alpha*A**T*B + beta*C  
        proc = transAB;
    } else if (trA === 't' && trB === 'c') {
        //Form  C := alpha*A**T*B**H + beta*C
        proc = transAconjB;
    } else /*if (trA === 't' && trB === 't')*/ {
        //  //  Form  C := alpha*A**T*B**T + beta*C
        proc = transAtransB;
    }
    // throw new Error('unreachable code');

    //console.log('name', proc.name);
    return proc(
        betaIsZero,
        betaIsOne,
        beta,
        alpha,
        <MatrixEComplex>a,
        <MatrixEComplex>b,
        <MatrixEComplex>c,
        n,
        m,
        k);
}

