/*
*>  -- Jacob Bogers, 03/2018, JS Port
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.

*>
*> CGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*/

import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../../f_func';

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

    if (m === 0 || n === 0 || ((alphaIsZero || k === 0) && betaIsOne)) {
        return;
    }

    // And when  alpha.eq.zero.
    if (alphaIsZero) {
        if (betaIsZero) {//fast shortcut
            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        }
        else if (!betaIsOne) {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    c.i[coorCJ + i] = beta.re * c.i[coorCJ + i] + beta.re * c.r[coorCJ + i];
                }
            }
        }
        return;
    }

    //    Start the operations.
    switch (true) {
        case (trA === 'n' && trB === 'n'):
            //Form  C := alpha*A*B + beta*C.
            return AB(beta, alpha, a, b, c, n, m, k);
        case (trA === 'n' && trB === 'c'):
            //Form  C := alpha*A*B**H + beta*C.
            return AconjB(beta, alpha, a, b, c, n, m, k);
        case (trA === 'n' && trB === 't'):
            //Form  C := alpha*A*B**T + beta*C
            return AtransB(beta, alpha, a, b, c, n, m, k);


        case (trA === 'c' && trB === 'n'):
            // Form  C := alpha*A**H*B + beta*C.   
            return conjAB(beta, alpha, a, b, c, n, m, k);
        case (trA === 'c' && trB === 'c'):
            //  Form  C := alpha*A**H*B**H + beta*C. 
            return conjAconjB(beta, alpha, a, b, c, n, m, k);
        case (trA === 'c' && trB === 't'):
            //  Form  C := alpha*A**H*B**T + beta*C 
            return conjAtransB(beta, alpha, a, b, c, n, m, k);


        case (trA === 't' && trB === 'n'):
            //Form  C := alpha*A**T*B + beta*C  
            return transAB(beta, alpha, a, b, c, n, m, k);
        case (trA === 't' && trB === 'c'):
            //Form  C := alpha*A**T*B**H + beta*C
            return transAconjB(beta, alpha, a, b, c, n, m, k);
        case (trA === 't' && trB === 't'):
            //  //  Form  C := alpha*A**T*B**T + beta*C
            return transAtransB(beta, alpha, a, b, c, n, m, k);
        default:
            throw new Error('unreachable code');
    }

}
