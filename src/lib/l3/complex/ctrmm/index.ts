/*
*>  -- Jacob Bogers, 03/2018, Javascript Port, jkfbogers@gmail.com
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*/

import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../../f_func';

import { AB } from './AB';
import { AtranB, AtranB as AconjB } from './AtranB';
import { BA } from './BA';
import { BtranA, BtranA as BconjA } from './BtranA';

const { max } = Math;

/*
*>
*> CTRMM  performs one of the matrix-matrix operations
*>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*/
export function ctrmm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    transA: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    m: number,
    n: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number): void {

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }

    const si = lowerChar(side);
    const ul = lowerChar(uplo);
    const trA = lowerChar(transA);
    const di = lowerChar(diag);

    const nrowA = si === 'l' ? m : n;

    const noconj = trA === 't';
    const nounit = diag === 'n';
    const upper = ul === 'u';

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;

    let info = 0;

    if (!'lr'.includes(si)) {
        info = 1;
    }
    else if (!'ul'.includes(ul)) {
        info = 2;
    }
    else if (!'ntc'.includes(trA)) {
        info = 3;
    }
    else if ('un'.includes(di)) {
        info = 4;
    }
    else if (m < 0) {
        info = 5;
    }
    else if (n < 0) {
        info = 6;
    }
    else if (lda < max(1, nrowA)) {
        info = 9;
    }
    else if (ldb < max(1, m)) {
        info = 11;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctrmm', info));
    }

    //     Quick return if possible.

    if (m === 0 || n === 0) return;

    //     And when  alpha.eq.zero.

    if (alphaIsZero) {
        for (let j = 1; j <= n; j++) {
            b.setCol(j, 1, m, 0);
        }
        return;
    }

    //     Start the operations.
    if (si === 'l' && trA === 'n') {
        //Form  B := alpha*A*B.
        return AB(nounit, upper, n, m, a, b, alpha);
    }
    if (si === 'l' && trA === 't') {
        //Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
        return AtranB(nounit, upper, noconj, n, m, a, b, alpha)
    }
    if (si === 'l' && trA === 'c') {
        //B := alpha*A**H*B.
        return AconjB(nounit, upper, noconj, n, m, a, b, alpha)
    }
    if (si === 'r' && trA === 'n') {
        // Form  B := alpha*B*A.
        return BA(nounit, upper, n, m, a, b, alpha);
    }
    if (si === 'r' && trA === 't') {
        //Form  B := alpha*B*A**T 
        return BtranA(nounit, upper, noconj, n, m, a, b, alpha);
    }
    if (si === 'r' && trA === 'c') {
        // B := alpha*B*A**H.
        return BconjA(nounit, upper, noconj, n, m, a, b, alpha);
    }
}
