
/*
*>  -- Jacob Bogers, Javascript Port, 03/2018, jkfbogers@gmail
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*/

import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../../f_func';
import { BinvA } from './BinvA';
import { BinvTranConjA } from './BinvTranConjA';
import { invAB } from './invAB';
import { invTranConjAB } from './invTranConjAB';

/*
*>
*> CTRSM  solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*>
*> The matrix X is overwritten on B.
*/

const { max } = Math;

export function ctrsm(
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

    const alphaIsZero = (alpha.re === 0 && alpha.im === 0);

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
    else if (!'un'.includes(di)) {
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
        throw new Error(errWrongArg('ctrsm', info));
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

    // start operations

    if (si === 'l' && trA === 'n') {
        //Form  B := alpha*inv( A )*B.
        //invA*B
        return invAB(nounit, upper, n, m, a, b, alpha);
    }

    if (si === 'l' && trA !== 'n') {
        //Form  B := alpha*inv( A**T )*B
        //or    B := alpha*inv( A**H )*B.
        invTranConjAB(nounit, upper, noconj, n, m, a, b, alpha);
    }

    if (si === 'r' && trA === 'n') {
        //   Form  B := alpha*B*inv( A ).
        //BinvA
        return BinvA(nounit, upper, n, m, a, b, alpha);
    }

    if (si === 'r' && trA !== 'n') {
        //Form  B := alpha*B*inv( A**T )
        // or    B := alpha*B*inv( A**H ).
        //BinvTranConjA
        return BinvTranConjA(nounit, upper, noconj, n, m, a, b, alpha);
    }

    throw new Error('unreachable code');

}
