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
    MatrixEComplex
} from '../../../f_func';

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
    const alphaIsOne = (alpha.re === 1 && alpha.im === 0);

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
        return invAB(
            nounit,
            upper,
            alphaIsOne,
            alphaIsZero,
            noconj,
            n,
            m,
            <MatrixEComplex>a,
            <MatrixEComplex>b,
            alpha
        );
    }

    if (si === 'l' && trA !== 'n') {
        //Form  B := alpha*inv( A**T )*B
        //or    B := alpha*inv( A**H )*B.
        return invTranConjAB(
            nounit,
            upper,
            alphaIsOne,
            alphaIsZero,
            noconj,
            n,
            m,
            <MatrixEComplex>a,
            <MatrixEComplex>b,
            alpha
        );
    }

    if (si === 'r' && trA === 'n') {
        //   Form  B := alpha*B*inv( A ).
        //BinvA
        return BinvA(
            nounit,
            upper,
            alphaIsOne,
            alphaIsZero,
            noconj,
            n,
            m,
            <MatrixEComplex>a,
            <MatrixEComplex>b,
            alpha
        );
    }

    /* if (si === 'r' && trA !== 'n') {*/
    //Form  B := alpha*B*inv( A**T )
    // or    B := alpha*B*inv( A**H ).
    //BinvTranConjA
    return BinvTranConjA(
        nounit,
        upper,
        alphaIsOne,
        alphaIsZero,
        noconj,
        n,
        m,
        <MatrixEComplex>a,
        <MatrixEComplex>b,
        alpha
    );
    //}

    //throw new Error('unreachable code');

}
