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

import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix, MatrixEComplex } from '../../../f_func';

import { AB } from './AB';
import { BA } from './BA';
import { BtranA, BtranA as BconjA } from './BtranA';
import { tranAB, tranAB as conjAB } from './tranAB';

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
    ldb: number,
): void {
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
    } else if (!'ul'.includes(ul)) {
        info = 2;
    } else if (!'ntc'.includes(trA)) {
        info = 3;
    } else if (!'un'.includes(di)) {
        info = 4;
    } else if (m < 0) {
        info = 5;
    } else if (n < 0) {
        info = 6;
    } else if (lda < max(1, nrowA)) {
        info = 9;
    } else if (ldb < max(1, m)) {
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

    let proc: (
        nounit: boolean,
        upper: boolean,
        noconj: boolean,
        n: number,
        m: number,
        a: MatrixEComplex,
        b: MatrixEComplex,
        alpha: Complex,
    ) => void;

    //     Start the operations.
    if (si === 'l' && trA === 'n') {
        proc = AB;
    } else if (si === 'l' && trA === 't') {
        //Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
        proc = tranAB;
    } else if (si === 'l' && trA === 'c') {
        //B := alpha*A**H*B.
        proc = conjAB;
    } else if (si === 'r' && trA === 'n') {
        // Form  B := alpha*B*A.
        proc = BA;
    } else if (si === 'r' && trA === 't') {
        //Form  B := alpha*B*A**T
        proc = BtranA;
    } else {
        /*if (si === 'r' && trA === 'c') {*/
        proc = BconjA;
    }

    proc(nounit, upper, noconj, n, m, <MatrixEComplex>a, <MatrixEComplex>b, alpha);
}
