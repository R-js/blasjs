/*
*>  -- Jacob Bogers, 03/2018
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> CTBSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )
*> diagonals.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*/

import {
    errMissingIm,
    errWrongArg,
    FortranArr,
    FortranArrEComplex,
    lowerChar,
    Matrix,
    MatrixEComplex

} from '../../../f_func';

import { normLower } from './norm-lower';
import { normUpper } from './norm-upper';
import { transLower } from './trans-lower';
import { transUpper } from './trans-upper';


export function ctbsv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number
): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }


    // faster then String.toLowerCase()
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (!'un'.includes(dg)) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (k < 0) {
        info = 5;
    }
    else if (lda < (k + 1)) {
        info = 7;
    }
    else if (incx === 0) {
        info = 9;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctbsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        if (ul === 'u') return normUpper(<FortranArrEComplex>x, incx, kx, <MatrixEComplex>a, nounit, n, k);
        return normLower(<FortranArrEComplex>x, incx, kx, <MatrixEComplex>a, nounit, n, k);
    }
    // a = transpose or transpose conjugate
    if (ul === 'u') return transUpper(noconj, <FortranArrEComplex>x, incx, kx, <MatrixEComplex>a, nounit, n, k);
    return transLower(noconj, <FortranArrEComplex>x, incx, kx, <MatrixEComplex>a, nounit, n, k);

}





