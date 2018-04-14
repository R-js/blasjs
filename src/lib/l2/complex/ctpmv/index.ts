
/*
*>  -- Jacob Bogers, 03/2018,
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> CTPMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix, supplied in packed form.
*/

import {
    errMissingIm,
    errWrongArg,
    FortranArr,
    FortranArrEComplex,
    lowerChar
} from '../../../f_func';

import { normLower } from './norm-lower';
import { normUpper } from './norm-upper';
import { transLower } from './trans-lower';
import { transUpper } from './trans-upper';

// SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
export function ctpmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    ap: FortranArr,
    x: FortranArr,
    incx: number): void {


    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('ap.i'));
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
    else if (incx === 0) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpmv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    //*     Start the operations. In this version the elements of AP are
    //*     accessed sequentially with one pass through AP.
    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;

    if (tr === 'n') {
        if (ul === 'u') {
            return normUpper(
                kx,
                nounit,
                <FortranArrEComplex>x,
                incx,
                <FortranArrEComplex>ap,
                n
            );
        }
        return normLower(
            kx,
            nounit,
            <FortranArrEComplex>x,
            incx,
            <FortranArrEComplex>ap,
            n
        );
    }

    if (ul === 'u') {
        return transUpper(
            kx,
            noconj,
            nounit,
            <FortranArrEComplex>x,
            incx,
            <FortranArrEComplex>ap,
            n);
    }

    return transLower(
        kx,
        noconj,
        nounit,
        <FortranArrEComplex>x,
        incx,
        <FortranArrEComplex>ap,
        n);
}

