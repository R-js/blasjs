/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

//*        Form  x := A**T*x  or  x := A**H*x.

import {
    errMissingIm,
    errWrongArg,
    FortranArr,
    FortranArrEComplex,
    lowerChar,
    Matrix,
    MatrixEComplex,
} from '../../../f_func';

import { normalLower } from './normal-lower';
import { normalUpper } from './normal-upper';
import { transLower } from './trans-lower';
import { transUpper } from './trans-upper';

const { max } = Math;

export function ctrmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

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
    else if (lda < max(1, n)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    let proc: (
        kx: number,
        noconj: boolean,
        nounit: boolean,
        x: FortranArrEComplex,
        incx: number,
        a: MatrixEComplex,
        n: number) => void;


    if (trans === 'n') {
        if (uplo === 'u') {
            proc = normalUpper;
        }
        else {
            proc = normalLower;
        }
    }
    else {
        if (uplo === 'u') {
            proc = transUpper;
        }
        else {
            proc = transLower;
        }
    };

    return proc(
        kx,
        noconj,
        nounit,
        <FortranArrEComplex>x,
        incx,
        <MatrixEComplex>a,
        n);

}
