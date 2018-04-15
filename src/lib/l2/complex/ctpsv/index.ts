/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/
/*
*>
*> CTPSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix, supplied in packed form.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
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
/*
*>
*> CTPSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix, supplied in packed form.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*/



export function ctpsv(
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
        ap: FortranArrEComplex,
        n: number) => void;

    if (tr === 'n') {
        if (ul === 'u') {
            proc = normUpper;
            //console.log(`normHigher`);
        }
        else {
            //console.log(`normLower`);
            proc = normLower;
        }
    }
    else {
        if (ul === 'u') {
            //console.log(`transUpper`);
            proc = transUpper;
        }
        else {
            //console.log(`transLower`);
            proc = transLower;
        }
    }

    return proc(
        kx,
        noconj,
        nounit,
        <FortranArrEComplex>x,
        incx,
        <FortranArrEComplex>ap,
        n);

}
