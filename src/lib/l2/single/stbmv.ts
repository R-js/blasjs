/*
 Jacob Bogers, 03/2018, jkfbogers@gmail.com
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

import { errWrongArg, FortranArr, Matrix } from '../../f_func';

//    STBMV  performs one of the matrix-vector operations
//    x := A*x,   or   x := A**T*x

/*
    where x is an n element vector and  A is an n by n unit, or non-unit,
    upper or lower triangular band matrix, with ( k + 1 ) diagonals.

    TRANS = 'N' or 'n'   x := A*x.
    TRANS = 'T' or 't'   x := A**T*x.
    TRANS = 'C' or 'c'   x := A**T*x.
*/

const { max } = Math;

export function stbmv(
    _uplo: 'U' | 'L',
    trans: 'N' | 'T' | 'C',
    diag: 'U' | 'N',
    n: number,
    k: number,
    A: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

    const ul = _uplo.toUpperCase()[0];
    const tr = trans.toUpperCase()[0];
    const dg = diag.toUpperCase()[0];

    let info = 0;

    if (ul !== 'U' && ul !== 'L') {
        info = 1;
    }
    else if (tr !== 'N' && tr !== 'T' && tr !== 'C') {
        info = 2;
    }
    else if (dg !== 'U' && dg !== 'N') {
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
        throw new Error(errWrongArg('stbmv', info));
    }

    if (n === 0) return;

    let kx = incx < 0 ? 1 : 1 - (n - 1) * incx;

    if (tr === 'N') {
        // Form  x := A*x.
        if (ul === 'U') {
            const kplus1 = k + 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    let L = kplus1 - j;
                    let coords = A.colOfEx(j);
                    for (let i = max(1, j - k); i <= j - 1; i++) {
                        x.r[ix - x.base] += temp * A.r[coords + L + i];
                        ix += incx;
                    }
                    if (dg === 'N') x.r[jx - x.base] *= A.r[coords + kplus1];
                }
                jx += incx;
                if (j > k) kx += incx;
            }
        }
    }//N
}
