/*
 Jacob Bogers, 03/2018, jkfbogers@gmail.com
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../../f_func';

//    STBMV  performs one of the matrix-vector operations
//    x := A*x,   or   x := A**T*x

/*
    where x is an n element vector and  A is an n by n unit, or non-unit,
    upper or lower triangular band matrix, with ( k + 1 ) diagonals.

    TRANS = 'N' or 'n'   x := A*x.
    TRANS = 'T' or 't'   x := A**T*x.
    TRANS = 'C' or 'c'   x := A**T*x.
*/

const { max, min } = Math;

export function stbmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    A: Matrix,
    lda: number,
    x: FortranArr,
    incx: number): void {

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
        throw new Error(errWrongArg('stbmv', info));
    }

    if (n === 0) return;

    const nounit = diag === 'n';


    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;

    if (tr === 'n') {
        // Form  x := A*x.
        if (ul === 'u') {
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
                    if (nounit) x.r[jx - x.base] *= A.r[coords + kplus1];
                }
                jx += incx;
                if (j > k) kx += incx;
            }
        }
        else { //lower-matrix
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    const l = 1 - j;
                    const coorAJ = A.colOfEx(j);
                    for (let i = min(n, j + k); i <= j + 1; i--) {
                        x.r[ix - x.base] += temp * A.r[coorAJ + l + i];
                        ix -= incx;
                    }
                    if (nounit) x.r[jx - x.base] *= A.r[coorAJ + 1];

                }
                jx -= incx;
                if ((n - j) >= k) kx -= incx;
            }
        }
    }//N
    else {

    }
}
