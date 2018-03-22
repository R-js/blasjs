/*
  Jacob Bogers, JS Port, 03/2018

  --Written on 22 - October - 1986.
    Jack Dongarra, Argonne National Lab.
    Jeremy Du Croz, Nag Central Office.
    Sven Hammarling, Nag Central Office.
    Richard Hanson, Sandia National Labs.
*/
const { max } = Math;

import { errWrongArg, FortranArr, Matrix } from '../../f_func';

/** 
SSYR2  performs the symmetric rank 2 operation
A:= alpha * x * y ** T + alpha * y * x ** T + A,

    where alpha is a scalar, x and y are n element vectors and A is an n
by n symmetric matrix.

 @param a is REAL array, dimension(LDA, N)
    
Before entry with  UPLO = 'U' or 'u', the leading n by n
upper triangular part of the array A must contain the upper
triangular part of the symmetric matrix and the strictly
lower triangular part of A is not referenced.

On exit, the upper triangular part of the array A is overwritten by the
upper triangular part of the updated matrix.

Before entry with UPLO = 'L' or 'l', the leading n by n
lower triangular part of the array A must contain the lower
triangular part of the symmetric matrix and the strictly
upper triangular part of A is not referenced.On exit, the
lower triangular part of the array A is overwritten by the
lower triangular part of the updated matrix.

*/
export function ssyr2(
    _uplo: 'U' | 'L',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    A: Matrix,
    lda: number): void {

    const ul = _uplo.toUpperCase()[0];

    let info = 0;
    if (ul !== 'U' && ul !== 'L') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 3;
    }
    else if (incy === 0) {
        info = 5;
    }
    else if (lda < max(1, n)) { //n can be 0?
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('ssyr2', info));
    }

    //     Quick return if possible.

    if (n === 0 || alpha === 0) return;

    //     Set up the start points in X and Y if the increments are not both
    //     unity
    const kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    const ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    let jx = kx;
    let jy = ky;

    // Start the operations. In this version the elements of A are
    //    accessed sequentially with one pass through the triangular part
    //    of A.

    const xb = x.base;
    const yb = y.base;

    if (ul === 'U') {
        //Form  A  when A is stored in the upper triangle.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0 || y.r[jy - y.base] !== 0) {
                let temp1 = alpha * y.r[jy - y.base];
                let temp2 = alpha * x.r[jx - x.base];
                let ix = kx;
                let iy = ky;
                const coords = A.colOfEx(j);
                for (let i = 1; i <= j; i++) {
                    A.r[coords + i] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
        }//for
    }
    else {
        //Form  A  when A is stored in the lower triangle.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - xb] !== 0 || y.r[jy - yb]) {
                let temp1 = alpha * y.r[jy - yb];
                let temp2 = alpha * x.r[jx - xb];
                let ix = jx;
                let iy = jy;
                const coords = A.colOfEx(j);
                for (let i = j; j <= n; j++) {
                    A.r[coords + i] += x.r[ix - xb] * temp1 + y.r[iy - yb] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
        }
    }
}
