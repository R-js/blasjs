/** 
   Jacob Bogers, 03/2018
   
   -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

/*
 SSYMV  performs the matrix-vector  operation

    y := alpha*A*x + beta*y,

 where alpha and beta are scalars, x and y are n element vectors and
 A is an n by n symmetric matrix.
*/

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function ssymv(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number): void {

    const lu = lowerChar(uplo);

    //input check

    let info = 0;

    if (!'ul'.includes(lu)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (lda < max(1, n)) {
        info = 5;
    }
    else if (incx === 0) {
        info = 7;
    }
    else if (incy === 0) {
        info = 9
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssymv', info));
    }

    // Quick return if possible.
    if (n === 0 || (alpha === 0 && beta === 1)) return;

    // Set up the start points in  X  and  Y.

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.

    //     First form  y := beta*y.


    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through the triangular part
    //     of A.
    //     First form  y := beta*y.

    if (beta !== 1) {

        let iy = ky;
        for (let i = 1; i <= n; i++) {
            y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
            iy += incy;
        }

    }

    if (alpha === 0) return;

    let jx = kx;
    let jy = ky;

    if (lu === 'u') {
        // Form  y  when A is stored in upper triangle.
        // First form  y := beta*y.  
        for (let j = 1; j <= n; j++) {
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            let ix = kx;
            let iy = ky;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= j - 1; i++) {
                y.r[iy - y.base] += temp1 * a.r[coords + i];
                temp2 += a.r[coords + i] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1 * a.r[coords + j] + alpha * temp2;
            jx += incx;
            jy += incy;
        }
    }
    else {
        // Form  y,  when A is stored in lower triangle.
        for (let j = 1; j <= n; j++) {
            //
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            const coords = a.colOfEx(j);
            y.r[jy - y.base] += temp1 * a.r[coords + j];
            let ix = jx;
            let iy = jy;
            //
            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                iy += incy;
                y.r[iy - y.base] += temp1 * a.r[coords + i];
                temp2 += a.r[coords + i] * x.r[ix - x.base];
            }
            //
            y.r[jy - y.base] += alpha * temp2;
            jx += incx;
            jy += incy;
            //
        }
    }
}
