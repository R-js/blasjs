import { errWrongArg, FortranArr, Matrix } from '../../f_func';

/*

Written on 03/2008 by 
Jacob Bogers, jkfbogers@gmail.com

 -- Written on 22-October-1986.
    Jack Dongarra, Argonne National Lab.
    Jeremy Du Croz, Nag Central Office.
    Sven Hammarling, Nag Central Office.
    Richard Hanson, Sandia National Labs.
*/

/*

SSBMV  performs the matrix-vector  operation

    y := alpha*A*x + beta*y,

Where alpha and beta are scalars, x and y are n element vectors and
A is an n by n symmetric band matrix, with k super-diagonals.
*/

const { max, min } = Math;

export function ssbmv(
    _uplo: string,
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number
): void {

    // test input params

    const ul = _uplo.toLocaleUpperCase()[0];

    let info = 0;
    if (ul !== 'U' && ul !== 'L') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (k < 0) {
        info = 3;
    }
    else if (lda < (k + 1)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    else if (incy === 0) {
        info = 11;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('ssbmv', info));
    }

    if (n === 0 || (alpha === 0 && beta === 1)) return;

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incx;

    // First form  y := beta*y.

    if (beta !== 1) {
        let iy = ky;
        //efficiency
        if (beta === 0 && incy === 1) {
            y.r.fill(0);
        }
        else {
            for (let i = 1; i <= n; i++) {
                y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
                iy += incy;
            }
        }
    }

    if (alpha === 0) return;
    if (ul === 'U') {

        let kplus1 = k + 1;
        let jx = kx;
        let jy = ky;

        for (let j = 1; j <= n; j++) {
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            let ix = kx;
            let iy = ky;
            let l = kplus1 - j;
            let coords = a.colOfEx(j);
            for (let i = max(1, j - k); i <= j - 1; i++) {
                y.r[iy - y.base] = temp1 * a.r[coords + l + i];
                temp2 = temp2 + a.r[coords + l + i] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1 * a.r[coords + j] + alpha * temp2;
            jx += incx;
            jy += incy;
            if (j > k) {
                kx += incx;
                ky += incy;
            }
        }
    }
    // ul !=== 'U'
    // Form  y  when lower triangle of A is stored.
    else {
        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            y.r[jy - y.base] += temp1 * a.r[1 + coorAJ];
            let l = 1 - j;
            let ix = jx;
            let iy = jy;
            for (let i = j + 1; i <= min(n, j + k); i++) {
                ix += incx;
                iy += incy;
                y.r[iy - y.base] += temp1 * a.r[coorAJ + l + i];
                temp2 += a.r[coorAJ + l + i] * x.r[ix - x.base];
            }
            y.r[jy - y.base] += alpha * temp2;
            jx += incx;
            jy += incy;
        }
    }
} 
