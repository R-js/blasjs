import { errWrongArg, FortranArr } from '../../f_func';

/*
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

/*
 SSPMV  performs the matrix-vector operation
    y := alpha*A*x + beta*y,
    where alpha and beta are scalars, x and y are n element vectors and
    A is an n by n symmetric matrix, supplied in packed form.
*/

export function sspmv(
    _uplo: string,
    n: number,
    alpha: number,
    ap: FortranArr, // a symmetric matrix in packed form
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number): void {

    const ul = _uplo.toUpperCase()[0];

    let info = 0;

    if (ul !== 'U' && ul !== 'L') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 6;
    }
    else if (incy === 0) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('SSPMV', info));
    }

    if (n === 0 || (alpha === 0 && beta === 1)) return;

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    //    First form  y := beta*y.

    if (beta !== 1) {
        if (beta === 0 && incy === 1) {
            y.r.fill(0);
        }
        else {
            let iy = ky;
            for (let i = 1; i <= n; i++) {
                y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
                iy += incy;
            }
        }
    }
    //
    if (alpha === 0) return;

    let kk = 1;

    if (ul === 'U') {
        // Form  y  when AP contains the upper triangle.

        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            let ix = kx;
            let iy = ky;
            for (let k = kk; k < kk + j - 2; k++) {
                y.r[iy - y.base] += temp1 * ap.r[k - ap.base];
                temp2 += ap.r[k - ap.base] * x.r[ix - x.base];
                ix += incx;
                iy += incy;
            }
            y.r[jy - y.base] += temp1 * ap.r[kk + j - 1 - ap.base] + alpha * temp2;
            jx += incx;
            jy += incy;
            kk += j;
        }

    }
    else {
        // Form  y  when AP contains the lower triangle.
        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp1 = alpha * x.r[jx - x.base];
            let temp2 = 0;
            y.r[jy - y.base] = temp1 * ap.r[kk - ap.base];
            let ix = jx;
            let iy = jy;
            for (let k = kk + 1; k < kk + n - j; k++) {
                ix += incx;
                iy += incy;
                y.r[iy - y.base] += temp1 * ap.r[k - ap.base];
                temp2 += ap.r[k - ap.base] * x.r[ix - x.base];
            }
            y.r[jy - y.base] += alpha * temp2;
            jx += incx;
            jy += incy;
            kk += (n - j + 1);
        }
    }
}

