/*
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

import { errWrongArg, FortranArr, Matrix2D } from '../../f_func';

const { max } = Math;

export function sger(
    m: number,
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    a: Matrix2D,
    lda: number): void {

    let err = 0;
    switch (true) {
        case (m < 0): err = 1; break;
        case (n < 0): err = 2; break;
        case (incx === 0): err = 5; break;
        case (incy === 0): err = 7; break;
        case (lda < max(1, m)): err = 9; break;
        default:
            err = 0;
    }

    if (err) {
        throw new Error(errWrongArg('sger', err));
    }

    if (m === 0 || n === 0 || alpha === 0) return;

    let jy = incy < 0 ? 1 - (n - 1) * incy : 1;
    let kx = incx < 0 ? 1 - (m - 1) * incx : 1;
    for (let j = 1; j <= n; j++) {
        if (y.r[jy - y.base] !== 0) {
            let temp = alpha * y.r[jy - y.base];
            let ix = kx;
            const coords = a.colOf(j);
            for (let i = 1; i <= m; i++) {
                a.r[coords + i - a.rowBase] += x.r[ix - x.base] * temp;
                ix += incx;
            }
        }
        jy += incy;
    }
}
