/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

import {
    Complex,
    errMissingIm,
    errWrongArg,
    FortranArr,
    isZero,
    Matrix
} from '../../f_func';

const { max } = Math;

export function cgeru(
    m: number,
    n: number,
    alpha: Complex,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    a: Matrix,
    lda: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }
    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }


    let info = 0;
    if (m < 0) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (incy === 0) {
        info = 7;
    }
    else if (lda < max(1, m)) {
        info = 9;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('cgeru', info));
    }

    const { re: AlphaRe, im: AlphaIm } = alpha;

    const alphaIsZero = isZero(alpha);

    //Quick return if possible.

    if (m === 0 || n === 0 || alphaIsZero) return;

    let jy = incy > 0 ? 1 : 1 - (n - 1) * incy;
    let kx = incx > 0 ? 1 : 1 - (m - 1) * incx;

    for (let j = 1; j <= n; j++) {
        if (!(y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0)) {
            let tempRe = AlphaRe * y.r[jy - y.base] - AlphaIm * y.i[jy - y.base];
            let tempIm = AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];
            let ix = kx;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                a.r[coords + i] += x.r[ix - x.base] * tempRe - x.i[ix - x.base] * tempIm;
                a.i[coords + i] += x.r[ix - x.base] * tempIm + x.i[ix - x.base] * tempRe;
                ix += incx;
            }
        }
        jy += incy;
    }
}
