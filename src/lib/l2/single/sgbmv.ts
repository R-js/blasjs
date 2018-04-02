/*
Ported to Javascript 03/2018
     Jacob Bogers, jkfbogers@gmail.com

Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

const { min, max } = Math;

export function sgbmv(
    trans: string,
    m: number,
    n: number,
    kl: number,
    ku: number,
    alpha: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number
): void {
    // lowerCase it all in a fast way

    const tr = lowerChar(trans);

    let info = 0;

    if (!'ntc'.includes(tr)) {
        info = 1;
    }
    else if (m < 0) { info = 2; }
    else if (n < 0) { info = 3; }
    else if (kl < 0) { info = 4; }
    else if (ku < 0) { info = 5; }
    else if (lda < (kl + ku + 1)) { info = 8; }
    else if (incx === 0) { info = 10; }
    else if (incy === 0) { info = 13; }

    if (info !== 0) {// error
        throw new Error(errWrongArg('sgbmv', info));
    }

    if (m === 0 || n === 0 || (alpha === 0 && beta === 1)) return;

    const lenx = tr === 'n' ? n : m;
    const leny = tr === 'n' ? m : n;

    let kx = incx > 0 ? 1 : 1 - (lenx - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (leny - 1) * incy;

    //Start the operations. In this version the elements of A are
    //accessed sequentially with one pass through the band part of A.

    //First form  y = beta*y.
    if (beta !== 1) {
        let iy = ky;
        for (let i = 1; i <= leny; i++) {
            y.r[iy - y.base] = beta === 0 ? 0 : beta * y.r[iy - y.base];
            iy += incy;
        }
    }


    if (alpha === 0) return;

    const kup1 = ku + 1;

    if (tr === 'n') {
        // FORM: y := alpha*A*x + y.
        let jx = kx;

        for (let j = 1; j <= n; j++) {
            let temp = alpha * x.r[jx - x.base];
            let iy = ky;
            let k = kup1 - j;
            const coorAJ = a.colOfEx(j);
            for (let i = max(1, j - ku); i <= min(m, j + kl); i++) {
                y.r[iy - y.base] += temp * a.r[coorAJ + k + i];
                iy += incy;
            }
            jx += incx;
            if (j > ku) ky += incy;
        }
    }
    else {
        // Form  y := alpha*A**T*x + y.
        // A**T = transpose(A), aka $$ A^{t} $$
        let jy = ky;

        for (let j = 1; j <= n; j++) {
            let temp = 0;
            let ix = kx;
            let k = kup1 - j;
            const coorAJ = a.colOfEx(j);
            for (let i = max(1, j - ku); i <= min(m, j + kl); i++) {
                temp += a.r[coorAJ + k + i] * x.r[ix - x.base];
                ix += incx;
            }
            y.r[jy - y.base] += alpha * temp;
            jy += incy;
            if (j > ku) kx += incx;
        }
    }
}


