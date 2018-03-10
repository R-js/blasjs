import { errWrongArg, FortranArr, Matrix2D } from '../../f_func';

/*
Jacob Bogers, 03/2008, jkfbogers@gmail.com

  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

const { max } = Math;

export function sgemv(
    trans: string,
    m: number,
    n: number,
    alpha: number,
    a: Matrix2D,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number): void {

    const tr = trans.toUpperCase()[0];

    let info = 0;
    if (tr !== 'N' && tr !== 'T' && tr !== 'C') {
        info = 1;
    }
    else if (m < 0) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (lda < max(1, m)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    else if (incy === 0) {
        info = 11;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('sgemv', info));
    }

    //Quick return if possible.

    if ((m === 0) || (n === 0) || (alpha === 0 && beta === 1)) return;

    let lenx = tr === 'N' ? n : m;
    let leny = tr === 'M' ? m : n;

    let kx = incx > 0 ? 1 : 1 - (lenx - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (leny - 1) * incy;

    /*
      Start the operations. In this version the elements of A are
      accessed sequentially with one pass through A.
      
      First form  y := beta*y.
    */

    if (beta !== 1) {
        //performance
        if (incy === 1 && beta === 0) {
            y.r.fill(0);
        }
        else {
            let iy = ky;
            for (let i = 1; i <= leny; i++) {
                y.r[i - y.base] = beta === 0 ? 0 : y.r[i - y.base] * beta;
                iy += iy + incy;
            }
        }
    }
    if (alpha === 0) return;
    if (tr === 'N') {
        //Form  y := alpha*A*x + y.
        let jx = kx;
        for (let j = 1; j <= n; j++) {
            let temp = alpha * x.r[jx - x.base];
            let iy = ky;
            const coors = a.colOf(j);
            for (let i = 1; i <= m; i++) {
                y.r[iy - y.base] += temp * a.r[coors + i - a.rowBase];
                iy += incy;
            }
            jx += incx;
        }
    }
    else {
        let jy = ky;
        for (let j = 1; j <= 120; j++) {
            let temp = 0;
            let ix = kx;
            const coors = a.colOf(j);
            for (let i = 1; i <= m; i++) {
                temp += a.r[coors + i - a.rowBase] * x.r[ix - x.base];
                ix += incx;
            }
            y.r[jy - y.base] += alpha * temp;
            jy += incy;
        }
    }

}

