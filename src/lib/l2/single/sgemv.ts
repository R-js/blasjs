import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

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
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: number,
    y: FortranArr,
    incy: number): void {

    const tr = lowerChar(trans);

    let info = 0;
    if (!'ntc'.includes(tr)) {
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

    if (m === 0 || n === 0 || (alpha === 0 && beta === 1)) return;

    let lenx = tr === 'n' ? n : m;
    let leny = tr === 'n' ? m : n;

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
            y.r.fill(0, 1 - y.base, 1 - y.base + leny);
        }
        else {
            let iy = ky;
            for (let i = 1; i <= leny; i++) {
                y.r[iy - y.base] = beta === 0 ? 0 : y.r[iy - y.base] * beta;
                iy += incy;
            }
        }
    }
    if (alpha === 0) return;
    if (tr === 'n') {
        //Form  y := alpha*A*x + y.
        let jx = kx;
        for (let j = 1; j <= n; j++) {
            let temp = alpha * x.r[jx - x.base];
            let iy = ky;
            const coorAJ = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                y.r[iy - y.base] += temp * a.r[coorAJ + i];
                iy += incy;
            }
            jx += incx;
        }
    }
    else {
        //Form  y:= alpha * A ** T * x + y.
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let temp = 0;
            let ix = kx;
            const coorAJ = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                temp += a.r[coorAJ + i] * x.r[ix - x.base];
                ix += incx;
            }
            y.r[jy - y.base] += alpha * temp;
            jy += incy;
        }
    }

}

