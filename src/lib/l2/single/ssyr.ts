import { errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

/*
  Jacob Bogers, JavaScript Port, 03/2018,  jkfbogers@gmail.com

  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

/* 
    SSYR   performs the symmetric rank 1 operation
    A := alpha*x*x**T + A,
*/

const { max } = Math;

export function ssyr(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    a: Matrix,
    lda: number
): void {

    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(uplo)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (lda < max(1, n)) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssyr', info));
    }

    // Quick return if possible

    if (n === 0 || alpha === 0) return;

    //Original had .LE.0 but the 0 case is already handled
    const kx = (incx < 0) ? 1 - (n - 1) * incx : 1;
    //let ky = (incy < 0) ? 1 - (n - 1) * incy : 1;

    let jx = kx;

    if (ul === 'u') {
        //Form  A  when A is stored in upper triangle. 
        for (let j = 1; j <= n; j++) {

            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = kx;
                const coords = a.colOfEx(j);
                for (let i = 1; i <= j; i++) {
                    a.r[coords + i] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
        }
    }
    else {
        //  Form  A  when A is stored in lower triangle.
        jx = kx;
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = jx;
                const coorAJ = a.colOfEx(j);
                let i = j
                for (; i <= n; i++) {
                    const delta = x.r[ix - x.base] * temp;
                    //   console.log(`i:${i},j:${j}, A_${i}.${j}=${a.r[coorAJ + i]}, x_${i}.x_${j}*alpha=${delta},A_N=${a.r[coorAJ + i] + delta}`);

                    a.r[coorAJ + i] = a.r[coorAJ + i] + delta;
                    // console.log(a.r[coorAJ + i]);
                    ix += incx;
                }
            }
            jx += incx;
        }
    }
}









