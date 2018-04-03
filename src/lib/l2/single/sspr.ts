import { errWrongArg, FortranArr, lowerChar } from '../../f_func';

/*
  Jacob Bogers, 03/2008, jkfbogers@gmail.com
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.

SSPR    performs the symmetric rank 1 operation
    A := alpha*x*x**T + A,
*/

export function sspr(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    ap: FortranArr): void {

    // test parameters

    let info = 0;
    let ul = lowerChar(uplo);

    if (!'ul'.includes(uplo)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('sspr', info));
    }

    if (n === 0 || alpha === 0) return;

    // corrected the comparison was incx <= 0
    const kx = (incx < 0) ? 1 - (n - 1) * incx : 1;

    let kk = 1;
    let jx = kx;

    if (ul === 'u') {
        //  Form  A  when upper triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = kx;
                for (let k = kk; k <= kk + j - 1; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
            kk += j;
        }

    }
    else {
        //  Form  A  when lower triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0) {
                let temp = alpha * x.r[jx - x.base];
                let ix = jx;
                for (let k = kk; k <= kk + n - j; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp;
                    ix += incx;
                }
            }
            jx += incx;
            kk += (n - j + 1);
        }
    }
}


