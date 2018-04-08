import { errWrongArg, FortranArr, lowerChar } from '../../f_func';

/*
  Jacob Bogers, jkfbogers@gmail.com, 03/2008
  -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

export function stpsv(
    uplo: 'u' | 'l',
    trans: 't' | 'n' | 'c',
    diag: 'u' | 'n',
    n: number,
    ap: FortranArr,
    x: FortranArr,
    incx: number
): void {

    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0

    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (!'un'.includes(dg)) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (incx === 0) {
        info = 7;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('stpsv', info));
    }

    //  Quick return if possible.

    if (n === 0) return;

    const nounit = dg === 'n';

    //Set up the start point in X if the increment is not unity. This
    // will be  ( N - 1 )*INCX  too small for descending loops.

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    //Start the operations. In this version the elements of AP are
    //accessed sequentially with one pass through AP.

    if (tr === 'n') {
        //  Form  x := inv( A )*x.
        if (ul === 'u') {
            let kk = (n * (n + 1)) / 2;
            let jx = kx + (n - 1) * incx;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx - x.base] !== 0) {
                    if (nounit) x.r[jx - x.base] /= ap.r[kk - ap.base];
                    const temp = x.r[jx - x.base];
                    let ix = jx;
                    for (let k = kk - 1; k >= kk - j + 1; k--) {
                        ix -= incx;
                        x.r[ix - x.base] -= temp * ap.r[k - ap.base];
                    }
                }
                jx -= incx;
                kk -= j;
            }
        }
        else {
            let kk = 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                if (x.r[jx - x.base] !== 0) {
                    if (nounit) x.r[jx - x.base] /= ap.r[kk - ap.base];
                    const temp = x.r[jx - x.base];
                    let ix = jx;
                    for (let k = kk + 1; k <= kk + n - j; k++) {
                        ix += incx;
                        x.r[ix - x.base] -= temp * ap.r[k - ap.base];
                    }
                }
                jx += incx;
                kk += (n - j + 1);
            }
        }
    }
    else {
        //  Form  x := inv( A**T )*x.
        if (ul === 'u') {
            let kk = 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                let temp = x.r[jx - x.base];
                let ix = kx;
                for (let k = kk; k <= kk + j - 2; k++) {
                    temp -= ap.r[k - ap.base] * x.r[ix - x.base];
                    ix += incx;
                }
                if (nounit) temp /= ap.r[kk + j - 1 - ap.base];
                x.r[jx - x.base] = temp;
                jx += incx;
                kk += j;
            }
        }
        else {
            let kk = n * (n + 1) / 2;
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx - x.base];
                let ix = kx;
                for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                    temp -= ap.r[k - ap.base] * x.r[ix - x.base];
                    ix -= incx;
                }
                if (nounit) temp /= ap.r[kk - n + j - ap.base];
                x.r[jx - x.base] = temp;
                jx -= incx;
                kk -= (n - j + 1);
            }//for
        }//if
    }//if
}//proc
