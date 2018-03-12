import { errWrongArg, FortranArr } from '../../f_func';

/*
Jacob Bogers, 03/2018, jkfbogers@gmail.com

--Written on 22 - October - 1986.
   Jack Dongarra, Argonne National Lab.
   Jeremy Du Croz, Nag Central Office.
   Sven Hammarling, Nag Central Office.
   Richard Hanson, Sandia National Labs.
*/

/*
 STPMV  performs one of the matrix-vector operations

    x := A*x,   or   x := A**T*x,

 where x is an n element vector and  A is an n by n unit, or non-unit,
 upper or lower triangular matrix, supplied in packed form.
*/

export function stpmv(
    _uplo: 'U' | 'L',
    trans: 'N' | 'T' | 'C',
    diag: 'U' | 'N',
    n: number,
    ap: FortranArr,
    x: FortranArr,
    incx: number
): void {

    const ul = _uplo.toUpperCase()[0];
    const tr = trans.toUpperCase()[0];
    const dg = trans.toUpperCase()[0];

    let info = 0;

    if (ul !== 'U' && ul !== 'L') {
        info = 1;
    }
    else if (tr !== 'N' && tr !== 'T' && tr !== 'C') {
        info = 2;
    }
    else if (dg !== 'U' && dg !== 'N') {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (incx === 0) {
        info = 7;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('stpmv', info));
    }

    //    Quick return if possible.

    if (n === 0) return;

    const nounit = dg === 'N';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    // Start the operations. In this version the elements of AP are
    // accessed sequentially with one pass through AP.

    if (tr === 'N') {
        if (ul === 'U') {
            let kk = 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    for (let k = kk; k <= kk + j - 2; k++) {
                        x.r[ix - x.base] += temp * ap.r[k - ap.base];
                        ix += incx;
                    }
                    if (nounit) x.r[jx - x.base] *= ap.r[kk + j - 1 - ap.base];
                }
                jx += incx;
                kk += j;
            }
        }
        else {
            let kk = n * (n + 1) / 2
            kx += (n - 1) * incx;
            let jx = kx;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx - x.base] !== 0) {
                    let temp = x.r[jx - x.base];
                    let ix = kx;
                    for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                        x.r[ix - x.base] += temp * ap.r[k - ap.base];
                        ix -= incx;
                    }
                    if (nounit) x.r[jx - x.base] *= ap.r[kk - n - j - ap.base];
                }
                jx -= incx;
                kk -= (n - j + 1);
            }
        }
    }
    else {
        // Form  x := A**T*x.
        if (ul === 'U') {
            let kk = (n * (n + 1)) / 2;
            let jx = kx + (n - 1) * incx;

            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx - x.base];
                let ix = jx;
                if (nounit) temp *= ap.r[kk - ap.base];
                for (let k = kk - 1; k >= kk - j + 1; k--) {
                    ix -= incx;
                    temp += ap.r[k - ap.base] * x.r[ix - x.base];
                }
                x.r[jx - x.base] = temp;
                jx -= incx;
                kk -= j;
            }
        }
        //ul !== u
        else {
            let kk = 1;
            let jx = kx;
            for (let j = 1; j <= n; j++) {
                let temp = x.r[jx - x.base];
                let ix = jx;
                if (nounit) temp *= ap.r[kk - ap.base];
                for (let k = kk + 1; k <= kk + n - j; k++) {
                    ix += incx;
                    temp += ap.r[k - ap.base] * x.r[ix - x.base];
                }
                x.r[jx - x.base] = temp;
                jx += incx;
                kk += (n - j - 1);
            }
        }
    }
}
