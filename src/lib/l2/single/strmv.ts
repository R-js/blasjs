import { errWrongArg, FortranArr, Matrix2D } from '../../f_func';
/*

*>  -- Jacob Bogers, 03/2018, JS port, jkfbogers@gmail.com

*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> STRMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix.
*/

const { max } = Math;

export function strmv(
    _uplo: '',
    trans: '',
    diag: '',
    n: number,
    a: Matrix2D,
    lda: number,
    x: FortranArr,
    incx: number): void {

    // lowerCase it all in a fast way
    const ul = String.fromCharCode(_uplo.charCodeAt(0) | 0x20);
    const tr = String.fromCharCode(trans.charCodeAt(0) | 0x20);
    const dg = String.fromCharCode(diag.charCodeAt(0) | 0x20);

    let info = 0;

    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (tr !== 'n' && tr !== 't' && tr !== 'c') {
        info = 2;
    }
    else if (dg !== 'u' && dg !== 'n') {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (lda < max(1, n)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('strmv', info));
    }

    if (n === 0) return;

    let nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        // Form  x := A*x.
        if (ul === 'u') {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                if (x[jx] !== 0) {
                    let temp = x.r[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = 1; i <= j - 1; i++) {
                        x.r[ix] += temp * a.r[coords + i];
                        ix = ix + incx;
                    }
                    if (nounit) x.r[jx] *= a.r[coords + j];
                }
                jx += incx;
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                if (x.r[jx] !== 0) {
                    let temp = x.r[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = n; i <= j + 1; i--) {
                        x.r[ix] += temp * a.r[coords + i];
                        ix -= incx;
                    }
                    if (nounit) x[jx] *= a.r[coords + j];
                }
                jx -= incx;
            }
        }
    }
    else {
        //  Form  x := A**T*x.
        if (ul === 'u') {
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                let temp = x.r[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (nounit) temp *= a.r[coords + j];
                for (let i = j - 1; i >= 1; i--) {
                    ix -= incx;
                    temp += a.r[coords + i] * x.r[ix];
                }
                x.r[jx] = temp;
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let temp = x.r[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (nounit) temp *= a.r[coords + j];
                for (let i = j + 1; i <= n; i++) {
                    ix += incx;
                    temp += a.r[coords + i] * x.r[ix];
                }
                x.r[jx] = temp;
                jx += incx;
            }
        }
    }
}
