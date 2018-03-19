/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/
import { errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';
/*
*>
*> CTBMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.
*/

const { max, min } = Math;

export function ctbmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    a: Matrix2D,
    lda: number,
    x: FortranArr,
    incx: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    // faster then String.toLowerCase()
    const ul: 'u' | 'l' = String.fromCharCode(uplo.charCodeAt(0) | 0X20) as any;
    const tr: 'n' | 't' | 'c' = String.fromCharCode(trans.charCodeAt(0) | 0X20) as any;
    const dg: 'u' | 'n' = String.fromCharCode(diag.charCodeAt(0) | 0X20) as any;

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
    else if (k < 0) {
        info = 5;
    }
    else if (lda < (k + 1)) {
        info = 7;
    }
    else if (incx === 0) {
        info = 9;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctbmv', info));
    }

    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        if (ul === 'u') {
            let kplus1 = k + 1;
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    let L = kplus1 - j;
                    const coords = a.colOfEx(j);
                    for (let i = max(1, j - k); i <= j - 1; i++) {
                        x.r[ix] += tempRe * a.r[coords + L + i] - tempIm * a.i[coords + L + i];
                        ix += incx;
                    }
                    if (nounit) {
                        x.r[jx] += x.r[jx] * a.r[coords + kplus1] - x.i[jx] * a.i[coords + kplus1];
                        x.i[jx] += x.r[jx] * a.i[coords + kplus1] + x.i[jx] * a.r[coords + kplus1];
                    }
                    jx += incx;
                    if (j > k) kx += incx;
                }
            }//for
        }
        //ul !== ['u']
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                const coords = a.colOfEx(j);
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    let L = 1 - j;
                    for (let i = min(n, j + k); i >= j + 1; i--) {
                        x.r[ix] += tempRe * a.r[coords + L + i] - tempIm * a.i[coords + L + i];
                        x.i[ix] += tempIm * a.i[coords + L + i] + tempRe * a.r[coords + L + i];
                        ix -= incx;
                    }
                    if (nounit) {
                        x.r[jx] = x.r[jx] * a.r[coords + 1] - x.i[jx] * a.i[coords + 1];
                        x.i[jx] = x.r[jx] * a.i[coords + 1] + x.i[jx] * a.r[coords + 1];
                    }
                }
                jx -= incx;
                if (n - j >= k) kx -= incx;
            }
        }
    } else {
        if (ul === 'u') {
            let kplus1 = k + 1;
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];

                kx -= incx;
                let ix = kx - x.base;
                let L = kplus1 - j;
                const coords = a.colOfEx(j);
                const extrI = max(1, j - k); // evaluate once!
                if (noconj) {
                    if (nounit) {
                        tempRe = tempRe * a.r[coords + kplus1] - tempIm * a.i[coords + kplus1];
                        tempIm = tempRe * a.i[coords + kplus1] + tempIm * a.r[coords + kplus1];
                    }
                    for (let i = j - 1; j >= extrI; i--) {
                        tempRe += a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * a.i[ix];
                        tempIm += a.r[coords + L + i] * x.i[ix] + a.i[coords + L + i] * x.r[ix];
                        ix -= incx;
                    }
                }
                else {
                    if (nounit) {
                        // (a+ib)*(c-id) =(ac-iad+ibc+bd) = (ac+bd)+i(-ad+bc)
                        tempRe += (tempRe * a.r[coords + kplus1] + tempIm * a.i[coords + kplus1]);
                        tempIm += (-tempRe * a.i[coords + kplus1] + tempIm * a.r[coords + kplus1]);
                    }

                    for (let i = j - 1; i >= extrI; i--) {
                        // (a-ib)*(c+id) = ac+iad-ibc+bd = (ac+bd)+i(ad-bc)
                        tempRe += a.r[coords + L + i] * x.r[ix] + a.i[coords + L + i] * x.i[ix];
                        tempIm += a.r[coords + L + i] * x.i[ix] - a.i[coords + L + i] * x.r[ix];
                        ix -= incx;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                kx += incx;
                let ix = kx - x.base;
                const L = 1 - j;
                const coords = a.colOfEx(j);
                const extrI = max(n, j + k); // evaluate once!
                if (noconj) {
                    if (nounit) {
                        tempRe = tempRe * a.r[coords + 1] - tempIm * a.i[coords + 1];
                        tempIm = tempRe * a.i[coords + 1] + tempIm * a.r[coords + 1];
                    }
                    for (let i = j + 1; i <= extrI; i++) {
                        tempRe += a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * x.i[ix];
                        tempIm += a.r[coords + L + i] * x.i[ix] + a.i[coords + L + i] * x.r[ix];
                        ix += incx;
                    }
                }
                else {
                    if (nounit) {
                        //(a+ib)*(c-id) = (ac+bd) + i(-ad+bc)
                        tempRe = tempRe * a.r[coords + 1] + tempIm * a.i[coords + 1];
                        tempIm = -tempRe * a.i[coords + 1] + tempIm * a.r[coords + 1];
                    }
                    for (let i = j + 1; i <= extrI; i++) {
                        //(a-ib)*(c+id) = ac+iad-ibc+bd = (ac+bd)+i(ad-bc)
                        tempRe += a.r[coords + L + i] * x.r[ix] + a.i[coords + L + i] * x.i[ix];
                        tempIm += a.r[coords + L + i] * x.i[ix] - a.i[coords + L + i] * x.r[ix];
                        ix += incx;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
            }
        }
    }
}
