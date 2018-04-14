
/*
*>  -- Jacob Bogers, 03/2018,
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> CTPMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix, supplied in packed form.
*/

import { errMissingIm, errWrongArg, FortranArr } from '../../../f_func';


export function ctpmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    ap: FortranArr,
    x: FortranArr,
    incx: number): void {


    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('ap.i'));
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
    else if (incx === 0) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpmv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    //*     Start the operations. In this version the elements of AP are
    //*     accessed sequentially with one pass through AP.

    if (tr === 'n') {
        //  Form  x:= A*x.
        if (ul === 'u') {
            let kk = 1;
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    for (let k = kk; k <= kk + j - 2; k++) {
                        const apk = k - ap.base;
                        x.r[ix] += tempRe * ap.r[apk] - tempIm * ap.i[apk];
                        x.i[ix] += tempRe * ap.i[apk] + tempIm * ap.r[apk];
                        ix += incx;
                    }
                    if (nounit) {
                        const apk2 = kk + j - 1 - ap.base;
                        x.r[jx] = x.r[jx] * ap.r[apk2] - x.i[jx] * ap.i[apk2];
                        x.i[jx] = x.r[jx] * ap.i[apk2] + x.i[jx] * ap.r[apk2];
                    }
                }
                jx += incx;
                kk += j;
            }
        }
        else {
            let kk = n * (n + 1) / 2;
            let kx = (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                        const apk = k - ap.base;
                        x.r[ix] += tempRe * ap.r[apk] - tempIm * ap.i[apk];
                        x.i[ix] += tempRe * ap.i[apk] + tempIm * ap.r[apk];
                        ix -= incx;
                    }
                    if (nounit) {
                        const apk2 = kk - n + j - ap.base;
                        x.r[jx] = x.r[jx] * ap.r[apk2] - x.i[jx] * ap.i[apk2];
                        x.i[jx] = x.r[jx] * ap.i[apk2] + x.i[jx] * ap.r[apk2];
                    }
                }
                jx -= incx;
                kk -= (n - j + 1);
            }
        }
    }
    else {
        //Form  x := A**T*x  or  x := A**H*x.
        if (ul === 'u') {
            let kk = n * (n + 1) / 2;
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = jx;
                if (noconj) {
                    if (nounit) {
                        const apk = kk - ap.base;
                        const tr = tempRe * ap.r[apk] - tempIm * ap.i[apk];
                        const ti = tempRe * ap.i[apk] + tempIm * ap.r[apk];
                        tempRe = tr;
                        tempIm = ti;

                    }
                    for (let k = kk - 1; k >= kk - j + 1; k--) {
                        ix -= incx;
                        const apk = k - ap.base;
                        tempRe += ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm += ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                    }
                }
                else {
                    if (nounit) {
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        const apk2 = kk - ap.base;
                        const tr = tempRe * ap.r[apk2] - tempIm * ap.i[apk2];
                        const ti = tempRe * ap.i[apk2] + tempIm * ap.i[apk2];
                        tempRe = tr;
                        tempIm = ti;
                    }
                    for (let k = kk - 1; k >= kk - j + 1; k--) {
                        ix -= incx;
                        const apk = k - ap.base;
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe += ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm += ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
                kk -= j;
            }
        }
        else {
            let kk = 1;
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = jx;
                if (noconj) {
                    if (nounit) {
                        const apkk = kk - ap.base;
                        const tr = tempRe * ap.r[apkk] - tempIm * ap.i[apkk];
                        const ti = tempRe * ap.i[apkk] + tempIm * ap.r[apkk];
                        tempRe = tr;
                        tempIm = ti;
                    }
                    for (let k = kk + 1; kk + n - j; k++) {
                        ix += incx;
                        const apk = k - ap.base;
                        tempRe = ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm = ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                    }
                }
                else {
                    if (nounit) {
                        const tr = tempRe;
                        const ti = tempIm;
                        const apkk = kk - ap.base;
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        tempRe = tr * ap.r[apkk] + ti * ap.i[apkk];
                        tempIm = -tr * ap.i[apkk] + ti * ap.r[apkk];
                    }
                    for (let k = kk + 1; k <= kk + n - j; k++) {
                        ix += incx;
                        //(a-ib)*(c+id) =(ac+bd)+i(ad-bc)
                        const apk = k - ap.base;
                        tempRe += ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm += ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
                kk += n - j + 1;
            }
        }
    }
}

