import { errMissingIm, errWrongArg, FortranArr, Matrix } from '../../f_func';

/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

//*        Form  x := A**T*x  or  x := A**H*x.

const { max } = Math;

export function ctrmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    a: Matrix,
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
    else if (lda < max(1, n)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;
    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;


    if (trans === 'n') {
        //  Form  x := A*x.
        if (ul === 'u') {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = 1; i <= j - 1; i++) {
                        x.r[ix] += tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                        x.i[ix] += tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                        ix += incx;
                    }
                    if (nounit) {
                        const tr = x.r[jx] * a.r[coords + j] - x.i[jx] * a.i[coords + j];
                        const ti = x.r[jx] * a.i[coords + j] + x.i[jx] * a.r[coords + j];
                        x.r[jx] = tr;
                        x.i[jx] = ti;
                    }
                }
                jx += incx;
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = kx - x.base;
                    const coords = a.colOfEx(j);
                    for (let i = n; i >= j + 1; i--) {
                        x.r[ix] += tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                        x.i[ix] += tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                        ix -= incx;
                    }
                    if (nounit) {
                        const tr = x.r[jx] * a.r[coords + j] - x.i[jx] * a.i[coords + j];
                        const ti = x.r[jx] * a.i[coords + j] + x.i[jx] * a.r[coords + j];
                        x.r[jx] = tr;
                        x.i[jx] = ti;
                    }
                }
                jx -= incx;
            }
        }
    }
    else {
        //       Form  x := A**T*x  or  x := A**H*x.

        if (ul === 'u') {
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (noconj) {
                    if (nounit) {
                        const tr = tempRe * a.r[coords + j] - tempIm * a.i[coords + j];
                        const ti = tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                        tempRe = tr;
                        tempIm = ti;
                    }
                    for (let i = j - 1; i >= 1; i--) {
                        ix -= incx;
                        tempRe += a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                        tempIm += a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
                    }
                } else {

                    if (nounit) {
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        const tr = tempRe * a.r[coords + j] + tempIm * a.i[coords + j];
                        const ti = -tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                        tempRe = tr;
                        tempIm = ti;
                    }
                    for (let i = j - 1; i >= 1; i--) {
                        ix -= incx;
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe += a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                        tempIm += a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = jx;
                const coords = a.colOfEx(j);
                if (noconj) {
                    if (nounit) {
                        //(a+ib)*(c+id) = (ac-bd)+i(ad+bc)
                        const tr = tempRe * a.r[coords + j] - tempIm * a.i[coords + j];
                        const ti = tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                        tempRe = tr;
                        tempIm = ti;

                    }
                    for (let i = j + 1; i <= n; i++) {
                        ix += incx;
                        tempRe += a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                        tempIm += a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
                    }
                } else {
                    if (nounit) {
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        const tr = tempRe * a.r[coords + j] + tempIm * a.i[coords + j];
                        const ti = -tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                        tempRe = tr;
                        tempIm = ti;
                    }
                    for (let i = j + 1; i <= n; i++) {
                        ix += incx;
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        tempRe += a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                        tempIm += -a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
            }
        }
    }
}
