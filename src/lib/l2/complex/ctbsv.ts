/*
*>  -- Jacob Bogers, 03/2018
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
*>
*> CTBSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )
*> diagonals.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*/

import { errMissingIm, errWrongArg, FortranArr, Matrix } from '../../f_func';

const { max, min } = Math;

export function ctbsv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    k: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number
): void {

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
        throw new Error(errWrongArg('ctbsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (trans === 'n') {
        if (ul === 'u') {
            let kplus1 = k + 1;
            kx += (n - 1) * incx;
            let jx = kx - x.base;

            for (let j = n; j >= 1; j--) {
                const coords = a.colOfEx(j);
                kx -= incx;
                const isXZero = x.r[jx] === 0 && x.i[jx] === 0;
                const extrI = max(1, j - k);
                if (isXZero) {
                    let ix = kx - x.base;
                    let L = kplus1 - j;
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = x.r[jx];
                        let xb = x.i[jx];
                        let xc = a.r[coords + kplus1];
                        let xd = a.i[coords + kplus1];

                        let DOM = xc * xc + xd * xd;

                        x.r[jx] = (xa * xc + xb * xd) / DOM;
                        x.i[jx] = (-xa * xd + xb * xc) / DOM;
                    }
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];

                    for (let i = j - 1; j >= extrI; j--) {
                        x.r[ix] -= tempRe * a.r[coords + L + i] - tempIm * a.i[coords + L + i];
                        x.i[ix] -= tempRe * a.i[coords + L + i] - tempIm * a.i[coords + L + i];
                        ix -= incx;
                    }
                }
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                kx += incx;
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                const extrI = min(n, j + k);
                if (!xIsZero) {
                    let ix = kx - x.base;
                    let L = 1 - j;
                    const coords = a.colOfEx(j);
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = x.r[jx];
                        let xb = x.i[jx];
                        let xc = a.r[coords + j];
                        let xd = a.i[coords + j];

                        let DOM = xc * xc + xd * xd;

                        x.r[jx] = (xa * xc + xb * xd) / DOM;
                        x.i[jx] = (-xa * xd + xb * xc) / DOM;
                    }
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    for (let i = j + 1; i <= extrI; i++) {
                        x.r[ix] -= tempRe * a.r[coords + L + i] - tempIm * a.i[coords + L + i];
                        x.i[ix] -= tempRe * a.i[coords + L + i] + tempIm * a.r[coords + L + i];
                        ix += incx;
                    }
                }
                jx += incx;
            }
        }
    }
    else {
        // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
        if (ul === 'u') {
            let kplus1 = k + 1;
            let jx = kx - x.base;

            for (let j = 1; j <= n; j++) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = kx - x.base;
                let L = kplus1 - j;
                const extrI = max(1, j - k);
                const coords = a.colOfEx(j);
                if (noconj) {
                    for (let i = extrI; i <= j - 1; i++) {
                        tempRe -= a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * x.i[ix];
                        tempIm -= a.r[coords + L + i] * x.i[ix] + a.r[coords + L + i] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {

                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + kplus1];
                        let xd = a.i[coords + kplus1];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                } else {
                    for (let i = extrI; i <= j - 1; j++) {
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc);
                        tempRe -= a.r[coords + L + i] * x.r[ix] + a.i[coords + L + i] * x.i[ix];
                        tempIm -= a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + kplus1];
                        //conj(a)
                        let xd = -a.i[coords + kplus1];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
                if (j > k) {
                    kx += incx;
                }
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;

            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = kx - x.base;
                let L = 1 - j;

                //
                const extrI = min(n, j + k);
                const coords = a.colOfEx(j);
                if (noconj) {
                    for (let i = extrI; i >= j + 1; j--) {
                        tempRe -= a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * x.i[ix];
                        tempIm -= a.r[coords + L + i] * x.i[ix] + a.i[coords + L + i] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + 1];
                        let xd = a.i[coords + 1];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                else {
                    for (let i = extrI; i >= j + 1; i--) {
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe -= a.r[coords + L + i] * x.r[ix] - a.i[coords + L + i] * x.i[ix];
                        tempIm -= a.r[coords + L + i] * x.i[ix] + a.i[coords + L + i] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + 1];
                        //conj (a)
                        let xd = -a.i[coords + 1];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
                if ((n - j) >= k) {
                    kx -= incx;
                }
            }
        }
    }
}
