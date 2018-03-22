import { errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';

/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/
/*
*> CTRSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*/

const { max } = Math;

export function ctrsv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
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
    else if (lda < max(1, n)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctrsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;
    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        if (ul === 'u') {
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                const coords = a.colOfEx(j);
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
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
                    let ix = jx;
                    for (let i = j - 1; i >= 1; i--) {
                        ix -= incx;
                        x.r[ix] -= tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                        x.i[ix] -= tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                    }
                }
                jx -= incx;
            }
        }
        else {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const coords = a.colOfEx(j);
                const isXZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!isXZero) {
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
                    let ix = jx;
                    for (let i = j + 1; i <= n; i++) {
                        ix += incx;
                        x.r[ix] -= tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                        x.i[ix] -= tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                    }
                }
                jx += incx;
            }
        }
    }
    else {
        if (ul === 'u') {
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let ix = kx - x.base;
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                const coords = a.colOfEx(j);
                if (noconj) {
                    for (let i = 1; i <= j - 1; i++) {
                        tempRe -= a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                        tempIm -= a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + j];
                        let xd = a.i[coords + j];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                else {
                    for (let i = 1; i <= j - 1; i++) {
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe -= a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                        tempIm -= a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {
                        // (a+ib)/(c-id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + j];
                        let xd = -a.i[coords + j];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
            }
        }
        else {
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                let ix = kx - x.base;
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                const coords = a.colOfEx(j);
                if (noconj) {
                    for (let i = n; i >= j + 1; i--) {
                        tempRe -= a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                        tempIm -= a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/(c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + j];
                        let xd = a.i[coords + j];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }

                } else {
                    for (let i = n; i >= j + 1; i--) {
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe -= a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                        tempIm -= a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/(c-id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = a.r[coords + j];
                        let xd = -a.i[coords + j];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
            }
        }
    }
}
