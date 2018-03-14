import { Complex, errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';

/*
*>  -- Jacob Bogers, 03/2018, JS port
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

//helpers

const { max } = Math;

/*
*>
*> CGEMV performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
*>
*>    y := alpha*A**H*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*/

/*
*>    TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>    TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>    TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
*/


export function cgemv(
    trans: 'n' | 't' | 'c',
    m: number,
    n: number,
    alpha: Complex,
    a: Matrix2D,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number): void {

    //checks

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }
    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    // dont use String.toUpperCase()[0] because slow
    const tr = String.fromCharCode(trans.charCodeAt(0) | 0x20);


    const betaIsZero = beta.re === 0 && beta.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;

    //stripp mine
    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;


    let info = 0;

    if (!(tr === 'n' || tr === 't' || tr === 'c')) {
        info = 1;
    }
    else if (m < 0) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (lda < max(1, m)) {
        info = 6;
    }
    else if (incx === 0) {
        info = 8;
    }
    else if (incy === 0) {
        info = 11;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('cgemv', info));
    }

    //*     Quick return if possible.

    if (m === 0 || n === 0 || (alphaIsZero && betaIsOne)) {
        return;
    }

    const noconj = tr === 't';

    const lenx = tr === 'n' ? n : m;
    const leny = tr === 'n' ? m : n;

    let kx = incx > 0 ? 1 : 1 - (lenx - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (leny - 1) * incy;


    /*     Start the operations. In this version the elements of A are
     *     accessed sequentially with one pass through A.
     *
     *     First form  y := beta*y.
     */
    if (betaIsOne) {
        let iy = ky;
        if (betaIsZero && incy === 1) {
            y.r.fill(0);
            y.i.fill(0);
        }
        else {
            for (let i = 1; i <= leny; i++) {
                y.r[iy - y.base] = betaIsZero ? 0 : BetaRe * y.r[iy - y.base] - BetaIm * y.i[iy - y.base];
                y.i[iy - y.base] = betaIsZero ? 0 : BetaRe * y.i[iy - y.base] + BetaIm * y.r[iy - y.base];
                iy += incy;
            }
        }
    }
    if (alphaIsZero) return;
    if (tr === 'n') {
        let jx = kx;
        for (let j = 1; j <= n; j++) {
            let tempRe = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let tempIm = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];
            let iy = ky;
            const coords = a.colOf(j) - a.rowBase;
            for (let i = 1; i <= m; i++) {
                y.r[iy - y.base] += tempRe * a.r[coords + i] - tempIm * a.i[coords + i];
                y.i[iy - y.base] += tempRe * a.i[coords + i] + tempIm * a.r[coords + i];
                iy += incy;
            }
            jx += incx;
        }
    }
    else {
        // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
        let jy = ky;

        for (let j = 1; j <= n; j++) {
            let tempRe = 0;
            let tempIm = 0;
            let ix = kx;
            const coords = a.colOf(j) - a.rowBase;
            if (noconj) {
                for (let i = 1; i <= m; i++) {
                    tempRe += a.r[coords + i] * x.r[ix - x.base] - a.i[coords + i] * x.i[ix - x.base];
                    tempIm += a.r[coords + i] * x.i[ix - x.base] + a.i[coords + i] * x.r[ix - x.base];
                    ix += incx;
                }
            } else {
                for (let i = 1; i <= m; i++) {
                    tempRe += a.r[coords + i] * x.r[ix - x.base] + a.i[coords + i] * x.i[ix - x.base];
                    tempIm += a.r[coords + i] * x.i[ix - x.base] - a.i[coords + i] * x.r[ix - x.base];
                    ix += incx;
                }
            }
            y.r[jy - y.base] += AlphaRe * tempRe - AlphaRe * tempIm;
            y.i[jy - y.base] += AlphaRe * tempIm + AlphaRe * tempRe;
            jy += incy;
        }
    }
}
