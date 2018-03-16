import { Complex, errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';

const { max } = Math;

/*  -- Jacob Bogers, 03/2018/, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

/*
CGERC  performs the rank 1 operation
*>
*>    A := alpha*x*y**H + A,
*>
*> where alpha is a scalar, x is an m element vector, y is an n element
*> vector and A is an m by n matrix.
*/

export function cgerc(
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

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const betaIsZero = beta.re === 0 && beta.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;

    //stripp mine
    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;


    let info = 0;

    if (m < 0) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (incy === 0) {
        info = 7;
    }
    else if (lda < max(1, m)) {
        info = 9;
    }

    if (info !== 0) {
        throw new Error(errWrongArg('cgerc', info));
    }

    // Quick return if possible.

    if (m === 0 || n === 0 || alphaIsZero) return;

    //Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through A.

    let jy = incy > 0 ? 1 : 1 - (n - 1) * incy;
    let kx = incx > 0 ? 1 : 1 - (m - 1) * incx;

    for (let j = 1; j <= n; j++) {
        if (y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0) {
            //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
            //CONJ( (a + bi) )(c+di)= (a*c+b*d), i(a*d-b*c)
            let tempRe = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
            let tempIm = AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];

            let ix = kx;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                a.r[coords] += x.r[ix - x.base] * tempRe - x.i[ix - x.base] * tempIm;
                a.i[coords] += x.r[ix - x.base] * tempIm + x.i[ix - x.base] * tempRe;
                ix += incx;
            }
            jy += incy;
        }
    }
}
