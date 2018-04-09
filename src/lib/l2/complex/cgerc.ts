import {
    Complex,
    errMissingIm,
    errWrongArg,
    FortranArr,
    isZero,
    Matrix
} from '../../f_func';

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
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    a: Matrix,
    lda: number): void {

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const alphaIsZero = isZero(alpha);

    //stripp mine
    const { re: AlphaRe, im: AlphaIm } = alpha;
    //const { re: BetaRe, im: BetaIm } = beta;


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
        if (!(y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0)) {
            // TEMP = ALPHA*DCONJG(Y(JY))
            //(a + bi) )*(c-di)= (a*c+b*d)+ i(-a*d+b*c)
            let tempRe = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
            let tempIm = -AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];
            //console.log(`j:${j}, alpha*conj(y[j])=(${tempRe},${tempIm})`);
            let ix = kx;
            const coords = a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)

                const re = x.r[ix - x.base] * tempRe - x.i[ix - x.base] * tempIm;
                const im = x.r[ix - x.base] * tempIm + x.i[ix - x.base] * tempRe;

                //console.log(`i,j:(${i},${j}), alpha*x[i]*conj(y[j])=(${re},${im})`);

                a.r[coords + i] += re;
                a.i[coords + i] += im;

                ix += incx;
            }
        }
        jy += incy;
    }
}
