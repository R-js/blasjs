/*
*>  -- Jacob Bogers, JS port , 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

import { Complex, errMissingIm, errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

// helpers

const { min, max } = Math;
/*
*>
*> CGBMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
*>
*>    y := alpha*A**H*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
*/

/*
*>    TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>    TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>    TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
*/
//ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
export function cgbmv(
    trans: 'n' | 't' | 'c',
    m: number,
    n: number,
    kl: number,
    ku: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number
): void {

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
    const tr = lowerChar(trans);

    const betaIsZero = beta.re === 0 && beta.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;

    //stripp mine
    const { re: AlphaRe, im: AlphaIm } = alpha;
    const { re: BetaRe, im: BetaIm } = beta;


    let info = 0;

    if (!'ntc'.includes(tr)) {
        info = 1;
    }
    else if (m < 0) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (kl < 0) {
        info = 4;
    }
    else if (ku < 0) {
        info = 5;
    }
    else if (lda < (kl + ku + 1)) {
        info = 8;
    }
    else if (incx === 0) {
        info = 10;
    }
    else if (incy === 0) {
        info = 13;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('cgbmv', info));
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



    //   Start the operations. In this version the elements of A are
    //   accessed sequentially with one pass through the band part of A.

    //   First form  y := beta*y.


    if (!betaIsOne) {
        let iy = ky;
        //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
        for (let i = 1; i <= leny; i++) {
            const re = betaIsZero ? 0 : BetaRe * y.r[iy - y.base] - BetaIm * y.i[iy - y.base];
            const im = betaIsZero ? 0 : BetaRe * y.i[iy - y.base] + BetaIm * y.r[iy - y.base];
            y.r[iy - y.base] = re;
            y.i[iy - y.base] = im;
            iy += incy;
        }

    }
    if (alphaIsZero) return;
    const kup1 = ku + 1;
    if (trans === 'n') { // not [t]ranspose or [c]onjugate
        //Form  y := alpha*A*x + y.
        let jx = kx;
        for (let j = 1; j <= n; j++) {

            //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
            let tempRe = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
            let tempIm = AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base];
            //console.log({ tempRe, tempIm });
            let iy = ky;
            let k = kup1 - j;
            const coorAJ = a.colOfEx(j);
            for (let i = max(1, j - ku); i <= min(m, j + kl); i++) {
                //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                //console.log(`i${i + k},j:${j}, a[i,j]=(${a.r[coorAJ + k + i]},${a.i[coorAJ + k + i]})`);
                const re = tempRe * a.r[coorAJ + k + i] - tempIm * a.i[coorAJ + k + i];
                const im = tempRe * a.i[coorAJ + k + i] + tempIm * a.r[coorAJ + k + i];
                //console.log(`i${i + k},j:${j}, a[i,j]=(${a.r[coorAJ + k + i]},${a.i[coorAJ + k + i]}), d=(${re},${im})`);
                //console.log({ re, im, ar: a.r[coords + k + i], ai: a.i[coords + k + i] });
                y.r[iy - y.base] += re;
                y.i[iy - y.base] += im;
                //console.log(`y(${iy})=(${y.r[iy - y.base]},${y.i[iy - y.base]})`);
                iy += incy;
            }
            jx += incx;
            if (j > ku) { ky += incy; }
        }
    }
    else {
        // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            let tempRe = 0;
            let tempIm = 0;
            let ix = kx;
            let k = kup1 - j;
            const coords = a.colOfEx(j);
            const istart = max(1, j - ku);
            const iend = min(m, j + kl);
            // NOTE: I am not going to merge the 2, keeping it readable this way
            if (noconj) {
                for (let i = istart; i <= iend; i++) {
                    //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                    tempRe += a.r[coords + k + i] * x.r[ix - x.base] - a.i[coords + k + i] * x.i[ix - x.base];
                    tempIm += a.r[coords + k + i] * x.i[ix - x.base] + a.i[coords + k + i] * x.r[ix - x.base];
                    //console.log(`i,j=(${k + i},${j}), a=(${a.r[coords + k + i]},${a.i[coords + k + i]})`);
                    ix += incx;
                }
            }
            // CONJUGATE
            else {
                for (let i = istart; i <= iend; i++) {
                    //CONJ( (a + bi) )(c+di)= (a*c+b*d)+i(a*d-b*c)
                    tempRe += a.r[coords + k + i] * x.r[ix - x.base] + a.i[coords + k + i] * x.i[ix - x.base];
                    tempIm += a.r[coords + k + i] * x.i[ix - x.base] - a.i[coords + k + i] * x.r[ix - x.base];
                    ix += incx;
                }
            }
            //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
            y.r[jy - y.base] += AlphaRe * tempRe - AlphaIm * tempIm;
            y.i[jy - y.base] += AlphaRe * tempIm + AlphaIm * tempRe;
            jy += incy;
            if (j > ku) kx += incx;
        }
    }
}
