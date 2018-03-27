import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../f_func';

/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com, Javascript Port
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*/

/*
*>
*> CSYRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*> in the second case.
*/

const { max } = Math;

export function csyrk(
    uplo: 'u' | 'l',
    trans: 't' | 'n',
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    beta: Complex,
    c: Matrix,
    ldc: number): void {

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    const tr = lowerChar(trans);
    const ul = lowerChar(uplo);
    const alphaIsZero = alpha.re === 0 && alpha.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;
    const betaIsZero = beta.re === 0 && beta.im === 0;
    const upper = ul === 'u';

    const nrowA = tr === 'n' ? n : k;

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'tn'.includes(tr)) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (k < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (lda < max(1, n)) {
        info = 10;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('csyrk', info))
    }

    //  Quick return if possible.
    if (n === 0 || ((alphaIsZero || k === 0) && betaIsOne)) return;

    //And when  alpha.eq.zero.

    if (alphaIsZero) {
        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = betaIsZero ? 0 : (betaIsOne ? c.r[coorCJ + i] : beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i]);
                c.i[coorCJ + i] = betaIsZero ? 0 : (betaIsOne ? c.r[coorCJ + i] : beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i]);
            }
        }
        return;
    }

    //*     Start the operations.
    if (tr === 'n') {
        // Form  C := alpha*A*A**T + beta*C.
        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = betaIsZero ? 0 : (betaIsOne ? c.r[coorCJ + i] : beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i]);
                c.i[coorCJ + i] = betaIsZero ? 0 : (betaIsOne ? c.i[coorCJ + i] : beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i])
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j] === 0;
                if (!aIsZero) {
                    let tempRe = alpha.re * a.r[coorAL + j] - alpha.im * a.i[coorAL + j];
                    let tempIm = alpha.re * a.i[coorAL + j] + alpha.im * a.r[coorAL + j];
                    for (let i = start; i <= stop; i++) {
                        c.r[coorCJ + i] += tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                        c.i[coorCJ + i] += tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
                    }
                }
            }
        }
    }
    else {
        //  Form  C := alpha*A**T*A + beta*C.

        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            const start = upper ? 1 : j;
            const stop = upper ? j : n;
            for (let i = start; i <= stop; i++) {
                const coorAI = a.colOfEx(i);
                let tempRe = 0;
                let tempIm = 0;
                for (let l = 1; l <= k; l++) {
                    tempRe += a.r[coorAI + l] * a.r[coorAJ + l] - a.i[coorAI + l] * a.i[coorAJ + l];
                    tempIm += a.r[coorAI + l] * a.i[coorAJ + l] + a.i[coorAI + l] * a.r[coorAJ + l];
                }
                //ALPHA*TEMP
                let re = alpha.re * tempRe - alpha.im * tempIm;
                let im = alpha.re * tempIm + alpha.im * tempRe;
                if (betaIsZero) {
                    re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }
    }
}

