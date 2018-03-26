import { Complex, errMissingIm, errWrongArg, lowerChar, Matrix } from '../../f_func';

/*
*>  -- Jacob Bogers, Javascript Port, 03/2018, jkfbogers@gmail.com
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*>
*>  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
*>     Ed Anderson, Cray Research Inc.
*>
*>
*> CHER2K  performs one of the hermitian rank 2k operations
*>
*>    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
*>
*> or
*>
*>    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
*>
*> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
*> hermitian matrix and  A and B  are  n by k matrices in the first case
*> and  k by n  matrices in the second case.
*/

const { max } = Math;

export function cher2k(
    uplo: 'u' | 'l',
    trans: 'n' | 'c',
    n: number,
    k: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: number,
    c: Matrix,
    ldc: number): void {

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    const tr = lowerChar(trans);
    const ul = lowerChar(uplo);

    const alphaIsZero = alpha.re === 0 && alpha.im === 0;

    const nrowA = trans === 'n' ? n : k;
    const upper = ul === 'u';

    let info = 0;

    if ('ul'.includes(ul)) {
        info = 1;
    }
    else if ('nc'.includes(tr)) {
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
    else if (ldb < max(1, nrowA)) {
        info = 9;
    }
    else if (ldc < max(1, n)) {
        info = 12;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('cher2k', info));
    }

    //     Quick return if possible.

    /*IF(
            (N.EQ.0).OR.(
                (
                    (ALPHA.EQ.ZERO).OR.(K.EQ.0)
                ).AND.(BETA.EQ.ONE)
            )
    ) RETURN*/

    if (n === 0 || ((alphaIsZero || k === 0) && beta === 1)) return;

    //* And when  alpha.eq.zero.
    if (alphaIsZero) {
        if (upper) {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, 1, j, 0);
                }
                else {
                    for (let i = 1; i <= j - 1; i++) {
                        c.r[coorCJ + i] *= beta;
                        c.i[coorCJ + i] *= beta;
                    }
                    c.r[coorCJ + j] = beta;
                    c.i[coorCJ + j] = 0;
                }
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                if (beta === 0) {
                    c.setCol(j, j, n, 0);
                }
                else {
                    const coorCJ = c.colOfEx(j);
                    c.r[coorCJ + j] *= beta;
                    c.i[coorCJ + j] = 0;
                    for (let i = j + 1; j <= n; i++) {
                        c.r[coorCJ + i] *= beta;
                        c.r[coorCJ + i] *= beta;
                    }
                }
            }
        }
        return;
    }

    //Start the operations.

    if (tr === 'n') {
        // Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
        if (upper) {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, 1, j, 0);
                }
                else if (beta !== 1) {
                    for (let i = 1; i <= j - 1; i++) {
                        c.r[coorCJ + i] *= beta;
                        c.i[coorCJ + i] *= beta;
                    }
                    c.r[coorCJ + j] *= beta;
                    c.i[coorCJ + j] = 0;
                }
                else {
                    c.i[coorCJ + j] = 0;
                }
                for (let l = 1; l <= k; l++) {
                    const coorAL = a.colOfEx(l);
                    const coorBL = b.colOfEx(l);
                    const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j];
                    const bisZero = b.r[coorBL + j] === 0 && b.i[coorBL + j];
                    if (!aIsZero || !bisZero) {
                        //TEMP1 = ALPHA * CONJG(B(J, L))
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        //
                        let temp1Re = alpha.re * b.r[coorBL + j] + alpha.im * b.i[coorBL + j];
                        let temp1Im = -alpha.re * b.i[coorBL + j] + alpha.im * b.r[coorBL + j];

                        let temp2Re = alpha.re * a.r[coorAL + j] - alpha.im * a.i[coorAL + j];
                        let temp2Im = alpha.re * a.i[coorAL + j] + alpha.im * a.r[coorAL + j];

                        temp2Im = -temp2Im;
                        for (let i = 1; i <= j; i++) {
                            //C(I, J) = C(I, J) + [A(I, L) * TEMP1] +
                            //+                                 B(I, L) * TEMP2
                            c.r[coorCJ + i] += a.r[coorAL + i] * temp1Re - a.i[coorAL + i] * temp1Im;
                            c.i[coorCJ + i] += a.r[coorAL + i] * temp1Im + a.i[coorAL + i] * temp1Re;

                            c.r[coorCJ + i] += b.r[coorBL + i] * temp2Re - b.i[coorBL + i] * temp2Im;
                            c.i[coorCJ + i] += b.r[coorBL + i] * temp2Im + b.i[coorBL + i] * temp2Re;
                        }
                        c.i[coorCJ + j] = 0;
                    }
                }
            }
        }
        //not upper
        else {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                if (beta === 0) {
                    c.setCol(j, j, n, 0);
                }
                else if (beta !== 1) {
                    for (let i = j + 1; i <= n; i++) {
                        c.r[coorCJ + i] *= beta;
                    }
                    c.r[coorCJ + j] *= beta;
                    c.i[coorCJ + j] = 0;
                }
                else {//beta =1
                    c.i[coorCJ + j] = 0;
                }
                for (let l = 1; l <= k; l++) {
                    const coorAL = a.colOfEx(l);
                    const coorBL = b.colOfEx(l);
                    const aIsZero = a.r[coorAL + j] === 0 && a.i[coorAL + j] === 0;
                    const bIsZero = b.r[coorBL + j] === 0 && b.i[coorBL + j] === 0;
                    if (!aIsZero || !bIsZero) {
                        //TEMP1 = ALPHA * CONJG(B(J, L))
                        //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                        let temp1Re = alpha.re * b.r[coorBL + j] + alpha.im * b.i[coorBL + j];
                        let temp1Im = -alpha.re * b.i[coorBL + j] + alpha.im * b.i[coorBL + j];

                        let temp2Re = alpha.re * a.r[coorAL + j] - alpha.im * a.i[coorAL + j];
                        let temp2Im = -(alpha.re * a.i[coorAL + j] + alpha.im * a.i[coorAL + j]);
                        for (let i = j; i <= n; i++) {
                            //DO 160 I = J + 1, N, note we use starting i=j (because of cleanup after)
                            // A(I, L) * TEMP1
                            c.r[coorCJ + i] += a.r[coorAL + i] * temp1Re - a.i[coorAL + i] * temp1Im;
                            c.i[coorCJ + i] += a.r[coorAL + i] * temp1Im + a.i[coorAL + i] * temp1Re;
                            // B(I, L) * TEMP2
                            c.r[coorCJ + i] += b.r[coorBL + i] * temp2Re - b.i[coorBL + i] * temp2Re;
                            c.i[coorCJ + i] += b.r[coorBL + i] * temp2Im + b.i[coorBL + i] * temp2Im;
                        }
                        c.i[coorCJ + j] = 0;
                    }
                }//for (k)
            }//for (j)
        }// end upper
    }//tr === 'c', conjugate
    else {
        // Form  C := alpha*A**H*B + conjg( alpha )*B**H*A + C.
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);

            const start = upper ? 1 : j;
            const stop = upper ? j : n;

            for (let i = start; i <= stop; i++) {
                const coorAI = a.colOfEx(i);
                const coorBI = b.colOfEx(i);
                let temp1Re = 0;
                let temp1Im = 0;
                let temp2Re = 0;
                let temp2Im = 0;

                for (let l = 1; l <= k; l++) {
                    //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                    //CONJG(A(L, I)) * B(L, J)
                    temp1Re += a.r[coorAI + l] * b.r[coorBJ + l] + a.i[coorAI + l] * b.i[coorBJ + l];
                    temp1Im += a.r[coorAI + l] * b.i[coorBJ + l] - a.i[coorAI + l] * b.r[coorBJ + l];
                    //CONJG(B(L, I)) * A(L, J)
                    //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                    temp2Re += b.r[coorBI + l] * a.r[coorAJ + l] + b.i[coorBI + l] * a.i[coorAJ + l];
                    temp2Im += b.r[coorBI + l] * a.i[coorAJ + l] - b.i[coorBI + l] * a.r[coorAJ + l];
                }
                // done more efficiently then in fortran
                // (ALPHA*TEMP1+CONJG(ALPHA)*TEMP2)

                let re = (alpha.re * temp1Re - alpha.im * temp1Im) + (alpha.re * temp2Re + alpha.im * temp2Im);
                let im = (i === j) ? 0 : (alpha.re * temp1Im + alpha.im * temp1Re) + (alpha.re * temp2Im - alpha.im * temp2Re);

                if (beta !== 0) {
                    re += beta * c.r[coorCJ + i];
                    im += (i === j) ? 0 : beta * c.i[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }
    }//tr === 'c'
}
