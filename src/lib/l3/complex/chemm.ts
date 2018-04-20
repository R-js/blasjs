

/*
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd
*/

/*
*>
*> CHEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*A*B + beta*C,
*>
*> or
*>
*>    C := alpha*B*A + beta*C,
*>
*> where alpha and beta are scalars, A is an hermitian matrix and  B and
*> C are m by n matrices.
*/

import {
    Complex,
    errMissingIm,
    errWrongArg,
    isOne,
    isZero,
    lowerChar,
    Matrix,
    mul_cxr,
    mul_rxr
} from '../../f_func';

const { max } = Math;

export function chemm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    m: number,
    n: number,
    alpha: Complex,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: Complex,
    c: Matrix,
    ldc: number): void {

    const si = lowerChar(side);
    const ul = lowerChar(uplo);

    const alphaIsZero = isZero(alpha);
    const betaIsOne = isOne(beta);
    const betaIsZero = isZero(beta);

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }
    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }

    const nrowA = (si === 'l') ? m : n;
    const upper = ul === 'u';

    // Test the input parameters.
    let info = 0;

    if (!'lr'.includes(si)) {
        info = 1;
    }
    else if (!'ul'.includes(ul)) {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (ldb < max(1, m)) {
        info = 9;
    }
    else if (ldc < max(1, m)) {
        info = 12;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chemm', info));
    }
    // Quick return if possible.

    if (m === 0 || n === 0 || (alphaIsZero && betaIsOne)) return;

    // And when  alpha.eq.zero.
    if (alphaIsZero) {
        if (betaIsZero) {
            //performance
            for (let j = 1; j <= n; j++) {
                c.setCol(j, 1, m, 0);
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
        return;
    }
    //Start the operations.
    if (si === 'l') {
        //Form  C:= alpha * A * B + beta * C.
        if (upper) {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    // console.log(`* ${i}, ${j}   (${b.r[coorBJ + i]}, ${b.i[coorBJ + i]}`);
                    const coorAI = a.colOfEx(i);
                    let temp1Re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    let temp1Im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    let temp2Re = 0;
                    let temp2Im = 0;
                    for (let k = 1; k <= i - 1; k++) {
                        c.r[k + coorCJ] += temp1Re * a.r[coorAI + k] - temp1Im * a.i[coorAI + k];
                        c.i[k + coorCJ] += temp1Re * a.i[coorAI + k] + temp1Im * a.r[coorAI + k];
                        //TEMP2 = TEMP2 + B(K,J)*DCONJG(A(K,I))
                        //(a+ib)*(c-id)=(ac+bd)+i(-ad+bc)
                        const re1 = b.r[coorBJ + k] * a.r[coorAI + k] + b.i[coorBJ + k] * a.i[coorAI + k];
                        const im1 = - b.r[coorBJ + k] * a.i[coorAI + k] + b.i[coorBJ + k] * a.r[coorAI + k];
                        //console.log(` ${k}, ${j}   (${b.r[coorBJ + k]}, ${b.i[coorBJ + k]}`);
                        temp2Re += re1;
                        temp2Im += im1;
                    }

                    /*IF(BETA.EQ.ZERO) THEN
                    C(I, J) = TEMP1 * REAL(A(I, I)) + ALPHA * TEMP2
                    ELSE
                    C(I, J) = BETA * C(I, J) + [TEMP1 * REAL(A(I, I)) + ALPHA * TEMP2]
                    END IF
                    */
                    // always do this
                    let re = temp1Re * a.r[coorAI + i] + (alpha.re * temp2Re - alpha.im * temp2Im);
                    let im = temp1Im * a.r[coorAI + i] + (alpha.re * temp2Im + alpha.im * temp2Re);
                    if (!betaIsZero) {
                        re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                        im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
        //!upper
        else {
            for (let j = 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const coorCJ = c.colOfEx(j);
                for (let i = m; i >= 1; i--) {
                    const coorAI = a.colOfEx(i);
                    let temp1Re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    let temp1Im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    let temp2Re = 0;
                    let temp2Im = 0;
                    for (let k = i + 1; k <= m; k++) {
                        //C(K, J) = C(K, J) + TEMP1 * A(K, I)
                        c.r[k + coorCJ] += temp1Re * a.r[coorAI + k] - temp1Im * a.i[coorAI + k];
                        c.i[k + coorCJ] += temp1Re * a.i[coorAI + k] + temp1Im * a.r[coorAI + k];

                        //TEMP2 = TEMP2 + B(K, J) * CONJG(A(K, I))
                        temp2Re += b.r[coorBJ + k] * a.r[coorAI + k] + b.i[coorBJ + k] * a.i[coorAI + k];
                        temp2Im += -b.r[coorBJ + k] * a.i[coorAI + k] + b.i[coorBJ + k] * a.r[coorAI + k];
                    }


                    /*
                     IF (BETA.EQ.ZERO) THEN
                          C(I,J) = TEMP1*REAL(A(I,I)) + ALPHA*TEMP2
                      ELSE
                          C(I,J) = BETA*C(I,J) + [TEMP1*REAL(A(I,I)) + ALPHA*TEMP2]
                      END IF
                    */


                    let re = temp1Re * a.r[coorAI + i] - 0 + (alpha.re * temp2Re - alpha.im * temp2Im);
                    let im = 0 + temp1Im * a.r[coorAI + i] + (alpha.re * temp2Im + alpha.im * temp2Re);
                    //console.log({ re, im });
                    if (!betaIsZero) {
                        re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                        im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                    }
                    c.r[coorCJ + i] = re;
                    c.i[coorCJ + i] = im;
                }
            }
        }
    }
    else {
        // Form  C := alpha*B*A + beta*C.

        for (let j = 1; j <= n; j++) {
            const coorAJ = a.colOfEx(j);
            const coorBJ = b.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            let temp1Re = alpha.re * a.r[coorAJ + j];
            let temp1Im = alpha.im * a.r[coorAJ + j];
            /*  IF (BETA.EQ.ZERO) THEN
                    DO 110 I = 1,M
                        C(I,J) = [TEMP1*B(I,J)]
110                     CONTINUE
                ELSE
                    DO 120 I = 1,M
                        C(I,J) = BETA*C(I,J) + [TEMP1*B(I,J)]
120                     CONTINUE
                END IF
            */
            for (let i = 1; i <= m; i++) {
                let re = temp1Re * b.r[coorBJ + i] - temp1Im * b.i[coorBJ + i];
                let im = temp1Re * b.i[coorBJ + i] + temp1Im * b.r[coorBJ + i];
                if (!betaIsZero) {
                    re += beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                    im += beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                }
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
            //DO 140 K = 1,J - 1
            for (let k = 1; k <= j - 1; k++) {
                const coorBK = b.colOfEx(k);
                const coorAK = a.colOfEx(k);
                /*
                 IF (UPPER) THEN
                    TEMP1 = ALPHA*A(K,J)
                ELSE
                    //(a+ib)*(c-id)=(ac+bd)+i(-ad+bc)
                    TEMP1 = ALPHA*CONJG(A(J,K))
                END IF
                */
                //
                const { re, im } = mul_cxr(
                    alpha,
                    upper ? a.r[coorAJ + k] : a.r[coorAK + j],
                    upper ? a.i[coorAJ + k] : -a.i[coorAK + j]
                );
                let temp1Re = re;
                let temp1Im = im;
                /*
                DO 130 I = 1,M
                    C(I,J) = C(I,J) + TEMP1*B(I,K)
130             CONTINUE
                */
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1Re * b.r[coorBK + i] - temp1Im * b.i[coorBK + i];
                    c.i[coorCJ + i] += temp1Re * b.i[coorBK + i] + temp1Im * b.r[coorBK + i];
                }
            }
            //DO 160 K = J + 1,N
            for (let k = j + 1; k <= n; k++) {
                const coorBK = b.colOfEx(k);
                const coorAK = a.colOfEx(k);
                /*
                 IF (UPPER) THEN
                    //(a+ib)*(c-id)=(ac+bd)+i(-ad+bc)
                    TEMP1 = ALPHA*CONJG(A(J,K))
                ELSE
                    TEMP1 = ALPHA*A(K,J)
                END IF
                */
                const { re, im } = mul_cxr(
                    alpha,
                    upper ? a.r[coorAK + j] : a.r[coorAJ + k],
                    upper ? -a.i[coorAK + j] : a.i[coorAJ + k]
                );
                let temp1Re = re;
                let temp1Im = im;

                /*
                DO 150 I = 1,M
                    C(I,J) = C(I,J) + TEMP1*B(I,K)
150             CONTINUE
                */
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1Re * b.r[coorBK + i] - temp1Im * b.i[coorBK + i];
                    c.i[coorCJ + i] += temp1Re * b.i[coorBK + i] + temp1Im * b.r[coorBK + i];
                }
            }
        }
    }
}
