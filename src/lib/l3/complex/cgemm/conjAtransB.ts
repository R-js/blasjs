import { Complex, MatrixEComplex } from '../../../f_func';
/*
Form  C := alpha*A**H*B**T + beta*C
*
              DO 310 J = 1,N
                  DO 300 I = 1,M
                      TEMP = ZERO
                      DO 290 L = 1,K
                          TEMP = TEMP + CONJG(A(L,I))*B(J,L)
  290                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  300             CONTINUE
  310         CONTINUE
*/
export function conjAtransB(
    betaIsZero: boolean,
    betaIsOne: boolean,
    beta: Complex,
    alpha: Complex,
    a: MatrixEComplex,
    b: MatrixEComplex,
    c: MatrixEComplex,
    n: number,
    m: number,
    k: number): void {


    // DO 310 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        // DO 300 I = 1,M
        for (let i = 1; i <= m; i++) {
            const coorAI = a.colOfEx(i);
            // TEMP = ZERO
            let tempRe = 0;
            let tempIm = 0;
            //DO 290 L = 1,K
            for (let l = 1; l <= k; l++) {
                const coorBL = b.colOfEx(l);
                // TEMP = TEMP + CONJG(A(L,I)) * B(J,L)
                // (a-ib)*(c+id) = (ac+bd) + i(ad-bc)
                tempRe += a.r[coorAI + l] * b.r[coorBL + j] + a.i[coorAI + l] * b.i[coorBL + j];
                tempIm += a.r[coorAI + l] * b.i[coorBL + j] - a.i[coorAI + l] * b.r[coorBL + j];
            }
            /*
            IF (BETA.EQ.ZERO) THEN
              C(I,J) = ALPHA*TEMP
            ELSE
              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
            END IF
            */
            // C(I,J) = ALPHA*TEMP
            c.r[coorCJ + i] = alpha.re * tempRe - alpha.im * tempIm;
            c.i[coorCJ + i] = alpha.re * tempIm + alpha.im * tempRe;
            if (!betaIsZero) {
                // C(I,J) = ALPHA*TEMP + beta*C(I,J)
                const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                c.r[coorCJ + i] += re;
                c.i[coorCJ + i] += im;
            }
        }
    }
}

