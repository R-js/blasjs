import { Complex, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';

// //Form  C := alpha*A**T*B**H + beta*C


/*
 DO 340 J = 1,N
                  DO 330 I = 1,M
                      TEMP = ZERO
                      DO 320 L = 1,K
                          TEMP = TEMP + A(L,I)*CONJG(B(J,L))
  320                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  330             CONTINUE
  340         CONTINUE
*/
export function transAconjB(
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

    //DO 340 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        // DO 330 I = 1,M
        for (let i = 1; i <= m; i++) {
            const coorAI = a.colOfEx(i);
            //TEMP = ZERO
            let tempRe = 0;
            let tempIm = 0;
            // DO 320 L = 1,K
            for (let l = 1; l <= k; l++) {
                const coorBL = b.colOfEx(l);
                //  TEMP = TEMP + A(L,I)*CONJG(B(J,L))
                //(a+ib)*(c-id) = ac+bd + i(-ad+bc)
                let { re, im } = mul_rxr(
                    a.r[coorAI + l],
                    a.i[coorAI + l],
                    b.r[coorBL + j],
                    -b.i[coorBL + j],
                );
                tempRe += re;
                tempIm += im;
            }
            /*
            IF (BETA.EQ.ZERO) THEN
              C(I,J) = ALPHA*TEMP
            ELSE
              C(I,J) = ALPHA*TEMP + BETA*C(I,J)
            END IF
            */
            // C(I,J) = ALPHA*TEMP
            let { re, im } = mul_cxr(alpha, tempRe, tempIm);
            if (!betaIsZero) {
                // C(I,J) = ALPHA*TEMP + beta*C(I,J)
                const { re: re1, im: im1 } = mul_cxr(
                    beta,
                    c.r[coorCJ + i],
                    c.i[coorCJ + i]
                );
                re += re1;
                im += im1;
            }
            c.r[coorCJ + i] = re;
            c.i[coorCJ + i] = im;
        }
    }
}
