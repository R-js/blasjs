import { Complex, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';


//Form  C := alpha*A**H*B + beta*C.


export function conjAB(
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

    for (let j = 1; j <= n; j++) {
        const coorBJ = b.colOfEx(j);
        const coorCJ = c.colOfEx(j);
        for (let i = 1; i <= m; i++) {
            const coorAI = a.colOfEx(i);
            let tempRe = 0;
            let tempIm = 0;
            for (let l = 1; l <= k; l++) {
                // TEMP = TEMP + CONJG(A(L,I))*B(L,J)
                //(a-ib)*(c+id) = a*c+bd + i(ad-bc)
                const { re, im } = mul_rxr(
                    a.r[coorAI + l],
                    -a.i[coorAI + l],
                    b.r[coorBJ + l],
                    b.i[coorBJ + l]
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
            let { re, im } = mul_cxr(
                alpha,
                tempRe,
                tempIm
            );

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
