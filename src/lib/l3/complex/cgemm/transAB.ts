import { Complex, MatrixEComplex } from '../../../f_func';

//Form  C := alpha*A**T*B + beta*C

export function transAB(
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

    // Form  C := alpha*A**T*B + beta*C
    //DO 150 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorBJ = b.colOfEx(j);
        const coorCJ = c.colOfEx(j);

        for (let i = 1; i <= m; i++) {
            const coorAI = a.colOfEx(i);
            let tempRe = 0;
            let tempIm = 0;
            for (let l = 1; i <= k; l++) {
                // TEMP = TEMP + A(L, I) * B(L, J)
                tempRe += a.r[coorAI + l] * b.r[coorBJ + l] - a.i[coorAI + l] * b.i[coorBJ + l];
                tempIm += a.r[coorAI + l] * b.i[coorBJ + l] + a.i[coorAI + l] * b.r[coorBJ + l];
            }
            /*
             IF(BETA.EQ.ZERO) THEN{
                    C(I, J) = ALPHA * TEMP
                } ELSE{
                    C(I, J) = ALPHA * TEMP + BETA * C(I, J)
            } END IF
            */
            c.r[coorCJ + i] = alpha.re * tempRe - alpha.im * tempIm;
            c.i[coorCJ + i] = alpha.re * tempIm + alpha.im * tempRe;
            if (!betaIsZero) {
                const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                c.r[coorCJ + i] += re;
                c.i[coorCJ + i] += im;
            }
        }
    }
}
