import { Complex, MatrixEComplex } from '../../../f_func';

//Form  C := alpha*A*B**T + beta*C
export function AtransB(
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

    //
    // Form  C := alpha*A*B**T + beta*C
    //
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        //  IF (BETA.EQ.ZERO) THEN
        if (betaIsZero) {
            c.setCol(j, 1, m, 0);
        }
        else if (!betaIsOne) {
            // DO 210 I = 1,M
            for (let i = 1; i <= m; i++) {
                // C(I,J) = BETA*C(I,J)
                const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }

        for (let l = 1; l <= k; l++) {
            const coorBL = b.colOfEx(l);
            const coorAL = a.colOfEx(l);

            // TEMP = ALPHA*B(J,L)
            //(a+ib)*(c+id) = (ac-bd)+i(ad+bc)
            let tempRe = alpha.re * b.r[coorBL + j] - alpha.im * b.i[coorBL + j];
            let tempIm = alpha.re * b.i[coorBL + j] + alpha.im * b.r[coorBL + j];
            for (let i = 1; i <= m; i++) {
                // // //   C(I,J) = C(I,J) + TEMP*A(I,L)
                c.r[coorCJ + i] += tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                c.i[coorCJ + i] += tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
            }
        }
    }
    //Form  C := alpha*B*A + beta*C.
}
