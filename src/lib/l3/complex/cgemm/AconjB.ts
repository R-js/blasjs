
import { Complex, errMissingIm, Matrix } from '../../../f_func';


// Form  C := alpha*A*B**H + beta*C.
/*
DO 200 J = 1,N
    IF (BETA.EQ.ZERO) THEN
        DO 160 I = 1,M
            C(I,J) = ZERO
160                 CONTINUE
    ELSE IF (BETA.NE.ONE) THEN
        DO 170 I = 1,M
            C(I,J) = BETA*C(I,J)
170                 CONTINUE
    END IF
    DO 190 L = 1,K
        TEMP = ALPHA*CONJG(B(J,L))
        DO 180 I = 1,M
            C(I,J) = C(I,J) + TEMP*A(I,L)
180                 CONTINUE
190             CONTINUE
200         CONTINUE
*/

export function AconjB(beta: Complex, alpha: Complex, a: Matrix, b: Matrix, c: Matrix, n: number, m: number, k: number): void {

    const betaIsZero = beta.re === 0 && beta.im === 0;
    const betaIsOne = beta.re === 1 && beta.im === 0;


    if (c.i === undefined) {
        throw new Error(errMissingIm('c.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    //  DO 120 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);

        //IF (BETA.EQ.ZERO) THEN
        if (betaIsZero) {
            c.setCol(j, 1, m, 0);
        }
        // ELSE IF (BETA.NE.ONE) THEN
        else if (!betaIsOne) {
            // DO 170 I = 1,M
            for (let i = 1; i <= m; i++) {
                // C(I,J) = BETA*C(I,J)
                const re = beta.re * c.r[coorCJ + i] - beta.im * c.i[coorCJ + i];
                const im = beta.re * c.i[coorCJ + i] + beta.im * c.r[coorCJ + i];
                c.r[coorCJ + i] = re;
                c.i[coorCJ + i] = im;
            }
        }

        // DO 190 L = 1,K
        for (let l = 1; l <= k; l++) {
            const coorBL = b.colOfEx(l);
            const coorAL = a.colOfEx(l);
            // TEMP = ALPHA*CONJG(B(J,L))
            //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
            let tempRe = alpha.re * b.r[coorBL + j] + alpha.re * b.i[coorBL + j];
            let tempIm = -alpha.re * b.i[coorBL + j] + alpha.re * b.r[coorBL + j];

            //DO 180 I = 1,M
            for (let i = 1; i <= m; i++) {
                //C(I,J) = C(I,J) + TEMP*A(I,L)
                c.r[coorCJ + i] += tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                c.i[coorCJ + i] += tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
            }
        }
    }
}
