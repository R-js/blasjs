import { Complex, errMissingIm, Matrix } from '../../../f_func';


export function AB(beta: Complex, alpha: Complex, a: Matrix, b: Matrix, c: Matrix, n: number, m: number, k: number): void {

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

    // DO 90 J = 1,N
    for (let j = 1; j <= n; j++) {
        const coorCJ = c.colOfEx(j);
        const coorBJ = b.colOfEx(j);
        // IF (BETA.EQ.ZERO) THEN
        if (betaIsZero) {
            c.setCol(j, 1, m, 0);
        }
        //IF (BETA.NE.ONE) THEN
        else if (!betaIsOne) {
            /// DO 60 I = 1,M
            for (let i = 1; i <= m; i++) {
                // C(I,J) = BETA*C(I,J)
                c.r[coorCJ + i] = beta.re * c.r[coorCJ + i] - beta.im * c.r[coorCJ + i];
                c.i[coorCJ + i] = beta.re * c.i[coorCJ + i] + beta.im * c.i[coorCJ + i];
            }
        }
        //DO 80 L = 1,K
        for (let l = 1; l <= k; l++) {
            //TEMP = ALPHA*B(L,J)
            const coorAL = a.colOfEx(l);
            let tempRe = alpha.re * b.r[coorBJ + l] - alpha.im * b.i[coorBJ + l];
            let tempIm = alpha.re * b.r[coorBJ + l] + alpha.im * b.i[coorBJ + l];
            //DO 70 I = 1,M
            for (let i = 1; i <= m; i++) {
                // C(I,J) = C(I,J) + TEMP*A(I,L)
                c.r[coorCJ + i] += tempRe * a.r[coorAL + i] - tempIm * a.i[coorAL + i];
                c.i[coorCJ + i] += tempRe * a.i[coorAL + i] + tempIm * a.r[coorAL + i];
            }
        }
    }
}
