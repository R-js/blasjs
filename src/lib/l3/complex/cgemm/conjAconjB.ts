import { Complex, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';

//  Form  C := alpha*A**H*B**H + beta*C.
export function conjAconjB(
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

    //DO 280 J = 1,N
    for (let j = 1; j <= n; j++) {
        //const coorBJ = b.colOfEx(j);
        const coorCJ = c.colOfEx(j);
        //DO 270 I = 1,M
        for (let i = 1; i <= m; i++) {
            //TEMP = ZERO
            const coorAI = a.colOfEx(i);
            let tempRe = 0;
            let tempIm = 0;
            //DO 260 L = 1,K
            for (let l = 1; l <= k; l++) {
                const coorBL = b.colOfEx(l);
                // TEMP = TEMP + CONJG(A(L,I))*CONJB(L,J)
                //(a-ib)*(c-id) = (a*c-bd) + i(-ad-bc)
                const { re, im } = mul_rxr(
                    a.r[coorAI + l],
                    -a.i[coorAI + l],
                    b.r[coorBL + j],
                    -b.i[coorBL + j]
                );
                tempRe += re;
                tempIm += im;
            }
            let { re, im } = mul_cxr(
                alpha,
                tempRe,
                tempIm
            );
            if (!betaIsZero) {
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
