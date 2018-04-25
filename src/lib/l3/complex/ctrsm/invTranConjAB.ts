import {
    Complex,
    div_rxr,
    MatrixEComplex,
    mul_cxr
} from '../../../f_func';

//Form  B := alpha*inv( A**T )*B
//or    B := alpha*inv( A**H )*B.


export function invTranConjAB(
    nounit: boolean,
    upper: boolean,
    alphaIsOne: boolean,
    alphaIsZero: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex
): void {


    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                const coorAI = a.colOfEx(i);
                const { re, im } = mul_cxr(
                    alpha,
                    b.r[coorBJ + i],
                    b.i[coorBJ + i]
                )
                let tempRe = re; //alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                let tempIm = im; //alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                for (let k = 1; k <= i - 1; k++) {
                    //A(K,I)*B(K,J)
                    //CONJG(A(K,I))*B(K,J) if noconj=false
                    const akiRe = a.r[coorAI + k];
                    const akiIm = noconj ? a.i[coorAI + k] : -a.i[coorAI + k];

                    const re = akiRe * b.r[coorBJ + k] - akiIm * b.i[coorBJ + k];
                    const im = akiRe * b.i[coorBJ + k] + akiIm * b.r[coorBJ + k];

                    tempRe -= re;
                    tempIm -= im;
                }
                if (nounit) {
                    //(a+ib)/(c+id)
                    // re= (ac+bd)/(c*c+d*d)
                    // im =(bc-ad)/(c*c+d*d)
                    // TEMP/A(I,I) | TEMP/CONJG(A(I,I))
                    const { re, im } = div_rxr(
                        tempRe,
                        tempIm,
                        a.r[coorAI + i],
                        a.i[coorAI + i]
                    )
                    /*const _c = a.r[coorAI + i];
                    const _d = noconj ? a.i[coorAI + i] : -a.i[coorAI + i];

                    const _a = tempRe;
                    const _b = tempIm;
                    const n = _c * _c + _d * _d;
                    */
                    tempRe = re; //(_a * _c + _b * _d) / n;
                    tempIm = im; //(_b * _c - _a * _d) / n;
                }
                b.r[coorBJ + i] = tempRe;
                b.i[coorBJ + i] = tempIm;
            }//i
        }//j
    }//upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let i = m; i >= 1; i--) {
                const coorAI = a.colOfEx(i);
                let tempRe = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                let tempIm = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                for (let k = i + 1; k <= m; k++) {
                    //A(K,I)*B(K,J)
                    //CONJG(A(K,I))*B(K,J) if noconj=false
                    const akiRe = a.r[coorAI + k];
                    const akiIm = noconj ? a.i[coorAI + k] : -a.i[coorAI + k];

                    const re = akiRe * b.r[coorBJ + k] - akiIm * b.i[coorBJ + k];
                    const im = akiRe * b.i[coorBJ + k] + akiIm * b.r[coorBJ + k];

                    tempRe -= re;
                    tempIm -= im;
                }
                if (nounit) {
                    //(a+ib)/(c+id)
                    // re= (ac+bd)/(c*c+d*d)
                    // im =(bc-ad)/(c*c+d*d)
                    // TEMP/A(I,I) | TEMP/CONJG(A(I,I))

                    const _c = a.r[coorAI + i];
                    const _d = noconj ? a.i[coorAI + i] : -a.i[coorAI + i];

                    const _a = tempRe;
                    const _b = tempIm;
                    const n = _c * _c + _d * _d;

                    tempRe = (_a * _c + _b * _d) / n;
                    tempIm = (_b * _c - _a * _d) / n;
                }
                b.r[coorBJ + i] = tempRe;
                b.i[coorBJ + i] = tempIm;
            }//i
        }//j
    }
}
