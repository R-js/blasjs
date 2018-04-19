import {
    Complex,
    MatrixEComplex,
    mul_cxr,
    mul_rxr
} from '../../../f_func';


export function AB(
    nounit: boolean,
    upper: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex): void {

    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let k = 1; k <= m; k++) {
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!bIsZero) {
                    // TEMP = ALPHA*B(K,J)
                    const coorAK = a.colOfEx(k);
                    const { re, im } = mul_cxr(alpha, b.r[coorBJ + k], b.i[coorBJ + k]);
                    let tempRe = re;
                    let tempIm = im;

                    for (let i = 1; i <= k - 1; i++) {
                        //TEMP*A(I,K)
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAK + i], a.i[coorAK + i]);
                        //B(I,J) = B(I,J) + ...
                        // console.log(`${j},${k},${i},\t (${re},${im})`);
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                        //                console.log(`${j},${k},${i},\t (${b.r[coorBJ + i]},${b.i[coorBJ + i]})`);

                    }
                    if (nounit) {
                        //TEMP = TEMP*A(K,K)
                        const { re, im } = mul_rxr(tempRe, tempIm, a.r[coorAK + k], a.i[coorAK + k]);
                        tempRe = re;
                        tempIm = im;
                    }
                    //B(K,J) = TEMP
                    b.r[coorBJ + k] = tempRe;
                    b.i[coorBJ + k] = tempIm;
                }//if b[*]!=zero
            }//for k
        }//for j
    }
    else {
        for (let j = 1; j <= n; j++) {
            for (let k = m; k >= 1; k--) {
                const coorAK = a.colOfEx(j);
                const coorBJ = b.colOfEx(j);
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!bIsZero) {
                    let tempRe = alpha.re * b.r[coorBJ + k] - alpha.im * b.i[coorBJ + k];
                    let tempIm = alpha.re * b.r[coorBJ + k] + alpha.im * b.i[coorBJ + k];
                    b.r[coorBJ + k] = tempRe;
                    b.i[coorBJ + k] = tempIm;
                    //IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                    if (nounit) {
                        let re = b.r[coorBJ + k] * a.r[coorAK + k] - b.i[coorBJ + k] * a.i[coorAK + k];
                        let im = b.r[coorBJ + k] * a.i[coorAK + k] + b.i[coorBJ + k] * a.i[coorAK + k];
                        b.r[coorBJ + k] = re;
                        b.i[coorBJ + k] = im;
                    }
                    for (let i = k + 1; i <= m; i++) {
                        //TEMP*A(I,K)
                        let re = tempRe * a.r[coorAK + i] - tempIm * a.i[coorAK + i];
                        let im = tempRe * a.i[coorAK + i] + tempIm * a.r[coorAK + i];
                        //B(I,J) = B(I,J) + TEMP*A(I,K)
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }
        }
    }
}
