import { Complex, div_rxr, MatrixEComplex, mul_cxr, mul_rxr } from '../../../f_func';

/*
Form  B := alpha*B*inv( A**T )
or    B := alpha*B*inv( A**H ).
*/

export function BinvTranConjA(
    nounit: boolean,
    upper: boolean,
    alphaIsOne: boolean,
    alphaIsZero: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: MatrixEComplex,
    b: MatrixEComplex,
    alpha: Complex): void {

    if (upper) {
        for (let k = n; k >= 1; k--) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            if (nounit) {
                const { re: tempRe, im: tempIm } = div_rxr(
                    1,
                    0,
                    a.r[coorAK + k],
                    noconj ? a.i[coorAK + k] : -a.i[coorAK + k]
                );
                for (let i = 1; i <= m; i++) {
                    const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];

                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
            for (let j = 1; j <= k - 1; j++) {
                const coorBJ = b.colOfEx(j);
                const aIsZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!aIsZero) {
                    const tempRe = a.r[coorAK + j];
                    const tempIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];
                    ///   console.log(`k:${k},j:${j}\t(${tempRe},${tempIm})`);

                    for (let i = 1; i <= m; i++) {
                        //  B(I,J) = B(I,J) - TEMP*B(I,K)
                        const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                        // console.log(`k:${i},j:${j}\t(${coorBJ + i},${im})`);

                    }//for
                }//if
            }//for
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBK + i] - alpha.im * b.i[coorBK + i];
                    const im = alpha.re * b.i[coorBK + i] + alpha.im * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        }//k
    }//upper
    else {
        for (let k = 1; k <= n; k++) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);

            if (nounit) {
                const akkRe = a.r[coorAK + k];
                const akkIm = noconj ? a.i[coorAK + k] : -a.i[coorAK + k];
                const n = akkRe * akkRe + akkIm * akkIm;
                // (1+i0)/(c+id), a=1,b=0

                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                let tempRe = akkRe / n;
                let tempIm = -akkIm / n;
                for (let i = 1; i <= m; i++) {
                    let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
            for (let j = k + 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const isAZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!isAZero) {
                    let tempRe = a.r[coorAK + j];
                    let tempIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];
                    for (let i = 1; i <= m; i++) {
                        const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                }
            }//j
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBK + i] - alpha.im * b.i[coorBK + i];
                    const im = alpha.re * b.i[coorBK + i] + alpha.im * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        }//k
    }
}
