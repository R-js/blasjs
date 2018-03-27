import { Complex, errMissingIm, Matrix } from '../../../f_func';

/*
Form  B := alpha*B*inv( A**T )
or    B := alpha*B*inv( A**H ).
*/

export function BinvTranConjA(
    nounit: boolean,
    upper: boolean,
    noconj: boolean,
    n: number,
    m: number,
    a: Matrix,
    b: Matrix,
    alpha: Complex): void {


    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }

    const alphaIsOne = alpha.re === 1 && alpha.im === 0;

    if (upper) {
        for (let k = n; k >= 1; k--) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            if (nounit) {
                // (1+i0)/(c+id), a=1,b=0
                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                const akkRe = a.r[coorAK + k];
                const akkIm = noconj ? a.i[coorAK + k] : -a.i[coorAK + k];

                const n = akkRe * akkRe + akkIm * akkIm;
                let tempRe = akkRe / n;
                let tempIm = -akkIm / n;
                for (let i = 1; i <= m; i++) {
                    const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
            for (let j = 1; j <= k - 1; j++) {
                const aIsZero = a.r[coorBK + j] === 0 && a.i[coorBK + j] === 0;
                if (!aIsZero) {
                    let tempRe = a.r[coorAK + j];
                    let tempIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];
                    for (let i = 1; i <= m; i++) {
                        //  B(I,J) = B(I,J) - TEMP*B(I,K)
                        const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBK + i] -= re;
                        b.r[coorBK + i] -= im;
                    }
                }
            }
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    let re = alpha.re * b.r[coorBK + i] - alpha.im * b.i[coorBK + i];
                    let im = alpha.re * b.i[coorBK + i] + alpha.im * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        }//k
    }
    else {

    }


}
