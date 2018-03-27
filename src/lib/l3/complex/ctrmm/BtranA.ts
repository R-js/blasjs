
import { Complex, errMissingIm, Matrix } from '../../../f_func';


export function BtranA(
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

    if (upper) {
        for (let k = 1; k <= n; k++) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            for (let j = 1; j <= k - 1; j++) {
                const coorBJ = b.colOfEx(j);
                const aIsZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!aIsZero) {
                    let ajkRe = a.r[coorAK + j];
                    let ajkIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];
                    //TEMP = ALPHA*A(J,K)
                    //TEMP = ALPHA*CONJG(A(J,K))
                    let tempRe = alpha.re * ajkRe - alpha.im * ajkIm;
                    let tempIm = alpha.re * ajkIm + alpha.im * ajkRe;
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) + TEMP*B(I,K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.i[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }//j
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            if (nounit) {
                let akkRe = a.r[coorAK + k];
                let akkIm = noconj ? a.i[coorAK + k] : -a.i[coorAK + k];

                let re = tempRe * akkRe - tempIm * akkIm;
                let im = tempRe * akkIm + tempIm * akkRe;

                tempRe = re;
                tempIm = im;
            }
            if (!(tempRe === 1 && tempIm === 0)) {
                for (let i = 1; i <= m; i++) {
                    let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }
        }//k
    }//upper
    else {
        for (let k = n; k >= 1; k--) {
            const coorAK = a.colOfEx(k);
            const coorBK = b.colOfEx(k);
            for (let j = k + 1; j <= n; j++) {
                const coorBJ = b.colOfEx(j);
                const aIsZero = a.r[coorAK + j] === 0 && a.i[coorAK + j] === 0;
                if (!aIsZero) {
                    const ajkRe = a.r[coorAK + j];
                    const ajkIm = noconj ? a.i[coorAK + j] : -a.i[coorAK + j];

                    const tempRe = alpha.re * ajkRe - alpha.im * ajkIm;
                    const tempIm = alpha.re * ajkIm + alpha.im * ajkRe;
                    for (let i = 1; i <= m; i++) {
                        //B(I,J) = B(I,J) + TEMP*B(I,K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }//for(j)
            let tempRe = alpha.re;
            let tempIm = alpha.im;

            if (nounit) {
                const akkRe = a.r[coorAK + k];
                const akkIm = noconj ? a.i[coorAK + k] : - a.i[coorAK + k];

                const re = tempRe * akkRe - tempIm * akkIm;
                const im = tempRe * akkIm + tempIm * akkRe;

                tempRe = re;
                tempIm = im;
            }
            if (!(tempRe === 1 && tempIm === 0)) {
                for (let i = 1; i <= m; i++) {
                    //  B(I,K) = TEMP*B(I,K)
                    const re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                    const im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                    b.r[coorBK + i] = re;
                    b.i[coorBK + i] = im;
                }
            }

        }//k
    }//upper
}
