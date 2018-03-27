
import { Complex, errMissingIm, Matrix } from '../../../f_func';

// Form  B := alpha*B*A.

export function BA(nounit: boolean, upper: boolean, n: number, m: number, a: Matrix, b: Matrix, alpha: Complex): void {

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }

    if (upper) {
        for (let j = n; j >= 1; j--) {
            const coorAJ = a.colOfEx(j);
            const coorBJ = b.colOfEx(j);
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            //IF (NOUNIT) TEMP = TEMP*A(J,J)
            if (nounit) {
                let re = tempRe * a.r[coorAJ + j] - tempIm * a.i[coorAJ + j];
                let im = tempRe * a.i[coorAJ + j] + tempIm * a.r[coorAJ + j];
                tempRe = re;
                tempIm = im;
            }
            for (let i = 1; i <= m; i++) {
                // B(I,J) = TEMP*B(I,J)
                // (a-ib)*(c+id) = (ab+bd)+i(ad-bc)
                let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                b.r[coorBJ + i] = re;
                b.i[coorBJ + i] = im;
            }
            //DO 190 K = 1,J - 1
            for (let k = 1; k <= j - 1; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    //TEMP = ALPHA * A(K, J)
                    let tempRe = alpha.re * a.r[coorAJ + k] - alpha.im * a.i[coorAJ + k];
                    let tempIm = alpha.re * a.i[coorAJ + k] + alpha.im * a.r[coorAJ + k];
                    for (let i = 1; i <= m; i++) {
                        // B(I, J) = B(I, J) + TEMP * B(I, K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }//k
        }//j
    }//upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            let tempRe = alpha.re;
            let tempIm = alpha.im;
            if (nounit) {
                //  TEMP = TEMP * A(J, J)
                let re = tempRe * a.r[coorAJ + j] - tempIm * a.i[coorAJ + j];
                let im = tempRe * a.i[coorAJ + j] + tempIm * a.r[coorAJ + j];
                tempRe = re;
                tempIm = im;
            }
            for (let i = 1; i <= m; i++) {
                //   B(I,J) = TEMP*B(I,J)
                let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                b.r[coorBJ + i] = re;
                b.i[coorBJ + i] = im;
            }
            for (let k = j + 1; k <= n; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    tempRe = alpha.re * a.r[coorAJ + k] - alpha.im * a.i[coorAJ + k];
                    tempIm = alpha.re * a.i[coorAJ + k] + alpha.im * a.r[coorAJ + k];
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) + TEMP*B(I,K)
                        let re = tempRe * b.r[coorBK + i] - tempIm * b.i[coorBK + i];
                        let im = tempRe * b.i[coorBK + i] + tempIm * b.r[coorBK + i];
                        b.r[coorBJ + i] += re;
                        b.i[coorBJ + i] += im;
                    }
                }
            }//k
        }
    }
}
