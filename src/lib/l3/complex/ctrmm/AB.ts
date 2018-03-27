import { Complex, errMissingIm, Matrix } from '../../../f_func';


export function AB(nounit: boolean, upper: boolean, n: number, m: number, a: Matrix, b: Matrix, alpha: Complex): void {

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }
    if (b.i === undefined) {
        throw new Error(errMissingIm('b.i'));
    }

    if (upper) {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            for (let k = 1; k <= m; k++) {
                const coorAK = a.colOfEx(k);
                let tempRe = alpha.re * b.r[coorBJ + k] - alpha.im * b.i[coorBJ + k];
                let tempIm = alpha.re * b.i[coorBJ + k] + alpha.im * b.r[coorBJ + k];
                for (let i = 1; i <= k - 1; i++) {
                    b.r[coorBJ + i] += tempRe * a.r[coorAK + i] - tempIm * a.i[coorAK + i];
                    b.i[coorBJ + i] += tempRe * a.i[coorAK + i] + tempIm * a.i[coorAK + i];
                }
                if (nounit) {
                    let re = tempRe * a.r[coorAK + k] - tempIm * a.i[coorAK + k];
                    let im = tempRe * a.i[coorAK + k] + tempIm * a.r[coorAK + k];
                    tempRe = re;
                    tempIm = im;
                }
                b.r[coorBJ + k] = tempRe;
                b.i[coorBJ + k] = tempIm;
            }
        }
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
