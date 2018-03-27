import { Complex, errMissingIm, Matrix } from '../../../f_func';
//Form  B := alpha*inv( A )*B.

export function invAB(
    nounit: boolean,
    upper: boolean,
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
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = m; k >= 1; k--) {
                const coorAK = a.colOfEx(k);
                const isBZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                if (!isBZero) {
                    if (nounit) {
                        //B(K,J) = B(K,J)/A(K,K)
                        // TEMP = ONE/A(J,J)
                        //(a+ib)/(c+id)
                        // re= (ac+bd)/(c*c+d*d)
                        // im =(bc-ad)/(c*c+d*d)
                        const _a = b.r[coorBJ + k];
                        const _b = b.i[coorBJ + k];
                        const _c = a.r[coorAK + k];
                        const _d = a.i[coorAK + k];
                        const n = _c * _c + _d * _d;

                        const re = (_a * _c + _b * _d) / n;
                        const im = (_b * _c - _a * _d) / n;

                        b.r[coorBJ + k] = re;
                        b.i[coorBJ + k] = im;
                    }
                    for (let i = 1; i <= k - 1; i++) {
                        // B(I,J) = B(I,J) - B(K,J)*A(I,K)
                        const _a = b.r[coorBJ + k];
                        const _b = b.i[coorBJ + k];
                        const _c = a.r[coorAK + i];
                        const _d = a.i[coorAK + i];
                        const n = _c * _c + _d * _d;

                        const re = (_a * _c + _b * _d) / n;
                        const im = (_b * _c - _a * _d) / n;
                        b.r[coorBJ + k] -= re;
                        b.i[coorBJ + k] -= im;
                    }
                }
            }//k
        }//j
    }//upper
    else {
        for (let j = 1; j <= n; j++) {
            const coorBJ = b.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = 1; k <= m; k++) {
                const bIsZero = b.r[coorBJ + k] === 0 && b.i[coorBJ + k] === 0;
                const coorAK = a.colOfEx(k);
                if (!bIsZero) {
                    if (nounit) {
                        //B(K,J) = B(K,J)/A(K,K)
                        // TEMP = ONE/A(J,J)
                        // (a+ib)/(c+id)
                        // re= (ac+bd)/(c*c+d*d)
                        // im =(bc-ad)/(c*c+d*d)
                        const _a = b.r[coorBJ + k];
                        const _b = b.i[coorBJ + k];
                        const _c = a.r[coorAK + k];
                        const _d = a.i[coorAK + k];
                        const n = _a * _a + _b * _b;
                        const re = (_a * _c + _b * _d) / n;
                        const im = (_b * _c - _a * _d) / n;
                        b.r[coorBJ + k] = re;
                        b.i[coorBJ + k] = im;
                    }
                    for (let i = k + 1; k <= m; i++) {
                        const re = b.r[coorBJ + k] * a.r[coorAK + i] - b.i[coorBJ + k] * a.i[coorAK + i];
                        const im = b.r[coorBJ + k] * a.i[coorAK + i] + b.i[coorBJ + k] * a.r[coorAK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + k] -= im;
                    }
                }//bNoZero
            }//k
        }//j
    }//upper
}
