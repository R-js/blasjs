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
            const coorAJ = a.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = 1; k <= j - 1; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) - A(K,J)*B(I,K)
                        let re = a.r[coorAJ + k] * b.r[coorBK + i] - a.i[coorAJ + k] * b.i[coorBK + i];
                        let im = a.r[coorAJ + k] * b.i[coorBK + i] + a.r[coorAJ + k] * b.r[coorBK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                }
            }
            if (nounit) {
                // TEMP = ONE/A(J,J)
                //(a+ib)/(c+id)
                // re= (ac+bd)/(c*c+d*d)
                // im =(bc-ad)/(c*c+d*d)
                // (1+i0)/(c+id), a=1,b=0
                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                let _c = a.r[coorAJ + j];
                let _d = a.i[coorAJ + j];
                n = _c * _c + _d * _d;
                let tempRe = _c / n;
                let tempIm = _d / n;
                for (let i = 1; i <= m; i++) {
                    let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                    let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }// nounit
        }//k
    }//upper
    else {
        for (let j = n; j >= 1; j--) {
            const coorBJ = b.colOfEx(j);
            const coorAJ = a.colOfEx(j);
            if (!alphaIsOne) {
                for (let i = 1; i <= m; i++) {
                    const re = alpha.re * b.r[coorBJ + i] - alpha.im * b.i[coorBJ + i];
                    const im = alpha.re * b.i[coorBJ + i] + alpha.im * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }
            for (let k = j + 1; k <= n; k++) {
                const coorBK = b.colOfEx(k);
                const aIsZero = a.r[coorAJ + k] === 0 && a.i[coorAJ + k] === 0;
                if (!aIsZero) {
                    for (let i = 1; i <= m; i++) {
                        // B(I,J) = B(I,J) - A(K,J)*B(I,K)
                        let re = a.r[coorAJ + k] * b.r[coorBK + i] - a.i[coorAJ + k] * b.i[coorBK + i];
                        let im = a.r[coorAJ + k] * b.i[coorBK + i] + a.r[coorAJ + k] * b.r[coorBK + i];
                        b.r[coorBJ + i] -= re;
                        b.i[coorBJ + i] -= im;
                    }
                }
            }
            if (nounit) {
                // TEMP = ONE/A(J,J)
                //(a+ib)/(c+id)
                // re= (ac+bd)/(c*c+d*d)
                // im =(bc-ad)/(c*c+d*d)
                // (1+i0)/(c+id), a=1,b=0
                // re= c/(cc+dd)
                // im =-d/(cc+dd)
                let _c = a.r[coorAJ + j];
                let _d = a.i[coorAJ + j];
                n = _c * _c + _d * _d;
                let tempRe = _c / n;
                let tempIm = _d / n;
                for (let i = 1; i <= m; i++) {
                    let re = tempRe * b.r[coorBJ + i] - tempIm * b.i[coorBJ + i];
                    let im = tempRe * b.i[coorBJ + i] + tempIm * b.r[coorBJ + i];
                    b.r[coorBJ + i] = re;
                    b.i[coorBJ + i] = im;
                }
            }// nounit
        }
    }
}


