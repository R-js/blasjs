import {
    FortranArrEComplex,
    MatrixEComplex
} from '../../../f_func';

export function transLower(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number
): void {

    let jx = kx - x.base;
    for (let j = 1; j <= n; j++) {
        let tempRe = x.r[jx];
        let tempIm = x.i[jx];
        let ix = jx;
        const coords = a.colOfEx(j);
        if (noconj) {
            if (nounit) {
                //(a+ib)*(c+id) = (ac-bd)+i(ad+bc)
                const tr = tempRe * a.r[coords + j] - tempIm * a.i[coords + j];
                const ti = tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                tempRe = tr;
                tempIm = ti;
            }
            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                tempRe += a.r[coords + i] * x.r[ix] - a.i[coords + i] * x.i[ix];
                tempIm += a.r[coords + i] * x.i[ix] + a.i[coords + i] * x.r[ix];
            }
        } else {
            if (nounit) {
                //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
                const tr = tempRe * a.r[coords + j] + tempIm * a.i[coords + j];
                const ti = -tempRe * a.i[coords + j] + tempIm * a.r[coords + j];
                tempRe = tr;
                tempIm = ti;
            }
            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                tempRe += a.r[coords + i] * x.r[ix] + a.i[coords + i] * x.i[ix];
                tempIm += a.r[coords + i] * x.i[ix] - a.i[coords + i] * x.r[ix];
            }
        }
        x.r[jx] = tempRe;
        x.i[jx] = tempIm;
        jx += incx;
    }
}
