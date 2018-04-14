import { div_rxr, FortranArrEComplex, MatrixEComplex, mul_rxr } from '../../../f_func';

const { max } = Math;

export function transUpper(
    noconj: boolean,
    x: FortranArrEComplex,
    incx: number,
    kx: number,
    a: MatrixEComplex,
    nounit: boolean,
    n: number,
    k: number) {

    let kplus1 = k + 1;
    let jx = kx;

    for (let j = 1; j <= n; j++) {
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        let ix = kx;
        let L = kplus1 - j;
        const extrI = max(1, j - k);
        const coorAJ = a.colOfEx(j);

        for (let i = extrI; i <= j - 1; i++) {
            const { re, im } = mul_rxr(
                a.r[coorAJ + L + i],
                (noconj ? a.i[coorAJ + L + i] : -a.i[coorAJ + L + i]),
                x.r[ix - x.base],
                x.i[ix - x.base]
            );
            tempRe -= re;
            tempIm -= im;
            ix += incx;
        }
        if (nounit) {
            const { re, im } = div_rxr(
                tempRe,
                tempIm,
                a.r[coorAJ + kplus1],
                (noconj ? a.i[coorAJ + kplus1] : -a.i[coorAJ + kplus1])
            );
            tempRe = re;
            tempIm = im;
        }
        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx += incx;
        if (j > k) {
            kx += incx;
        }
    }
}
