import { div_rxr, FortranArrEComplex, isZeroE, MatrixEComplex, mul_rxr } from '../../../f_func';
const { max } = Math;

export function normUpper(
    kx: number,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    noconj: boolean,
    nounit: boolean,
    n: number,
    k: number) {


    let kplus1 = k + 1;
    kx += (n - 1) * incx;
    let jx = kx;

    for (let j = n; j >= 1; j--) {
        const coorAJ = a.colOfEx(j);
        kx -= incx;
        const isXZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        const extrI = max(1, j - k);
        if (!isXZero) {
            let ix = kx;
            let L = kplus1 - j;
            if (nounit) {
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[kplus1 + coorAJ],
                    a.i[kplus1 + coorAJ]
                )
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];

            for (let i = j - 1; i >= extrI; i--) {
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    a.r[coorAJ + i + L],
                    a.i[coorAJ + i + L]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
                ix -= incx;
            }
        }
        jx -= incx;
    }
}
