import { div_rxr, FortranArrEComplex, isZeroE, MatrixEComplex, mul_rxr } from '../../../f_func';

const { min } = Math;

export function normLower(
    kx: number,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    noconj: boolean,
    nounit: boolean,
    n: number,
    k: number) {

    /*if (!x.i || !a.i) {
        return;
    }*/

    let jx = kx;
    for (let j = 1; j <= n; j++) {
        kx += incx;
        const xIsZero = isZeroE(x.r[jx - x.base], x.i[jx - x.base]);
        const extrI = min(n, j + k);
        if (!xIsZero) {
            let ix = kx;
            let L = 1 - j;
            const coorAJ = a.colOfEx(j);
            if (nounit) {
                // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[coorAJ + 1],
                    a.i[coorAJ + 1]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            for (let i = j + 1; i <= extrI; i++) {
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    a.r[coorAJ + L + i],
                    a.i[coorAJ + L + i]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
                ix += incx;
            }
        }
        jx += incx;
    }
}
