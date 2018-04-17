import {
    div_rxr,
    FortranArrEComplex,
    MatrixEComplex,
    mul_rxr
} from '../../../f_func'

export function normalLower(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number
) {

    let jx = kx;
    for (let j = 1; j <= n; j++) {
        const coorAJ = a.colOfEx(j);
        const isXZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!isXZero) {
            if (nounit) {
                //(a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    a.r[coorAJ + j],
                    a.i[coorAJ + j]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = jx;
            for (let i = j + 1; i <= n; i++) {
                ix += incx;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    a.r[coorAJ + i],
                    a.i[coorAJ + i]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
            }
        }
        jx += incx;
    }
}
