
import {
    div_rxr,
    FortranArrEComplex,
    MatrixEComplex,
    mul_rxr
} from '../../../f_func';

export function transLower(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    n: number
) {
    kx += (n - 1) * incx;
    let jx = kx;
    for (let j = n; j >= 1; j--) {
        let ix = kx;
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        const coorAJ = a.colOfEx(j);
        for (let i = n; i >= j + 1; i--) {
            //TEMP = TEMP - CONJG(A(I,J))*X(IX)
            const { re, im } = mul_rxr(
                a.r[coorAJ + i],
                (noconj ? a.i[coorAJ + i] : - a.i[coorAJ + i]),
                x.r[ix - x.base],
                x.i[ix - x.base]
            );
            tempRe -= re;
            tempIm -= im;
            ix -= incx;
        }
        if (nounit) {
            // (a+ib)/(c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
            const { re, im } = div_rxr(
                tempRe,
                tempIm,
                a.r[coorAJ + j],
                (noconj ? a.i[coorAJ + j] : -a.i[coorAJ + j])
            )
            tempRe = re;
            tempIm = im;
        }
        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx -= incx;
    }
}
