
import {
    div_rxr,
    FortranArrEComplex,
    MatrixEComplex,
    mul_rxr
} from '../../../f_func';

export function transUpper(
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
        let ix = kx;
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        const coords = a.colOfEx(j);

        for (let i = 1; i <= j - 1; i++) {
            const { re, im } = mul_rxr(
                a.r[coords + i],
                (noconj ? a.i[coords + i] : -a.i[coords + i]),
                x.r[ix - x.base],
                x.i[ix - x.base]
            );
            tempRe -= re;
            tempIm -= im;
            ix += incx;
        }
        // IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
        if (nounit) {
            const { re, im } = div_rxr(
                tempRe,
                tempIm,
                a.r[coords + j],
                (noconj ? a.i[coords + j] : -a.i[coords + j])
            );
            tempRe = re;
            tempIm = im;
        }

        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx += incx;
    }
}
