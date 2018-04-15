import {
    div_rxr,
    FortranArrEComplex,
    mul_rxr
} from '../../../f_func';

export function transUpper(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {
    //console.log({ incx });

    let kk = 1;
    let jx = kx;

    for (let j = 1; j <= n; j++) {
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        let ix = kx;
        //
        for (let k = kk; k <= kk + j - 2; k++) {
            const { re, im } = mul_rxr(
                ap.r[k - ap.base],
                noconj ? ap.i[k - ap.base] : -ap.i[k - ap.base],
                x.r[ix - x.base],
                x.i[ix - x.base]
            );
            tempRe -= re;
            tempIm -= im;
            ix += incx;
        }
        //
        if (nounit) {
            const apkk = kk + j - 1 - ap.base;
            const { re, im } = div_rxr(
                tempRe,
                tempIm,
                ap.r[apkk],
                noconj ? ap.i[apkk] : -ap.i[apkk]
            );
            tempRe = re;
            tempIm = im;
        }
        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx += incx;
        kk += j;
    }
}
