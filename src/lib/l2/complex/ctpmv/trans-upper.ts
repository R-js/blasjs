import { FortranArrEComplex, mul_rxr } from '../../../f_func';

export function transUpper(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {

    //Form  x := A**T*x  or  x := A**H*x.
    let kk = n * (n + 1) / 2;
    let jx = kx + (n - 1) * incx;
    for (let j = n; j >= 1; j--) {
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        let ix = jx;
        if (nounit) {
            //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
            const apk2 = kk - ap.base;
            const { re, im } = mul_rxr(
                tempRe,
                tempIm,
                ap.r[apk2],
                noconj ? ap.i[apk2] : -ap.i[apk2]
            );
            tempRe = re;
            tempIm = im;
        }
        for (let k = kk - 1; k >= kk - j + 1; k--) {
            ix -= incx;
            const apk = k - ap.base;
            //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
            const { re, im } = mul_rxr(
                ap.r[apk],
                noconj ? ap.i[apk] : -ap.i[apk],
                x.r[ix - x.base],
                x.i[ix - x.base]
            )
            tempRe += re;
            tempIm += im;
        }

        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx -= incx;
        kk -= j;
    }
}
