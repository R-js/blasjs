import { FortranArrEComplex, mul_rxr } from '../../../f_func';

export function normLower(
    kx: number,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {

   
    let kk = n * (n + 1) / 2;
    kx = kx + (n - 1) * incx;
    let jx = kx;
    for (let j = n; j >= 1; j--) {
        const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!xIsZero) {
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = kx;
            for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                const apk = k - ap.base;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    ap.r[apk],
                    ap.i[apk]
                );
                x.r[ix - x.base] += re;
                x.i[ix - x.base] += im;
                ix -= incx;
            }
            if (nounit) {
                const apk2 = kk - n + j - ap.base;
                const { re, im } = mul_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    ap.r[apk2],
                    ap.i[apk2]
                );
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;
            }
        }
        jx -= incx;
        kk -= (n - j + 1);
    }
}
