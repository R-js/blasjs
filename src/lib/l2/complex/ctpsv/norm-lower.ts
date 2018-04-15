import {
    div_rxr,
    FortranArrEComplex,
    mul_rxr
} from '../../../f_func';

export function normLower(
    kx: number,
    noconj: boolean,
    nounit: boolean,
    x: FortranArrEComplex,
    incx: number,
    ap: FortranArrEComplex,
    n: number
) {
    let kk = 1;
    let jx = kx;
    for (let j = 1; j <= n; j++) {
        const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
        if (!xIsZero) {
            if (nounit) {
                // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                const apkk = kk - ap.base;
                const { re, im } = div_rxr(
                    x.r[jx - x.base],
                    x.i[jx - x.base],
                    ap.r[apkk],
                    ap.i[apkk]
                );
                //console.log(`jx:${jx},x(jx)=(${x.r[jx - x.base]},${x.i[jx - x.base]})x/p=(${re},${im})`);
                x.r[jx - x.base] = re;
                x.i[jx - x.base] = im;

            }
            let tempRe = x.r[jx - x.base];
            let tempIm = x.i[jx - x.base];
            let ix = jx;
            for (let k = kk + 1; k <= kk + n - j; k++) {
                ix += incx;
                const apk = k - ap.base;
                const { re, im } = mul_rxr(
                    tempRe,
                    tempIm,
                    ap.r[apk],
                    ap.i[apk]
                );
                x.r[ix - x.base] -= re;
                x.i[ix - x.base] -= im;
            }
        }
        jx += incx;
        kk += (n - j + 1);
    }
}
