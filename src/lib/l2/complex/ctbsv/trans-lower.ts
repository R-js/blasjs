import { div_rxr, FortranArrEComplex, MatrixEComplex, mul_rxr } from '../../../f_func';

const { min } = Math;
export function transLower(
    kx: number,
    x: FortranArrEComplex,
    incx: number,
    a: MatrixEComplex,
    noconj: boolean,
    nounit: boolean,
    n: number,
    k: number): void {


    
    kx += (n - 1) * incx;
    let jx = kx;

    for (let j = n; j >= 1; j--) {
        let tempRe = x.r[jx - x.base];
        let tempIm = x.i[jx - x.base];
        let ix = kx;
        let L = 1 - j;

        //
        const extrI = min(n, j + k);
        const coorAJ = a.colOfEx(j);
        const sign = noconj ? 1 : -1;
        for (let i = extrI; i >= j + 1; i--) {
            //TEMP = TEMP - DCONJG(A(L+I,J))*X(IX)
            const { re, im } = mul_rxr(
                a.r[coorAJ + L + i],
                sign * a.i[coorAJ + L + i],
                x.r[ix - x.base],
                x.i[ix - x.base]
            );
            tempRe -= re;
            tempIm -= im;
            ix -= incx;
        }
        if (nounit) {
            const { re, im } = div_rxr(
                tempRe,
                tempIm,
                a.r[coorAJ + 1],
                sign * a.i[coorAJ + 1]
            );
            tempRe = re;
            tempIm = im;
        }
        x.r[jx - x.base] = tempRe;
        x.i[jx - x.base] = tempIm;
        jx -= incx;
        if ((n - j) >= k) {
            kx -= incx;
        }
    }
}
