/*
*>  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/
import { errMissingIm, errWrongArg, FortranArr } from '../../f_func';
/*
*>
*> CTPSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix, supplied in packed form.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*/



export function ctpsv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    ap: FortranArr,
    x: FortranArr,
    incx: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('ap.i'));
    }


    // faster then String.toLowerCase()
    const ul: 'u' | 'l' = String.fromCharCode(uplo.charCodeAt(0) | 0X20) as any;
    const tr: 'n' | 't' | 'c' = String.fromCharCode(trans.charCodeAt(0) | 0X20) as any;
    const dg: 'u' | 'n' = String.fromCharCode(diag.charCodeAt(0) | 0X20) as any;

    let info = 0;
    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (tr !== 'n' && tr !== 't' && tr !== 'c') {
        info = 2;
    }
    else if (dg !== 'u' && dg !== 'n') {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (incx === 0) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (tr === 'n') {
        // Form  x := inv( A )*x.
        if (ul === 'u') {
            let kk = n * (n - 1) / 2;
            let jx = kx + (n - 1) * incx - x.base;
            for (let j = n; j >= 1; j--) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (xIsZero) {
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk - ap.base;
                        let xa = x.r[jx];
                        let xb = x.i[jx];
                        let xc = ap.r[apkk];
                        let xd = ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        x.r[jx] = (xa * xc + xb * xd) / DOM;
                        x.i[jx] = (-xa * xd + xb * xc) / DOM;
                    }
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = jx;
                    for (let k = kk - 1; k >= kk - j + 1; k--) {
                        ix -= incx;
                        const apk = k - ap.base;
                        x.r[ix] -= tempRe * ap.r[apk] - tempIm * ap.i[apk];
                        x.i[ix] -= tempRe * ap.r[apk] + tempIm * ap.i[apk];
                    }
                }
                jx -= incx;
                kk -= j;
            }
        }
        else {
            let kk = 1;
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                const xIsZero = x.r[jx] === 0 && x.i[jx] === 0;
                if (!xIsZero) {
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk - ap.base;
                        let xa = x.r[jx];
                        let xb = x.i[jx];
                        let xc = ap.r[apkk];
                        let xd = ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        x.r[jx] = (xa * xc + xb * xd) / DOM;
                        x.i[jx] = (-xa * xd + xb * xc) / DOM;
                    }
                    let tempRe = x.r[jx];
                    let tempIm = x.i[jx];
                    let ix = jx;
                    for (let k = kk + 1; k <= kk + n - j; k++) {
                        ix += incx;
                        const apk = k - ap.base;
                        x.r[ix] -= tempRe * ap.r[apk] - tempIm * ap.i[apk];
                        x.i[ix] -= tempRe * ap.i[apk] + tempIm * ap.r[apk];
                    }
                }
                jx += incx;
                kk += (n - j + 1);
            }
        }
    }
    else {
        // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
        if (ul === 'u') {
            let kk = 1;
            let jx = kx - x.base;
            for (let j = 1; j <= n; j++) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = kx - x.base;
                if (noconj) {
                    for (let k = kk; k <= kk + j - 2; k++) {
                        const apk = k - ap.base;
                        tempRe -= ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm -= ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk + j - 1 - ap.base;
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = ap.r[apkk];
                        let xd = ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                else {
                    for (let k = kk; k <= kk + j - 2; k--) {
                        //(a-ib)*(c+id)=(ac+bd)+i(ad-bc)
                        const apk = k - ap.base;
                        tempRe -= ap.r[apk] * x.r[ix] + ap.i[apk] * x.i[ix];
                        tempIm -= ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                        ix += incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk + j - 1 - ap.base;
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = ap.r[apkk];
                        //conj
                        let xd = -ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx += incx;
                kk += j;
            }
        }
        else {
            let kk = (n * (n + 1)) / 2;
            kx += (n - 1) * incx;
            let jx = kx - x.base;
            for (let j = n; j >= 1; j--) {
                let tempRe = x.r[jx];
                let tempIm = x.i[jx];
                let ix = kx - x.base;
                if (noconj) {
                    for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                        const apk = k - ap.base;
                        tempRe -= ap.r[apk] * x.r[ix] - ap.i[apk] * x.i[ix];
                        tempIm -= ap.r[apk] * x.i[ix] + ap.i[apk] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk - n - j - ap.base;
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = ap.r[apkk];
                        let xd = ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                else {
                    for (let k = kk; k >= kk - (n - (j + 1)); k--) {
                        const apk = k - ap.base;
                        //(a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                        tempRe -= ap.r[apk] * x.r[ix] + ap.i[apk] * x.i[ix];
                        tempIm -= ap.r[apk] * x.i[ix] - ap.i[apk] * x.r[ix];
                        ix -= incx;
                    }
                    if (nounit) {
                        // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
                        const apkk = kk - n + j - ap.base;
                        let xa = tempRe;
                        let xb = tempIm;
                        let xc = ap.r[apkk];
                        //conj
                        let xd = -ap.i[apkk];

                        let DOM = xc * xc + xd * xd;

                        tempRe = (xa * xc + xb * xd) / DOM;
                        tempIm = (-xa * xd + xb * xc) / DOM;
                    }
                }
                x.r[jx] = tempRe;
                x.i[jx] = tempIm;
                jx -= incx;
                kk -= (n - j + 1);
            }
        }
    }
}
