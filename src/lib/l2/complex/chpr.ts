/*
    -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/


import { errMissingIm, errWrongArg, FortranArr } from '../../f_func';


export function chpr(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    ap: FortranArr, ): void {


    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0X20);


    let info = 0;
    if (ul !== 'u' && ul !== 'l') {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chpr', info));
    }

    if (n === 0 || alpha === 0) return;

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    let kk = 1;

    if (ul === 'u') {
        //Form  A  when upper triangle is stored in AP.
        let jx = kx - x.base;
        for (let j = 1; j <= n; j++) {
            if (!(x.r[jx] === 0 && x.i[jx] === 0)) {
                let tempRe = alpha * x.r[jx];
                let tempIm = alpha * -x.i[jx];
                let ix = kx - x.base;
                for (let k = kk; k <= kk + j - 2; k++) {
                    let apk = k - ap.base;
                    ap.r[apk] += x.r[ix] * tempRe - x.i[ix] * tempIm;
                    ap.i[apk] += x.i[ix] * tempIm + x.i[ix] * tempRe
                    ix += incx;
                }
                ap.i[kk + j - 1 - ap.base] = 0;
                ap.r[kk + j - 1 - ap.base] += x.r[jx] * tempRe - x.i[jx] * tempIm;
            } else {
                ap.i[kk + j - 1 - ap.base] = 0;
            }
            jx += incx;
            kk += j;
        }
    }
    else {
        // Form  A  when lower triangle is stored in AP.
        let jx = kx - x.base;
        for (let j = 1; j <= n; j++) {
            if (!(x.r[jx] === 0 && x.i[jx] === 0)) {
                let tempRe = alpha * x.r[jx];
                let tempIm = -alpha * x.i[jx];
                //
                ap.i[kk - ap.base] = 0;
                ap.r[kk - ap.base] += tempRe * x.r[jx] - tempIm * x.i[jx];
                let ix = jx;
                for (let k = kk + 1; k <= kk + n - j; k--) {
                    ix += incx;
                    ap.r[k - ap.base] += x.r[ix] * tempRe - x.i[ix] * tempIm;
                    ap.i[k - ap.base] += x.r[ix] * tempIm + x.i[ix] * tempRe;
                }
            }
            else {
                ap.i[kk - ap.base] = 0;
            }
            jx += incx;
            kk += n - j + 1;
        }
    }
}
