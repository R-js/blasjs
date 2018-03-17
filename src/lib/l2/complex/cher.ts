/*  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

import { errMissingIm, errWrongArg, FortranArr, Matrix2D } from '../../f_func';

const { max } = Math;

export function cher(
    uplo: 'u' | 'l',
    n: number,
    alpha: number, //acc to normal pattern, this would be complex
    x: FortranArr,
    incx: number,
    a: Matrix2D,
    lda: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }
    if (a.i === undefined) {
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
        info = 3;
    }
    else if (lda < max(1, n)) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('cher', info));
    }



    if (n === 0 || alpha === 0) return;
    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    if (ul === 'u') {
        //*        Form  A  when A is stored in upper triangle.
        let jx = kx - x.base;
        for (let j = 1; j <= n; j++) {
            const coords = a.colOfEx(j);
            if (x.r[jx] !== 0 || x.i[jx] !== 0) {
                let tempRe = alpha * x.r[jx];
                let tempIm = alpha * x.i[jx];
                let ix = kx - x.base;

                for (let i = 1; i <= j - 1; i++) {
                    a.r[coords + i] += x.r[ix] * tempRe - x.i[ix] * tempIm;
                    a.i[coords + i] += x.r[ix] * tempIm + x.i[ix] * tempRe;
                    ix += incx;
                }
                a.i[coords + j] = 0;
                a.r[coords + j] += x.r[jx] * tempRe - x.i[jx] * tempIm;
            }
            else {
                a.i[coords + j] = 0;
            }
            jx += incx;
        }
    }
    else {
        //    * Form  A  when A is stored in lower triangle.
        let jx = kx - x.base;
        for (let j = 1; j <= n; j++) {
            const coords = a.colOfEx(j);
            if (!(x.r[jx] === 0 && x.i[jx] === 0)) {
                let tempRe = alpha * x.r[jx];
                let tempIm = -alpha * x.i[jx];
                a.i[coords + j] = 0;
                let ix = kx - x.base;
                for (let i = 1; i <= j - 1; i++) {
                    a.r[coords + i] += x.r[ix] * tempRe - x.i[ix] * tempIm;
                    a.i[coords + i] += x.r[ix] * tempIm + x.i[ix] * tempRe;
                    ix += incx;
                }
                a.r[coords + j] += x.r[jx] * tempRe - x.i[jx] * tempIm;
                a.i[coords + j] = 0;
            }
            else {
                a.i[coords + j] = 0;
            }
            jx += incx;
        }
    }
}
