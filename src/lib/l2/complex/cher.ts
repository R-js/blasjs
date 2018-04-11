/*  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*/

import { errMissingIm, errWrongArg, FortranArr, lowerChar, Matrix } from '../../f_func';

const { max } = Math;

export function cher(
    uplo: 'u' | 'l',
    n: number,
    alpha: number, //this is the ONLY level 2, complex function that has a REAL/DOUBLE PREC alpha  
    x: FortranArr,
    incx: number,
    a: Matrix,
    lda: number): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }
    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(ul)) {
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
        let jx = kx;
        for (let j = 1; j <= n; j++) {
            // console.log(`(${x.r[jx - x.base]},${-x.i[jx - x.base]})`);

            const coords = a.colOfEx(j);
            //console.log(`j=(${j})`);
            const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
            if (!xIsZero) {

                let tempRe = alpha * x.r[jx - x.base];
                let tempIm = alpha * -x.i[jx - x.base];

                let ix = kx;

                for (let i = 1; i <= j - 1; i++) {
                    //  console.log(`(${x.r[ix - x.base]},${x.i[ix - x.base]})*(${tempRe},${tempIm})`);
                    const re = x.r[ix - x.base] * tempRe - x.i[ix - x.base] * tempIm;
                    const im = x.r[ix - x.base] * tempIm + x.i[ix - x.base] * tempRe;
                    a.r[coords + i] += re;
                    a.i[coords + i] += im;
                    //console.log(`${i},${j},(${a.r[coords + i]},${a.i[coords + i]})`);

                    //console.log(`(i,j)=(${i},${j}), X*temp=(${re},${im})`);

                    ix += incx;
                }
                a.i[coords + j] = 0;
                a.r[coords + j] += (x.r[jx - x.base] * tempRe - x.i[jx - x.base] * tempIm);
            }
            else {
                a.i[coords + j] = 0;
            }
            jx += incx;
        }
    }
    else {
        //    * Form  A  when A is stored in lower triangle.
        let jx = kx
        for (let j = 1; j <= n; j++) {
            const coords = a.colOfEx(j);
            if (!(x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0)) {
                //   TEMP = ALPHA*CONJG(X(JX))
                let tempRe = alpha * x.r[jx - x.base];
                let tempIm = -alpha * x.i[jx - x.base];

                //A(J,J) = REAL(A(J,J)) + REAL(TEMP*X(JX))
                a.i[coords + j] = 0;
                a.r[coords + j] += tempRe * x.r[jx - x.base] - tempIm * x.i[jx - x.base];

                //IX = JX

                let ix = jx;
                for (let i = j + 1; i <= n; i++) {
                    ix += incx;
                    a.r[coords + i] += x.r[ix - x.base] * tempRe - x.i[ix - x.base] * tempIm;
                    a.i[coords + i] += x.r[ix - x.base] * tempIm + x.i[ix - x.base] * tempRe;

                }
            }
            else {
                a.i[coords + j] = 0;
            }
            jx += incx;
        }
    }
}
