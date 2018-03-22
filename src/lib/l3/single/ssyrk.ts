/*
*> -- Jacob Bogers, 03/2018, JS port, jkfbogers@gmail.com
*> -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*>
*/

/*
*> SSYRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*> in the second case.
*/


import { errWrongArg, Matrix } from '../../f_func';

const { max } = Math;

export function ssyrk(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    n: number,
    k: number,
    alpha: number,
    a: Matrix,
    lda: number,
    beta: number,
    c: Matrix,
    ldc: number): void {

    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0X20);
    const tr = String.fromCharCode(trans.charCodeAt(0) | 0X20);

    const nrowA = (tr === 'n') ? n : k;

    //Test the input parameters.

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (!'ntc'.includes(tr)) {
        info = 2;
    }
    else if (n < 0) {
        info = 3;
    }
    else if (k < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (ldc < max(1, n)) {
        info = 10;
    }

    if (info) {
        throw new Error(errWrongArg('ssyrk', info));
    }

    //*     Quick return if possible.

    if (n === 0 || alpha === 0 || (k === 0 && beta === 1)) return;

    //*      And when  alpha.eq.zero.
    if (alpha === 0) {

        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? n : j;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = beta ? beta * c.r[coorCJ + i] : 0;
            }
        }
        return;
    }

    //*     Start the operations.

    if (tr === 'n') {
        //Form  C := alpha*A*A**T + beta*C.

        for (let j = 1; j <= n; j++) {
            const coorCJ = c.colOfEx(j);
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            for (let i = start; i <= stop; i++) {
                c.r[coorCJ + i] = beta === 0 ? 0 : (beta !== 1 ? beta * c.r[coorCJ + i] : c.r[coorCJ + i]);
            }
            for (let l = 1; l <= k; l++) {
                const coorAL = a.colOfEx(l);
                if (a.r[coorAL + j] !== 0) {
                    let temp = alpha * a.r[coorAL + j];
                    for (let i = start; i <= stop; i++) {
                        c.r[coorCJ + i] += temp * a.r[coorAL + i];
                    }
                }
            }
        }
    }
    else {
        // Form  C := alpha*A**T*A + beta*C.
        for (let j = 1; j <= n; j++) {
            const start = ul === 'u' ? 1 : j;
            const stop = ul === 'u' ? j : n;
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            for (let i = start; i <= stop; i++) {
                let temp = 0;
                const coorAI = a.colOfEx(i);
                for (let l = 1; l <= k; l++) {
                    temp += a.r[coorAI + l] * a.r[coorAJ + l];
                }
                c.r[coorCJ + i] = alpha * temp + (beta === 0 ? 0 : (beta !== 1 ? beta * c.r[coorCJ + i] : c.r[coorCJ + i]))
            }
        }
    }
}
