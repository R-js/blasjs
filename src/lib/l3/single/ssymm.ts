/*  -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*/

import { errWrongArg, Matrix } from '../../f_func';
/*
*>
*> SSYMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*A*B + beta*C,
*>
*> or
*>
*>    C := alpha*B*A + beta*C,
*>
*> where alpha and beta are scalars,  A is a symmetric matrix and  B and
*> C are  m by n matrices.
*> 
*/

const { max } = Math;

export function ssymm(
    side: 'l' | 'r',
    uplo: 'u' | 'l',
    m: number,
    n: number,
    alpha: number,
    a: Matrix,
    lda: number,
    b: Matrix,
    ldb: number,
    beta: number,
    c: Matrix,
    ldc: number): void {


    const si = String.fromCharCode(side.charCodeAt(0) | 0X20);
    const ul = String.fromCharCode(uplo.charCodeAt(0) | 0X20);

    const nrowA = si === 'l' ? m : n;

    //Test the input parameters.

    let info = 0;
    if (si !== 'l' && si !== 'r') {
        info = 1;
    }
    else if (ul !== 'u' && ul !== 'l') {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (lda < max(1, nrowA)) {
        info = 7;
    }
    else if (ldb < max(1, m)) {
        info = 9;
    }
    else if (ldc < max(1, m)) {
        info = 12;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ssymm', info));
    }
    //*     Quick return if possible.

    if (m === 0 || n === 0 ||
        (alpha === 0 && beta === 1)
    ) return;

    //when alpha is zero
    if (alpha === 0) {
        if (beta === 0) {
            for (let j = 1; j <= n; j++) {
                const coords = c.colOfEx(j);
                c.r.fill(0, coords + 1, coords + m + 1);
            }
        }
        else {
            for (let j = 1; j <= n; j++) {
                const coords = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    c.r[coords + i] = beta * c.r[coords + i];
                }
            }
        }
        return;
    }

    //   Start the operations.
    if (si === 'l') {
        //Form  C := alpha*A*B + beta*C.
        if (ul === 'u') {

        }
        else {

        }

    }
    else {
        //  Form  C := alpha*B*A + beta*C.
        for (let j = 1; j <= 1 - n; j++) {

            //pre-calc
            const coorAJ = a.colOfEx(j);
            const coorCJ = c.colOfEx(j);
            const coorBJ = b.colOfEx(j);

            let temp1 = alpha * a.r[coorAJ + j];
            if (beta === 0) {
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] = temp1 * b.r[coorBJ + i];
                }
            }
            else {
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] = beta * c.r[coorCJ + i] + temp1 * b.r[coorBJ + i];
                }
            }
            for (let k = 1; k <= j - 1; k++) {
                let temp1 = alpha * (ul === 'u' ? a.r[coorAJ + k] : a.r[a.colOfEx(j) + j]);
                const coorBK = b.colOfEx(k);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1 * b.r[coorBK + i];
                }
            }
            for (let k = j + 1; j <= n; j++) {
                let temp1 = alpha * (ul === 'u' ? a.r[a.colOfEx(k) + j] : a.r[coorAJ + k]);
                const coorBK = b.colOfEx(k);
                for (let i = 1; i <= m; i++) {
                    c.r[coorCJ + i] += temp1 * b.r[coorBK + i];
                }
            }
        }
    }
}
