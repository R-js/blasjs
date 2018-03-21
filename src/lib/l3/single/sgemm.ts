import { errWrongArg, Matrix2D } from '../../f_func';

/*  -- Jacob Bogers, 03/2008, JS port, jkfbogers@gmail.cmom
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*/

const { max } = Math;

export function sgemm(
    transA: 'n' | 't' | 'c',
    transB: 'n' | 't' | 'c',
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: Matrix2D,
    lda: number,
    b: Matrix2D,
    ldb: number,
    beta: number,
    c: Matrix2D,
    ldc: number): void {

    // faster then String.toLowerCase()
    const trA: 'n' | 't' | 'c' = String.fromCharCode(transA.charCodeAt(0) | 0X20) as any;
    const trB: 'n' | 't' | 'c' = String.fromCharCode(transB.charCodeAt(0) | 0X20) as any;

    const notA = trA === 'n';
    const notB = trB === 'n';

    const nrowA = notA ? m : k;
    //ncolA is never used, I checked in the original F-code
    // also checked online http://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html#gafe51bacb54592ff5de056acabd83c260
    // ncolA is not used
    //const ncolA = notA ? k : n;

    const nrowB = notB ? k : n;

    let info = 0;
    if (!notA && trA !== 'c' && trA !== 't') {
        info = 1;
    }
    else if (!notB && trB !== 'c' && trB !== 't') {
        info = 2;
    }
    else if (m < 0) {
        info = 3;
    }
    else if (n < 0) {
        info = 4;
    }
    else if (k < 0) {
        info = 5;
    }
    else if (lda < max(1, nrowA)) {
        info = 8;
    }
    else if (ldb < max(1, nrowB)) {
        info = 10;
    }
    else if (ldb < max(1, m)) {
        info = 13;
    }
    // ok?
    if (info !== 0) {
        throw new Error(errWrongArg('sgemm', info));
    }

    //*     Quick return if possible.
    if (m === 0 || n === 0 ||
        (
            (alpha === 0 || k === 0) && beta === 1
        )
    ) return;


    if (alpha === 0) {
        if (beta === 0) {
            c.r.fill(0); //fast
        }
        else {
            // I have to typecast it this way for TS compiler not to nag!!
            (c.r as Float32Array).set((c.r as Float32Array).map(v => v * beta), 0);
        }
        return;
    }

    //*     Start the operations.
    if (notB) {
        if (notA) {
            //Form  C := alpha*A*B + beta*C.
            for (let j = 1; j <= n; j++) {
                if (beta === 0) {
                    c.r.fill(0);
                }
                else if (beta !== 1) {
                    (c.r as Float32Array).set((c.r as Float32Array).map(v => v * beta), 0);
                }
                const coorB = b.colOfEx(j);
                const coorC = c.colOfEx(j);
                for (let L = 1; L <= k; L++) {
                    let temp = alpha * b.r[coorB + L];
                    const coorA = a.colOfEx(L);
                    for (let i = 1; i <= m; i++) {
                        c.r[coorC + i] += temp * a.r[coorA + i]
                    }
                }
            }

        }
        else {
            //     Form  C := alpha*A**T*B + beta*C
            for (let j = 1; j <= n; j++) {
                const coorB = b.colOfEx(j);
                const coorC = c.colOfEx(j);
                for (let i = 1; i <= m; i++) {
                    let temp = 0;
                    const coorA = a.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        temp += a.r[L + coorA] * b.r[coorB + L];
                    }
                    if (beta === 0) {
                        c.r[coorC + i] = alpha * temp;
                    }
                    else {
                        c.r[coorC + i] = alpha * temp + beta * c.r[coorC + i];
                    }
                }
            }
        }
    }
    else {
        if (notA) {
            // Form  C := alpha*A*B**T + beta*C
            for (let j = 1; j <= n; j++) {
                if (beta === 0) {
                    c.r.fill(0);
                } else if (beta !== 1) {
                    (c.r as Float32Array).set((c.r as Float32Array).map(v => v * beta), 0);
                }
                const coorC = c.colOfEx(j);
                for (let L = 1; L <= k; L++) {
                    let temp = alpha * b.r[(L - b.colBase) * b.nrRows - b.rowBase + j];
                    const coorA = a.colOfEx(L);
                    for (let i = 1; i <= m; i++) {
                        c.r[coorC + i] += temp * a.r[coorA + i];
                    }
                }
            }
        }
        else {
            //  Form  C := alpha*A**T*B**T + beta*C
            for (let j = 1; j <= n; j++) {
                for (let i = 1; i <= m; i++) {
                    let temp = 0;
                    const coorA = a.colOfEx(i);
                    for (let L = 1; L <= k; L++) {
                        //TODO: write out colOfEx
                        temp += a.r[coorA + L] * b.r[b.colOfEx(L) + j];
                    }

                    const coorC = c.colOfEx(j);
                    if (beta === 0) {
                        c.r[coorC + i] = alpha * temp;
                    }
                    else {
                        c.r[coorC + i] = alpha * temp + beta * c.r[coorC + i];
                    }
                }
            }
        }
    }
}
