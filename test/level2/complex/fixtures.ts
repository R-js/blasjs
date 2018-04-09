import { complex, fortranArrComplex64 as arr64 } from '../../../src/lib/f_func';
import { bandmatrix_nxm_ku_kl, vector } from './matrices';

const pI = Infinity;
const nI = -Infinity;
const { PI, sin, cos, abs, sqrt } = Math;

const cospi = x => cos(PI * x);
const sinpi = x => sin(PI * x);

export const fixture = {
    // CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    cgbmv: {
        case0: {
            desc: 'y = alpha*A*x + beta*y, m=6,n=9,kl=1, ku=1, alpha(1.5,0.25)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                kl: 3,
                ku: 2,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: complex(2.5, +0.5),
                alpha: complex(0, 0),
                //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: bandmatrix_nxm_ku_kl(8, 6, 6, 3, 2),
                x: vector(8),
                y: vector(6),
            },
            expect: {
                y: [
                    complex(-2.0085962414741516, 1.6058503985404968),
                    complex(-0.18816410377621651, -0.60837343707680702),
                    complex(1.8727443814277649, -0.72995787858963013),
                    complex(2.9619127362966537, -0.49696569144725800),
                    complex(-0.13906472176313400, 2.5643529072403908),
                    complex(-0.15649497509002686, -0.74832186102867126)
                ]
            },
        },
    },
    // CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    /*
     *trans!='ntc'
     *m < 0
     *n < 0
     *kl < 0
     *ku < 0
     *lda < (kl + ku + 1)
     incx === 0
     incy === 0
    */
    cgbmvErrors: {
        case0: {
            desc: 'trans != ("n","t","c")',
            input: {
                trans: 'x',
                m: 6,
                n: 8,
                kl: 4,
                ku: 4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case1: {
            desc: 'm<0',
            input: {
                trans: 't',
                m: -6,
                n: 8,
                kl: 4,
                ku: 4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case2: {
            desc: 'sgbmv, m<0',
            input: {
                trans: 't',
                m: 6,
                n: -8,
                kl: 4,
                ku: 4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case3: {
            desc: 'kl<0',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: -4,
                ku: 4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            },
        },
        case4: {
            desc: 'ku<0',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: 4,
                ku: -4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case5: {
            desc: 'lda < (kl + ku + 1)',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: 4,
                ku: 4,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case6: {
            desc: 'incx=0',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: 4,
                ku: 4,
                lda: 10,
                incx: 0,
                incy: 1,
                beta: 2.5,
                alpha: 1.5,
            }
        },
        case7: {
            desc: 'incy=0',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: 4,
                ku: 4,
                lda: 10,
                incx: 1,
                incy: 0,
                beta: 2.5,
                alpha: 1.5,
            }
        },
    },

}
