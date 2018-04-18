import {
    complex,
    fortranArrComplex64 as arr64,
    fortranMatrixComplex64
} from '../../../../src/lib/f_func';

import {
    bandmatrix_nxm_ku_kl,
    diagonal_nxn,
    matrix_mxn,
    vector
} from '../../../matrices';

const pI = Infinity;
const nI = -Infinity;
const { PI, sin, cos, abs, sqrt } = Math;

const cospi = x => cos(PI * x);
const sinpi = x => sin(PI * x);


export const fixture = {
    // C := alpha*op( A )*op( B ) + beta*C,
    // CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    ctrmm: {
        case0: {
            desc: 'si="l", ul="u",trA="n", diag="n", alpha=(0.2,0.8)',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                a: matrix_mxn(6, 6), // only uses 4x4!!
                lda: 6, //NBupper 4x4 of A is referenced
                b: matrix_mxn(6, 6), //m*n
                ldb: 4
            },
            expect: {
                b: [
                    complex(0, 0)
                ]
            },
        },
    },
    ctrmmErrors: {
        case0: {
            desc: 'a has no imaginary part',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [0],
                b: [complex(0, 0)]
            }
        },
        case2: {
            desc: 'b has no imaginary part',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [0]
            }
        },
        case3: {
            desc: 'sid!="lr"',
            input: {
                side: 'x', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                side: 'l', //A*B
                uplo: 'x',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case5: {
            desc: 'transA',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: -1,
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case6: {
            desc: 'm<0',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: -1,
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 4,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case7: {
            desc: 'k<0',
            input: {
                trA: 'c',
                trB: 'c',
                m: -6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            }
        },
        case8: {
            desc: 'lda< max(1,nrowA)',
            input: {
                trA: 'n',
                trB: 'n',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 4, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            }
        },
        case9: {
            desc: 'ldb< max(1,nrowB)',
            input: {
                trA: 'n',
                trB: 'n',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 3, // ldb >= K
                ldc: 8, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            }
        },
        case10: {
            desc: 'ldc< max(1,nrowC)',
            input: {
                trA: 'n',
                trB: 'n',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 5, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            }
        },

    },
}
