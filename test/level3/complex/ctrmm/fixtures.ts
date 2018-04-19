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
                a: (() => {
                    const m = matrix_mxn(6, 6); //m*
                    m.r[0] = 0;
                    m.i[0] = 0;
                    return m;
                })(), // only uses 4x4!!
                lda: 6, //NBupper 4x4 of A is referenced
                b: (() => {
                    const m = matrix_mxn(6, 6); //m*
                    m.r[0] = 0;
                    m.i[0] = 0;
                    return m;
                })(),
                ldb: 4
            },
            expect: {
                b: [
                    complex(-0.41405065468799945, -1.5191217029575041),
                    complex(-0.69490876638628207, -1.4153085671446701),
                    complex(-0.83419064198640946, -0.50555104393530281),
                    complex(0.34645774818088582, 0.28057501053159395),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-0.94240802378391497, 1.8194772369531083),
                    complex(-2.9153658882833628, -2.3126831220695263),
                    complex(-1.5721870875427124, 0.75999293033223214),
                    complex(0.79851111695362287, 0.41539800883503275),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.30059598777530616, -4.7622991935580838E-004),
                    complex(-0.40030704404913392, 0.64982639160988953),
                    complex(0.33213353316073702, 0.20926732060823242),
                    complex(1.2837599485019321E-002, -0.19049002371540899),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.76325512471628953, 0.79546193461507575),
                    complex(-0.83792405458526675, 0.51076886233492236),
                    complex(-3.6959496736776953E-002, -5.4805674326484782E-002),
                    complex(0.12676149970562242, 6.4043945950648051E-002),
                    complex(0.13333636522293091, -2.2239003181457520),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-1.6508024981725062, -0.64423594455874311),
                    complex(0.37551185028990675, 7.4690390701194542E-002),
                    complex(-0.10496880817482218, 0.15546792320308533),
                    complex(-0.11660175095677679, -0.20949314958488327),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148E-002, 0.24226348102092743),
                    complex(0.67662846150460598, 0.40970986840795259),
                    complex(0.18500013285725259, 1.1724732013706247),
                    complex(0.51842448361866889, 1.3901729555213557E-002),
                    complex(-0.15617033770313232, -0.15971929721539613),
                    complex(0.72675073146820068, 1.9156390801072121E-002),
                    complex(1.1519117355346680, 0.25733837485313416),
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
