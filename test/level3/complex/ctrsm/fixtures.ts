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
    ctrsm: {
        case0: {
            desc: '(trivial) m=0',
            input: {
                side: 'l',
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                ldb: 4,
                m: 0,
                n: 6,
                lda: 6,
                alpha: complex(0.2, 0.8),
                a: (() => {
                    const m = matrix_mxn(6, 6);
                    return m;
                })(),
                b: (() => {
                    const m = matrix_mxn(6, 6);
                    return m;
                })(),

            },
            expect: {
                b: (() => {
                    const m = matrix_mxn(6, 6);
                    return m.toArr();
                })(),
            },
        },
        case1: {
            desc: '(trivial) n=0',
            input: {
                side: 'l',
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                ldb: 4,
                m: 4,
                n: 0,
                lda: 6,
                alpha: complex(0.2, 0.8),
                a: (() => {
                    const m = matrix_mxn(6, 6);
                    //do some stuff here
                    return m;
                })(),
                b: (() => {
                    const m = matrix_mxn(6, 6);
                    //do some stuff here
                    return m;
                })(),
            },
            expect: {
                b: (() => {
                    const m = matrix_mxn(6, 6);
                    // do stuff here
                    return m.toArr();
                })(),
            },
        },
    },
    ctrsmErrors: {
        case0: {
            desc: 'a has no imaginary part',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                //
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
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8), // B = alpha * A * B
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
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8), // B = alpha * A * B
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
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case5: {
            desc: 'transA!="ntc"',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'x',
                diag: 'n',
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case6: {
            desc: 'diag!="un"',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'x',
                m: 4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case7: {
            desc: 'm<0',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: -4,
                n: 6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case8: {
            desc: 'n<0',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4,
                n: -6,
                lda: 6,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case9: {
            desc: 'lda< max(1,nrowA)',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4,
                n: 6,
                lda: 3,
                ldb: 4,
                alpha: complex(0.2, 0.8),
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            }
        },
        case10: {
            desc: 'ldb< max(1,m)',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'n',
                m: 4, // A = m*m marrix , lda >= m
                n: 6, // rows in B
                lda: 6, //NBupper 4x4 of A is referenced
                ldb: 3,
                alpha: complex(0.2, 0.8),  // B = alpha * A * B
                //dummies, its not a data test
                a: [complex(0, 0)],
                b: [complex(0, 0)]
            },
        },
    },
}
