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
            desc: 'y = alpha*A*x + beta*y, m=6,n=8,kl=1, ku=1, alpha(0,0), beta(2.5,0.5)',
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
        case1: {
            desc: 'y = alpha*A*x + beta*y, m=6,n=8,kl=1, ku=1, alpha(0.2,0.8), beta(2.5,0.5)',
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
                alpha: complex(0.2, 0.8),
                //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: bandmatrix_nxm_ku_kl(8, 6, 6, 3, 2),
                x: vector(8),
                y: vector(6),
            },
            expect: {
                y: [
                    complex(-1.0943454405994295, 1.0642741157616680),
                    complex(-0.38628598141225878, -1.2014790478399264),
                    complex(0.88995188622798471, -1.1623543999107093),
                    complex(7.9118205360687366E-002, 0.76470894087791352),
                    complex(-5.7862862678327798E-002,
                        0.62930340373117366),
                    complex(-0.36689070841052379, -0.18794836059332087)
                ]
            },
        },
        case2: {
            desc: 'trans="t",m=6,n=8,kl=1, ku=1, alpha(0.2,0.8), beta(2.5,0.5)',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                kl: 3,
                ku: 2,
                lda: 6,
                incx: 1,
                incy: 1,
                beta: complex(0, 0),
                alpha: complex(0.2, 0.8),
                //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: bandmatrix_nxm_ku_kl(8, 6, 6, 3, 2),
                x: vector(6),
                y: vector(8),
            },
            expect: {
                y: [
                    complex(-0.38810574366012285, -1.6303050664096645),
                    complex(2.3103061857743366, -0.97417280945264095),
                    complex(-4.2076520090095242, -0.50591381502069677),
                    complex(-0.81853765257644051, 1.4503261671516843),
                    complex(2.6631576242668459E-002, -0.63622673851522416),
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(0.0000000000000000, 0.0000000000000000),
                ]
            },
        },
        case3: {
            desc: 'trans="c",incx=-1, incy=-1,m=6,n=8,kl=1, ku=1, alpha(0.2,0.8), beta(2.5,0.5)',
            input: {
                trans: 'c',
                m: 6,
                n: 8,
                kl: 3,
                ku: 2,
                lda: 6,
                incx: -1,
                incy: -1,
                beta: complex(1, 0),
                alpha: complex(0.2, 0.8),
                //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: bandmatrix_nxm_ku_kl(8, 6, 6, 3, 2),
                x: vector(6),
                y: vector(8),
            },
            expect: {
                y: [
                    complex(-0.64901006221771240, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(0.76170726548233147, -0.74829467121041837),
                    complex(0.49362246616457550, -0.11820536025789097),
                    complex(-0.14226157594261568, -3.1888940811820188),
                    complex(-2.0779711022358427, 1.5114922978334300),
                    complex(-2.4265906255376652, 0.38932194147661164),
                ]
            },
        },
        case4: {
            desc: '(trivial, alpha=1, beta=0)',
            input: {
                trans: 'c',
                m: 6,
                n: 8,
                kl: 3,
                ku: 2,
                lda: 6,
                incx: -1,
                incy: -1,
                beta: complex(1, 0),
                alpha: complex(0, 0),
                //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: bandmatrix_nxm_ku_kl(8, 6, 6, 3, 2),
                x: vector(6),
                y: vector(8),
            },
            expect: {
                y: [
                    complex(-0.64901006221771240, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                    complex(-0.91206836700439453, 1.2560187578201294),
                    complex(-1.4375861883163452, 0.64667439460754395),
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
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
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case8: {
            desc: 'x has no imaginary',
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
                x: [0],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case9: {
            desc: 'y has no imaginary',
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
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            }
        },
        case10: {
            desc: 'y has no imaginary',
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
                x: [complex(0, 0)],
                a: [0],
                y: [complex(0, 0)],
            }
        },
    },

}
