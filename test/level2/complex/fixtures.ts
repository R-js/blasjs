import { complex, fortranArrComplex64 as arr64 } from '../../../src/lib/f_func';
import { bandmatrix_nxm_ku_kl, matrix_nxm, vector } from './matrices';

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
    // CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    cgemv: {
        case0: {
            desc: 'trans=n, inc(1,1), m=6, n=8, alpha(0.8,0.2), beta(0.2,0.8)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0.8, 0.2),
                beta: complex(0.2, 0.8),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1
            },
            expect: {
                y: [
                    complex(-4.0955083795571205, 1.1357517892058278),
                    complex(0.93896042218602016, 1.0838676654345458),
                    complex(-1.5853445176637875, 1.6962449347899724),
                    complex(-0.25555151836345036, -4.3800522220206233),
                    complex(-1.2051793139138998, -2.9157429210651786),
                    complex(0.38280021798680286, 0.19114175133433564),
                ]
            },
        },
        case1: {
            desc: 'trans=n, inc(-1,-1), m=6, n=8, alpha(0.8,0.2), beta(0,0)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0.8, 0.2),
                beta: complex(0, 0),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: -1,
                y: vector(6),
                incy: -1
            },
            expect: {
                y: [
                    complex(2.6017477475494872, -0.62451989822044895),
                    complex(1.8074832963400977, -2.2773562366746307),
                    complex(-4.2649596122103990, 2.1071421255650922),
                    complex(0.27753542357045247, -1.9840473842815003),
                    complex(2.2299686271484305, 0.24845760476943476),
                    complex(0.20069066517622025, -1.0213684617581498),
                ]
            },
        },
        case2: {
            desc: 'trans=n, inc(1,1), m=6, n=8, alpha(0.8,0.2), beta(1,0)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0.8, 0.2),
                beta: complex(1, 0),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1
            },
            expect: {
                y: [
                    complex(-3.9970026807046444, 2.2726735819034305),
                    complex(0.66801290992029672, 1.0035901828539755),
                    complex(-1.3938842013734591, 0.82508814026657729),
                    complex(0.29003966083788602, -5.5960115682700131),
                    complex(-0.29257262206414292, -2.2331706205478370),
                    complex(6.7974918648228133E-002, 6.4722209339979742E-002),
                ]
            },
        },
        case3: {
            desc: 'trans=t, inc(1,1), m=6, n=8, alpha(0.8,0.2), beta(1,0)',
            input: {
                trans: 't',
                m: 6,
                n: 8,
                alpha: complex(0.8, 0.2),
                beta: complex(1, 0),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: vector(8),
                incy: 1
            },
            expect: {
                y: [
                    complex(-1.0722266311585562, 1.4276441972981582),
                    complex(4.0697286775452355, -1.3751872740072097),
                    complex(-1.7667910641302869, -1.3146598668049647),
                    complex(2.7682588305027580, 0.43286241964465866),
                    complex(0.61287748489518667, -1.2227666045246042),
                    complex(-5.6604442768105556E-002, 1.3467393639005534),
                    complex(-2.2131934192050711, 4.1552655661274187),
                    complex(-2.9062411441513247, 2.3044710975088387),
                ],
            },
        },
        case4: {
            desc: 'trans=c, inc(1,1), m=6, n=8, alpha(0.8,0.2), beta(1,0)',
            input: {
                trans: 'c',
                m: 6,
                n: 8,
                alpha: complex(0.8, 0.2),
                beta: complex(1, 0),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: vector(8),
                incy: 1
            },
            expect: {
                y: [
                    complex(2.2502367460895774, 2.0732672240604177),
                    complex(1.5386279370583802, 1.2491393341728396),
                    complex(3.5354757160882668, 0.53153277666431009),
                    complex(-0.65968123771844311, -0.51832467085111289),
                    complex(-8.1751634328745854E-002, 0.52362043122181512),
                    complex(-1.8050991363430615, -0.94441463917334822),
                    complex(-0.25656634156597757, 1.6792353267609732),
                    complex(-2.5932578406614151, -3.5422230487806061),
                ]
            },
        },
        case5: {
            desc: 'trans=n, inc(1,1), m=6, n=8, alpha(0,0), beta(0.2,0.8)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                beta: complex(0.2, 0.8),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1
            },
            expect: {
                y: [
                    complex(-0.74751576107018813, -0.36477962083786153),
                    complex(0.15177874642280853, -0.13923813908156957),
                    complex(0.47267537821345584, 0.44634650418678667),
                    complex(0.55537789695527984, 0.79697925443916384),
                    complex(-0.76883520735554689, 0.31441456547997237),
                    complex(0.19707170422482045, -0.14935848339436175),
                ]
            },
        },
        case6: {
            desc: '(trivial: alpha=0, beta=1)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                beta: complex(1, 0),
                a: matrix_nxm(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1
            },
            expect: {
                y: vector(6).toArr()
            },
        },
    },
    cgemvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [0],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case1: {
            desc: 'y has no imaginary part',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            }
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            }
        },
        case3: {
            desc: 'trans !="ntc"',
            input: {
                trans: 'x',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case4: {
            desc: 'n<0',
            input: {
                trans: 'n',
                m: 6,
                n: -8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case5: {
            desc: 'm<0',
            input: {
                trans: 'n',
                m: -6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case6: {
            desc: 'lda<max(1,m)',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 4,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case7: {
            desc: 'incx=0',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 0,
                beta: complex(2.5, +0.5),
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
        case8: {
            desc: 'incy=0',
            input: {
                trans: 'n',
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                beta: complex(2.5, +0.5),
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            }
        },
    }
}

