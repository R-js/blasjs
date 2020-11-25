/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import { complex } from '../../test-helpers';

import { bandmatrix_nxm_ku_kl, matrix_mxn, vector } from '../../matrices';

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
                    complex(2.9619127362966537, -0.496965691447258),
                    complex(-0.139064721763134, 2.5643529072403908),
                    complex(-0.15649497509002686, -0.74832186102867126),
                ],
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
                    complex(-1.0943454405994295, 1.064274115761668),
                    complex(-0.38628598141225878, -1.2014790478399264),
                    complex(0.88995188622798471, -1.1623543999107093),
                    complex(7.9118205360687366e-2, 0.76470894087791352),
                    complex(-5.7862862678327798e-2, 0.62930340373117366),
                    complex(-0.36689070841052379, -0.18794836059332087),
                ],
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
                    complex(2.6631576242668459e-2, -0.63622673851522416),
                    complex(0, 0),
                    complex(0, 0),
                    complex(0, 0),
                ],
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
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(0.76170726548233147, -0.74829467121041837),
                    complex(0.4936224661645755, -0.11820536025789097),
                    complex(-0.14226157594261568, -3.1888940811820188),
                    complex(-2.0779711022358427, 1.51149229783343),
                    complex(-2.4265906255376652, 0.38932194147661164),
                ],
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
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                    complex(-0.91206836700439453, 1.2560187578201294),
                    complex(-1.4375861883163452, 0.64667439460754395),
                ],
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1,
            },
            expect: {
                y: [
                    complex(-4.0955083795571205, 1.1357517892058278),
                    complex(0.93896042218602016, 1.0838676654345458),
                    complex(-1.5853445176637875, 1.6962449347899724),
                    complex(-0.25555151836345036, -4.3800522220206233),
                    complex(-1.2051793139138998, -2.9157429210651786),
                    complex(0.38280021798680286, 0.19114175133433564),
                ],
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: -1,
                y: vector(6),
                incy: -1,
            },
            expect: {
                y: [
                    complex(2.6017477475494872, -0.62451989822044895),
                    complex(1.8074832963400977, -2.2773562366746307),
                    complex(-4.264959612210399, 2.1071421255650922),
                    complex(0.27753542357045247, -1.9840473842815003),
                    complex(2.2299686271484305, 0.24845760476943476),
                    complex(0.20069066517622025, -1.0213684617581498),
                ],
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1,
            },
            expect: {
                y: [
                    complex(-3.9970026807046444, 2.2726735819034305),
                    complex(0.66801290992029672, 1.0035901828539755),
                    complex(-1.3938842013734591, 0.82508814026657729),
                    complex(0.29003966083788602, -5.5960115682700131),
                    complex(-0.29257262206414292, -2.233170620547837),
                    complex(6.7974918648228133e-2, 6.4722209339979742e-2),
                ],
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: vector(8),
                incy: 1,
            },
            expect: {
                y: [
                    complex(-1.0722266311585562, 1.4276441972981582),
                    complex(4.0697286775452355, -1.3751872740072097),
                    complex(-1.7667910641302869, -1.3146598668049647),
                    complex(2.768258830502758, 0.43286241964465866),
                    complex(0.61287748489518667, -1.2227666045246042),
                    complex(-5.6604442768105556e-2, 1.3467393639005534),
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: vector(8),
                incy: 1,
            },
            expect: {
                y: [
                    complex(2.2502367460895774, 2.0732672240604177),
                    complex(1.5386279370583802, 1.2491393341728396),
                    complex(3.5354757160882668, 0.53153277666431009),
                    complex(-0.65968123771844311, -0.51832467085111289),
                    complex(-8.1751634328745854e-2, 0.52362043122181512),
                    complex(-1.8050991363430615, -0.94441463917334822),
                    complex(-0.25656634156597757, 1.6792353267609732),
                    complex(-2.5932578406614151, -3.5422230487806061),
                ],
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1,
            },
            expect: {
                y: [
                    complex(-0.74751576107018813, -0.36477962083786153),
                    complex(0.15177874642280853, -0.13923813908156957),
                    complex(0.47267537821345584, 0.44634650418678667),
                    complex(0.55537789695527984, 0.79697925443916384),
                    complex(-0.76883520735554689, 0.31441456547997237),
                    complex(0.19707170422482045, -0.14935848339436175),
                ],
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
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(8),
                incx: 1,
                y: vector(6),
                incy: 1,
            },
            expect: {
                y: vector(6).toArr(),
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
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
            },
        },
    },
    // CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    cgerc: {
        case0: {
            desc: 'inc(1,1), m=6, n=8, alpha(0.2,0.8)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: (() => {
                    const v = vector(8);
                    v.s(4)(0, 0);

                    return v;
                })(),
                incy: 1,
            },
            expect: {
                a: [
                    complex(1.4664377569938685, 1.8060944675544837),
                    complex(-0.53225092996643109, -0.45634091440603552),
                    complex(1.3676711769934335, 0.5836481726123065),
                    complex(1.5273647960307939, -1.225424538318252),
                    complex(1.1563959711297374, 2.1474949691654932),
                    complex(-1.7831775159348064, 0.50551385861110387),
                    complex(-0.759411695954748, -0.57340502179859332),
                    complex(-0.28224278033580208, -0.78213264706952679),
                    complex(-0.16007534447744504, -1.1160014767921838),
                    complex(2.1635202157291209, -1.0386515298872647),
                    complex(0.7861956692836497, -1.7700218301387962),
                    complex(-0.78970753495767643, 1.2175962397684426),
                    complex(-1.4891467992751752, 0.27223155051811232),
                    complex(-0.12951010544349723, -0.2553245327718181),
                    complex(-0.17490710560589115, 0.76336937172774166),
                    complex(-0.38122952832127499, 0.38852989041687136),
                    complex(-0.39195398537726411, 2.3235693548980874),
                    complex(-0.69758972915325101, -0.8108153197236978),
                    complex(0.43568331003189087, -5.4877474904060364e-2),
                    complex(-1.2375384569168091, 0.25014132261276245),
                    complex(-0.22426788508892059, 0.61824327707290649),
                    complex(0.37739565968513489, -0.17262350022792816),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-0.52825871591104467, 1.0515473739625048),
                    complex(0.38661084813518909, -0.18238536895940438),
                    complex(1.5987282580693853, -1.3477283962554565),
                    complex(0.183471503102556, -0.55494690365844279),
                    complex(-1.0816686910303539, -3.246252501847624e-3),
                    complex(-7.3848983873964338e-2, 2.4312089369662865e-2),
                    complex(-4.7085681056201645e-2, -1.5882927279464436),
                    complex(-0.52236195772798288, 0.42419417151368499),
                    complex(-0.61206211258427001, 0.32620723072294699),
                    complex(-0.93465874966505846, 0.12460202594829736),
                    complex(0.73057521316282659, -0.22989490997416154),
                    complex(1.1698956214925516, 0.3292739186846686),
                    complex(1.2157758131530922, 2.5345520058413213),
                    complex(-0.74283140469563169, -0.3898756075997003),
                    complex(1.3678114595993076, 0.32901162588415334),
                    complex(0.21513199437223562, -0.15203527634896519),
                    complex(2.8540439631030279, 1.0935453027260844),
                    complex(0.19340614948484969, -1.6512506222549119),
                    complex(0.38594072116787248, 7.9234395057932172e-2),
                    complex(-1.1402800589724982, -0.19270504442970846),
                    complex(-1.5574412841395437, -0.95309580634241586),
                    complex(-1.3486081371207643, 0.8997782769273901),
                    complex(-0.25519137193589581, 0.80878144081759595),
                    complex(0.77664318903912266, -0.71173479570372455),
                ],
            },
        },

        case1: {
            desc: 'inc(-1,-1), m=6, n=8, alpha(0.8,0.2)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: -1,
                y: vector(8),
                incy: -1,
            },
            expect: {
                a: [
                    complex(0.88306036814678868, 1.0794348476983751),
                    complex(0.98235736678213703, -0.38432509127926684),
                    complex(1.0467817773918699, -0.26657089566477299),
                    complex(0.88155860714158418, -1.2266749209248053),
                    complex(0.10640467108044993, 1.8599184957333317),
                    complex(-0.70122531996791548, 1.5685475206682837),
                    complex(-1.2959069761255018, -0.56408458111912407),
                    complex(0.80142042293998772, -0.15313942732686359),
                    complex(0.48871110895462522, -2.5910351676300931),
                    complex(2.5341606321914707, -2.0663782886284809),
                    complex(0.4502751274012311, -1.6274243463177331),
                    complex(-0.57539383024900737, 2.4281348267336553),
                    complex(-1.1296731508234275, 0.9039826486670447),
                    complex(-0.28563710120446772, -0.47637998868149889),
                    complex(-0.58440221723440966, 0.32545120885206402),
                    complex(-0.59026260782734252, -0.29890814896817697),
                    complex(0.27274978786794912, 2.499617561186918),
                    complex(-0.70330023119372154, -0.95853341627041833),
                    complex(0.3151081552306404, -0.2728288665553249),
                    complex(-1.0346078437479442, 1.0618637752882218),
                    complex(0.65015746875979019, 0.17912169866508931),
                    complex(0.89035450285278694, -0.57970274468132166),
                    complex(1.6339225000576024e-2, -2.3952402088207974),
                    complex(0.33303756158219378, -0.57079592781056654),
                    complex(0.22244130892187952, 0.27685893490110408),
                    complex(-0.47458924311264428, 1.2988589637548592e-2),
                    complex(1.3633048638346255, 0.16949264392950569),
                    complex(-0.3575631755758607, 0.57362994646994336),
                    complex(-1.0591575896740089, -0.90467331741663981),
                    complex(-0.62343016692152187, -0.47254182333647243),
                    complex(-4.1375179015731151e-2, -1.4405746313997232),
                    complex(-1.1870657309731962, 0.24814596522485438),
                    complex(-0.40302903307820248, 1.0136452701079954),
                    complex(-0.52516363803653998, 0.56252018882397503),
                    complex(0.88670220892379703, -8.8394540644807507e-3),
                    complex(0.81042197304080377, -0.30247717946426378),
                    complex(1.0014621084444231, 1.3240134188761086),
                    complex(-0.40691086281321309, -0.53247309142076349),
                    complex(0.99717104313695781, 1.3567383846253696),
                    complex(-0.43365445905983463, 1.3229984144889442),
                    complex(1.7703807598272381, 0.46455208298342121),
                    complex(0.72990142965560356, -1.6605710629343813),
                    complex(-0.69601147479901848, -0.98379926699924769),
                    complex(-9.0288758923210644e-2, 9.4871429002452989e-2),
                    complex(-0.91163509525033404, -0.95184542373586234),
                    complex(-1.0277187375192007, 1.7499973452044695),
                    complex(-1.7697996686844639, 0.73676561769082727),
                    complex(1.3600205778862025, 1.4924824152384053e-2),
                ],
            },
        },
        case2: {
            desc: '(trivial case)alpha=(0,0)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: -1,
                y: vector(8),
                incy: -1,
            },
            expect: {
                a: [
                    complex(1.2629542350769043, 0.9921603798866272),
                    complex(-0.32623335719108582, -0.42951309680938721),
                    complex(1.3297992944717407, 1.2383041381835938),
                    complex(1.272429347038269, -0.2793462872505188),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-0.92856705188751221, -0.4527839720249176),
                    complex(-0.29472044110298157, -0.83204329013824463),
                    complex(-5.7671726681292057e-3, -1.1665705442428589),
                    complex(2.4046533107757568, -1.0655906200408936),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-1.147657036781311, 0.83204710483551025),
                    complex(-0.28946158289909363, -0.22732868790626526),
                    complex(-0.29921510815620422, 0.26613736152648926),
                    complex(-0.41151082515716553, -0.37670272588729858),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.43568331003189087, -5.4877474904060364e-2),
                    complex(-1.2375384569168091, 0.25014132261276245),
                    complex(-0.22426788508892059, 0.61824327707290649),
                    complex(0.37739565968513489, -0.17262350022792816),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-5.7106774300336838e-2, 0.35872888565063477),
                    complex(0.50360798835754395, -1.1045478284358978e-2),
                    complex(1.0857694149017334, -0.94064915180206299),
                    complex(-0.69095385074615479, -0.11582532525062561),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-0.23570655286312103, -1.4250984191894531),
                    complex(-0.54288828372955322, 0.36594113707542419),
                    complex(-0.43331032991409302, 0.24841265380382538),
                    complex(-0.64947164058685303, 6.5288178622722626e-2),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                    complex(0.9921603798866272, 1.2629542350769043),
                    complex(-0.42951309680938721, -0.32623335719108582),
                    complex(1.2383041381835938, 1.3297992944717407),
                    complex(-0.2793462872505188, 1.272429347038269),
                    complex(1.7579030990600586, 0.41464143991470337),
                    complex(0.56074607372283936, -1.5399500131607056),
                    complex(-0.4527839720249176, -0.92856705188751221),
                    complex(-0.83204329013824463, -0.29472044110298157),
                    complex(-1.1665705442428589, -5.7671726681292057e-3),
                    complex(-1.0655906200408936, 2.4046533107757568),
                    complex(-1.5637820959091187, 0.76359343528747559),
                    complex(1.1565370559692383, -0.79900926351547241),
                ],
            },
        },
    },
    cgercErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [0],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'y has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                beta: complex(2.5, +0.5),
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case3: {
            desc: 'n<0',
            input: {
                m: 6,
                n: -8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'm<0',
            input: {
                m: -6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'lda<max(1,m)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 2,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 0,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'incy=0',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    cgeru: {
        case0: {
            desc: 'inc(1,1), m=6, n=8, alpha(0.2,0.8)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: 1,
                y: (() => {
                    const v = vector(8);
                    v.s(4)(0, 0);

                    return v;
                })(),
                incy: 1,
            },
            expect: {
                a: [
                    complex(2.0297612143617068, 0.6517175810502388),
                    complex(-0.31722765173412737, -0.22195137259592293),
                    complex(0.67838525870385646, 1.3135933588490314),
                    complex(0.29660221113116736, -0.36776314720236125),
                    complex(0.67085048022165195, 0.96019479554600029),
                    complex(-1.5525255484272074, 0.80984860603561071),
                    complex(-0.91956134643055376, -0.24522224781145335),
                    complex(-0.34337267365494262, -0.84876825882173301),
                    complex(3.588471620912579e-2, -1.3235207357780114),
                    complex(2.5134190087092048, -1.2824797785023652),
                    complex(0.92423348688558515, -1.43247915314201),
                    complex(-0.85528057562333193, 1.131075604438585),
                    complex(-1.7990710725491952, 0.90733632550094778),
                    complex(-0.24780969402183861, -0.38427887944141781),
                    complex(0.20431807046274741, 0.36177464242009244),
                    complex(0.29590044861983555, -8.3330600887344008e-2),
                    complex(-0.12482089968205246, 2.9767875702135211),
                    complex(-0.82448777054323985, -0.97825149550145041),
                    complex(0.43568331003189087, -5.4877474904060364e-2),
                    complex(-1.2375384569168091, 0.25014132261276245),
                    complex(-0.22426788508892059, 0.61824327707290649),
                    complex(0.37739565968513489, -0.17262350022792816),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(0.19910226600661174, -0.43897941786342354),
                    complex(0.66424803995565351, 0.12025746448274954),
                    complex(0.70872505335330216, -0.4052261083371993),
                    complex(-1.4056841951938006, 0.55246203420080509),
                    complex(-1.7086030755539241, -1.5362834602014517),
                    complex(0.22396790866494898, 0.4172678909133698),
                    complex(-0.2482820881296228, -1.1759958868766818),
                    complex(-0.59915959583741274, 0.34047968554477093),
                    complex(-0.36587699745669194, 6.5500265815802733e-2),
                    complex(-0.49508001953504838, -0.18172001358543421),
                    complex(0.90399246920586351, 0.19416080069351449),
                    complex(1.0875160462414502, 0.22057782778244783),
                    complex(2.1321159056388295, 0.65676437050062852),
                    complex(-0.39305997571479284, -8.6017025087555266e-3),
                    complex(0.24657229610721754, 1.5163899086758028),
                    complex(-1.7869097919659473, 1.2430948361600873),
                    complex(2.0642227791536065, -0.83779758149610695),
                    complex(0.56860026305061906, -1.1561991079710021),
                    complex(0.85772800210885958, -0.88756420944139058),
                    complex(-0.96019658037878808, 3.5978134848148047e-3),
                    complex(-2.1347229948999118, -0.34176167823825904),
                    complex(-2.3793802908792001, 1.6180756075113232),
                    complex(-0.66183906951100635, -0.1855906437216317),
                    complex(0.96981580269622225, -0.45685234565599903),
                ],
            },
        },

        case1: {
            desc: 'inc(-1,-1), m=6, n=8, alpha(0.2,0.8)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: -1,
                y: vector(8),
                incy: -1,
            },
            expect: {
                a: [
                    complex(1.0762329818038883, 1.3343172977461006),
                    complex(0.57570966920702649, -1.3786971758184945),
                    complex(1.6009623633433945e-2, 0.45172643491916009),
                    complex(0.30427689638121636, -0.61534079282064869),
                    complex(0.28648814967415992, 2.0562213536478549),
                    complex(-0.22943803902692839, 0.60174891616896098),
                    complex(-0.9207128625597325, -6.9033066835214152e-2),
                    complex(1.1599238990566474e-2, -2.0844823115490549),
                    complex(-1.5133306773835578, -1.1959050551210406),
                    complex(1.4129214686993805, -0.87900000583683158),
                    complex(0.8000465563820699, -1.2461504412267883),
                    complex(0.34094626223672986, 0.55034719139296251),
                    complex(-1.2120527260745289, 0.79528655776482393),
                    complex(-0.1122198451614308, -5.2324278013822889e-2),
                    complex(-0.14482348710439957, 1.9129169318332417e-2),
                    complex(-0.34407749269976445, -0.55961511387532126),
                    complex(0.19595214975851932, 2.415903075218004),
                    complex(-0.90449663826714266, -0.54623657520065649),
                    complex(0.6129250477695537, 0.12012693498838201),
                    complex(-1.6615422282715144, -0.47117343241138221),
                    complex(-0.93899822953656631, 1.2865306365243372),
                    complex(3.5129813670364829e-4, 0.36279954323693553),
                    complex(0.29397641682104048, -2.0925973753786433),
                    complex(1.0603985434998502, -2.0613227196364949),
                    complex(9.7284846751467785e-2, 0.11172069344247793),
                    complex(-0.21112235609010166, 0.65724188116707161),
                    complex(2.031141746226159, -0.29589192058188163),
                    complex(1.6457423030846297e-2, 0.17754679974932891),
                    complex(-1.175833606265771, -1.0318578636387787),
                    complex(-0.92910096497981542, 0.15384662106908503),
                    complex(-0.16827322040571996, -1.6080108071774757),
                    complex(-0.91993264527798446, 0.901364180540288),
                    complex(0.27410094386290806, 0.54178477880377995),
                    complex(-0.14593846196790139, 0.16092545951632581),
                    complex(0.7684026203454557, -0.1377938007340804),
                    complex(0.50049769976678371, 0.33262759551857168),
                    complex(0.93588906777876768, 1.237492783546251),
                    complex(-0.26887304521127764, -0.19493041442397727),
                    complex(1.3470698361170415, 1.1129101360102691),
                    complex(-0.23769439837326381, 1.1154791555031165),
                    complex(1.7092508665080977, 0.39791647123121499),
                    complex(0.5697517791797978, -1.3323882889472414),
                    complex(-0.46535950729141939, -0.67946451957474086),
                    complex(-0.57583424983129605, -1.09242874461704),
                    complex(-2.1423976801499607, -9.4184032619971658e-2),
                    complex(-1.7170046558087777, 2.4799425314411945),
                    complex(-1.5547763904521601, 0.97115515950093989),
                    complex(1.923344035254041, -1.1394520623518609),
                ],
            },
        },
        case2: {
            desc: '(trivial case)alpha=(0,0)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                a: matrix_mxn(6, 8, 6),
                lda: 6,
                x: vector(6),
                incx: -1,
                y: vector(8),
                incy: -1,
            },
            expect: {
                a: [
                    complex(1.2629542350769043, 0.9921603798866272),
                    complex(-0.32623335719108582, -0.42951309680938721),
                    complex(1.3297992944717407, 1.2383041381835938),
                    complex(1.272429347038269, -0.2793462872505188),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-0.92856705188751221, -0.4527839720249176),
                    complex(-0.29472044110298157, -0.83204329013824463),
                    complex(-5.7671726681292057e-3, -1.1665705442428589),
                    complex(2.4046533107757568, -1.0655906200408936),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-1.147657036781311, 0.83204710483551025),
                    complex(-0.28946158289909363, -0.22732868790626526),
                    complex(-0.29921510815620422, 0.26613736152648926),
                    complex(-0.41151082515716553, -0.37670272588729858),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.43568331003189087, -5.4877474904060364e-2),
                    complex(-1.2375384569168091, 0.25014132261276245),
                    complex(-0.22426788508892059, 0.61824327707290649),
                    complex(0.37739565968513489, -0.17262350022792816),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-5.7106774300336838e-2, 0.35872888565063477),
                    complex(0.50360798835754395, -1.1045478284358978e-2),
                    complex(1.0857694149017334, -0.94064915180206299),
                    complex(-0.69095385074615479, -0.11582532525062561),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-0.23570655286312103, -1.4250984191894531),
                    complex(-0.54288828372955322, 0.36594113707542419),
                    complex(-0.43331032991409302, 0.24841265380382538),
                    complex(-0.64947164058685303, 6.5288178622722626e-2),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                    complex(0.9921603798866272, 1.2629542350769043),
                    complex(-0.42951309680938721, -0.32623335719108582),
                    complex(1.2383041381835938, 1.3297992944717407),
                    complex(-0.2793462872505188, 1.272429347038269),
                    complex(1.7579030990600586, 0.41464143991470337),
                    complex(0.56074607372283936, -1.5399500131607056),
                    complex(-0.4527839720249176, -0.92856705188751221),
                    complex(-0.83204329013824463, -0.29472044110298157),
                    complex(-1.1665705442428589, -5.7671726681292057e-3),
                    complex(-1.0655906200408936, 2.4046533107757568),
                    complex(-1.5637820959091187, 0.76359343528747559),
                    complex(1.1565370559692383, -0.79900926351547241),
                ],
            },
        },
    },
    cgeruErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [0],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'y has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                beta: complex(2.5, +0.5),
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case3: {
            desc: 'n<0',
            input: {
                m: 6,
                n: -8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'm<0',
            input: {
                m: -6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'lda<max(1,m)',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 2,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 0,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'incy=0',
            input: {
                m: 6,
                n: 8,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    chbmv: {
        case0: {
            desc: 'uplo=u, inc(1,1), n=6, k=3, alpha(0.2,0.8), beta(0.3,-0.7)',
            input: {
                uplo: 'u',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(0.3, -0.7),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).upperBand(3),
                lda: 6,
                incx: 1,
                incy: 1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.2534610371023731, 9.4703331729354634e-2),
                    complex(-0.57417853027717758, -6.9983511992388614e-3),
                    complex(-0.24071061122272883, -0.24386226459928881),
                    complex(-1.9639075687356289e-3, 0.89891271445795673),
                    complex(1.5460009072147165, 0.25069302006213029),
                    complex(-1.0003110232032266, -0.49852299275348511),
                ],
            },
        },
        case1: {
            desc: 'uplo=u, inc(-1,-1), n=6, k=3, alpha(0.2,0.8), beta(0,0)',
            input: {
                uplo: 'u',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(0, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).upperBand(3),
                lda: 6,
                incx: -1,
                incy: -1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.2736547002392706, -1.2409869471212633),
                    complex(-1.5069594023144455, 0.16300769978973884),
                    complex(1.3847341260375996, 0.25254847949457943),
                    complex(0.41526650594237807, -0.14984859659953415),
                    complex(-1.1204385690487695, -0.34592548946596602),
                    complex(0.86249408001175165, -0.56909457127832508),
                ],
            },
        },
        case2: {
            desc: 'uplo=u, inc(-1,-1), n=6, k=3, alpha(0,0), beta(0,0)',
            input: {
                uplo: 'u',
                n: 6,
                k: 3,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).upperBand(3),
                lda: 6,
                incx: -1,
                incy: -1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0)],
            },
        },
        case3: {
            desc: 'trival, beta=(1,0), alpha=(0,0)',
            input: {
                uplo: 'u',
                n: 6,
                k: 3,
                alpha: complex(0, 0),
                beta: complex(1, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).upperBand(3),
                lda: 6,
                incx: -1,
                incy: -1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })().toArr(),
            },
        },
        case4: {
            desc: 'upl="u",inc(-1,-1), beta=(1,0), alpha=(0.2,0.8)',
            input: {
                uplo: 'u',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).upperBand(3),
                lda: 6,
                incx: -1,
                incy: -1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.922664762456983, -0.46884477526152213),
                    complex(-1.5069594023144455, 0.16300769978973884),
                    complex(2.0488698205413836, -0.17226181084202949),
                    complex(0.41526650594237807, -0.14984859659953415),
                    complex(-0.97666708455456019, 0.65106137653134821),
                    complex(0.74474048489799738, -0.84487259666704273),
                ],
            },
        },
        case5: {
            desc: 'uplo="l", beta=(1,0), alpha=(0.2,0.8)',
            input: {
                uplo: 'l',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).lowerBand(3),
                lda: 6,
                incx: 1,
                incy: 1,
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(0.79654686947656517, 0.93230872508886198),
                    complex(3.8235731030044517e-2, 3.6497121577365221),
                    complex(0.38312451728225538, -3.0675595198095778),
                    complex(-1.6321936347511552, -2.1997081054941319),
                    complex(3.1857501193367059, 7.729882066971816e-2),
                    complex(0.28503740823625501, 1.0786387593676938),
                ],
            },
        },
    },
    chbmvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                y: [complex(0, 0)],
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'y has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case3: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                k: 3,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'lda<(k+1)',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                k: 4,
                lda: 4,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'k<0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                k: -4,
                lda: 4,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                k: 4,
                lda: 5,
                incx: 0,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'incy=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                k: 4,
                lda: 5,
                incx: 1,
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: complex(0, 0),
                k: 4,
                lda: 5,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    chemv: {
        case0: {
            desc: 'uplo=u, inc(1,1), n=6, k=3, alpha(0.2,0.8), beta(0.3,-0.7)',
            input: {
                uplo: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                incy: 1,
                alpha: complex(0.2, 0.8),
                beta: complex(0.3, -0.7),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.1978771580109713, 0.19256169948001889),
                    complex(-0.6440639123273405, 0.12304172454730536),
                    complex(-0.24071061122272883, -0.24386226459928881),
                    complex(-1.9639075687356289e-3, 0.89891271445795673),
                    complex(1.5927182864486373, 0.26326773843256362),
                    complex(-1.0745919479739299, -0.10110828501746819),
                ],
            },
        },
        case1: {
            desc: 'uplo=u, inc(-1,-1), n=6, alpha(-0.12,0.88), beta(-0.43,0.57)',
            input: {
                uplo: 'u',
                n: 6,
                lda: 6,
                incx: -1,
                incy: -1,
                alpha: complex(-0.12, 0.88),
                beta: complex(-0.43, 0.57),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-0.53727175687523987, -2.3112656523209361),
                    complex(-1.5860492969147173, -0.43291167355929333),
                    complex(1.2402782186527825, 1.367667879363657),
                    complex(0.47520593481607637, 1.6258004364783057e-2),
                    complex(-1.1206760713901169, -0.92993994108923006),
                    complex(1.669952422368981, -0.23413210623927322),
                ],
            },
        },
        case2: {
            desc: 'uplo=l, inc(-1,-1), n=6, alpha(-0.12,0.88), beta(-0.43,0.57)',
            input: {
                uplo: 'l',
                n: 6,
                lda: 6,
                incx: -1,
                incy: -1,
                alpha: complex(-0.12, 0.88),
                beta: complex(-0.43, 0.57),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-0.82109209135990002, -2.8288933164643106),
                    complex(1.8186627440037402, -1.9201765078819051),
                    complex(-4.8025707150157597, 2.6701570029481525),
                    complex(2.3815496358250265, 0.2074289110744601),
                    complex(-1.7405950102646242, 2.9574110901705071),
                    complex(3.6227184273789095, 3.3816196677830903),
                ],
            },
        },
        case3: {
            desc: 'trivial,  alpha(0,0), beta(0,0)',
            input: {
                uplo: 'l',
                n: 6,
                lda: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0), complex(0, 0)],
            },
        },
        case4: {
            desc: 'trivial,alpha=(0,0), beta=(1,0)',
            input: {
                uplo: 'l',
                n: 6,
                lda: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                beta: complex(1, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v.toArr();
                })(),
            },
        },
        case5: {
            desc: 'uplo="l", inc(-1,-1), alpha=(0,1), beta=(1,0)',
            input: {
                uplo: 'l',
                n: 6,
                lda: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 1),
                beta: complex(1, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: vector(6),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.7089394711323604, -1.5002952495366637),
                    complex(1.7368180095711416, -2.4188575862567556),
                    complex(-4.3244119692369942, 2.6519589765813314),
                    complex(2.688457354507543, -0.13089314190862256),
                    complex(-0.59245062296078643, 4.8521145635536458),
                    complex(4.1985782666268578, 2.9198979594665389),
                ],
            },
        },
    },
    chemvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                y: [complex(0, 0)],
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'y has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [0],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case3: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                alpha: complex(0.2, 0.8),
                beta: complex(1, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'lda<n',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                lda: 5,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                lda: 6,
                incx: 0,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'incy=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: complex(0, 0),
                lda: 6,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    //ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
    cher: {
        case0: {
            desc: 'uplo=u, incx=1, n=6, alpha=0.2, lda=6, x={})',
            input: {
                uplo: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                alpha: 0.2, // NOTE: MUST BE REAL
                a: (() => {
                    const m = matrix_mxn(6, 8);
                    //console.log(m.r);
                    const m1 = m.slice_used_for_test(1, 6, 1, 6);
                    //console.log(m1.r);
                    const m2 = m1.setLower(0);
                    //console.log(m2.i.length);
                    //process.exit(1);
                    return m2;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.4664377569938685, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.92856705188751221, -0.83204329013824463),
                    complex(-0.29472044110298157, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-1.2994659767673966, -0.17990848227932066),
                    complex(-0.28946158289909363, 0.26613736152648926),
                    complex(-0.17490710560589112, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.43568331003189087, 0.25014132261276245),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.37739565968513489, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(7.8194520501826237e-2, 0.14056783081885876),
                    complex(0.50360798835754395, -0.94064915180206299),
                    complex(1.0201601128637381, -0.26046736155203748),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-1.0816686910303539, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.26300986834671181, 0.31196009025279664),
                    complex(-0.54288828372955322, 0.24841265380382538),
                    complex(-0.42552053432548659, 0.11192357318741848),
                    complex(-0.64947164058685303, 1.9156390801072121e-2),
                    complex(0.6683753949148884, 0.24178842029114961),
                    complex(1.1698956214925516, 0.0),
                ],
            },
        },
        case1: {
            desc: 'uplo=u, incx=-1, n=6, alpha=0.2, lda=6, x={})',
            input: {
                uplo: 'u',
                n: 6,
                lda: 6,
                incx: -1,
                alpha: 0.2, // NOTE: MUST BE REAL
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.2809381210347879, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.98694238844082449, -0.81649333557626003),
                    complex(-9.1789827934116691e-2, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-1.147657036781311, -0.22732868790626526),
                    complex(-0.28946158289909363, 0.26613736152648926),
                    complex(-0.29921510815620422, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.4434731056204973, 0.20350592804806661),
                    complex(-1.3031477589548044, 0.76288531337431831),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.50170366223544804, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-5.7106774300336838e-2, -1.1045478284358978e-2),
                    complex(0.50360798835754395, -0.94064915180206299),
                    complex(1.0857694149017334, -0.11582532525062561),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-1.2845993041992188, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.26300986834671181, 0.41992218389805175),
                    complex(-0.40758698892739015, 9.6799344700607637e-2),
                    complex(-0.43331032991409302, 6.5288178622722626e-2),
                    complex(-0.80128058057293872, -2.8263814825872494e-2),
                    complex(0.72675073146820068, 0.25733837485313416),
                    complex(1.3553952574516321, 0.0),
                ],
            },
        },
        case2: {
            desc: '(trivial) n=0',
            input: {
                uplo: 'u',
                n: 6,
                lda: 6,
                incx: -1,
                alpha: 0, // NOTE: MUST BE REAL
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.2629542350769043, -0.42951309680938721),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.92856705188751221, -0.83204329013824463),
                    complex(-0.29472044110298157, -1.1665705442428589),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-1.147657036781311, -0.22732868790626526),
                    complex(-0.28946158289909363, 0.26613736152648926),
                    complex(-0.29921510815620422, -0.37670272588729858),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.43568331003189087, 0.25014132261276245),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.37739565968513489, -2.223900318145752),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-5.7106774300336838e-2, -1.1045478284358978e-2),
                    complex(0.50360798835754395, -0.94064915180206299),
                    complex(1.0857694149017334, -0.11582532525062561),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-1.2845993041992188, 0.24226348102092743),
                    complex(0.0, 0.0),
                    complex(-0.23570655286312103, 0.36594113707542419),
                    complex(-0.54288828372955322, 0.24841265380382538),
                    complex(-0.43331032991409302, 6.5288178622722626e-2),
                    complex(-0.64947164058685303, 1.9156390801072121e-2),
                    complex(0.72675073146820068, 0.25733837485313416),
                    complex(1.151911735534668, 0.0),
                ],
            },
        },
        case3: {
            desc: 'uplo=l, incx=-1, n=6, alpha=1.234, lda=6, x={})',
            input: {
                uplo: 'l',
                n: 6,
                lda: 6,
                incx: -1,
                alpha: 1.234, // NOTE: MUST BE REAL
                a: (() => {
                    const m = matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6);
                    const m2 = m.setUpper(0);
                    // console.log(m2.r); debug
                    // process.exit(1);
                    return m2;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.3739148068679738, 0.0),
                    complex(-0.68640916889390469, 1.1423609224868452),
                    complex(1.3297992944717407, -0.2793462872505188),
                    complex(1.320492383840858, 2.0456434716758203),
                    complex(0.41464143991470337, 0.56074607372283936),
                    complex(-1.7084114627576494, -0.78584701720584593),
                    complex(0.0, 0.0),
                    complex(0.95736139059139025, 0.0),
                    complex(-5.7671726681292057e-3, -1.0655906200408936),
                    complex(1.9998439338703389, -2.4562234231403806),
                    complex(0.76359343528747559, 1.1565370559692383),
                    complex(3.5799691038578052e-2, 1.7675011834827594),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.29921510815620422, 0.0),
                    complex(-0.41151082515716553, 2.4413645267486572),
                    complex(0.25222346186637878, -0.79533910751342773),
                    complex(-0.89192110300064087, -5.4877474904060364e-2),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(1.1443760038382795, 0.0),
                    complex(0.13333636522293091, -1.2636144161224365),
                    complex(-0.13247161795193985, 0.65131154232107824),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-1.2845993041992188, 0.0),
                    complex(4.6726170927286148e-2, -1.4250984191894531),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(2.4074050140643375, 0.0),
                ],
            },
        },
    },
    cherErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: 0.2,
                lda: 6,
                incx: 1,
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'A has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: 1,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case3: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                alpha: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'lda<n',
            input: {
                uplo: 'l',
                n: 6,
                alpha: 1.2,
                lda: 5,
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: 0.2,
                lda: 6,
                incx: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: 0.1,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    //CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    chpr2: {
        case0: {
            desc: 'uplo=u, incxy(1,1), n=6, alpha=(0.2,0.8), x={6},y={6})',
            input: {
                uplo: 'u',
                n: 6,
                incx: 1,
                incy: 1,
                alpha: complex(0.2, 0.8),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: [
                    complex(1.6699212789108326, 0.0),
                    complex(-0.759411695954748, -0.95266433991192034),
                    complex(-0.29472044110298157, 0.0),
                    complex(-1.4512749167534824, -0.13248827665237595),
                    complex(-0.44376975470840946, 0.21556829407581407),
                    complex(-5.0599103055578076e-2, 0.0),
                    complex(0.43568331003189087, 0.25014132261276245),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.37739565968513489, 0.0),
                    complex(0.21349581530398931, 0.2921811399220765),
                    complex(0.52621022235371806, -0.73440941757238531),
                    complex(0.95455081082574245, -0.40510939785344929),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-0.87873807786148905, 0.0),
                    complex(-0.29031318383030258, 0.25797904343016909),
                    complex(-0.53358655517175724, 0.18735347000462102),
                    complex(-0.41773073873688016, 0.15855896775211434),
                    complex(-0.64947164058685303, 1.9156390801072121e-2),
                    complex(0.61000005836157611, 0.22623846572916503),
                    complex(1.1878795074504351, 0.0),
                ],
            },
        },
        case1: {
            desc: 'uplo=l, incxy(-1,-1), n=6, alpha=(0.2,0.8), x={6},y={6})',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0.2, 0.8),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).packedLower(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: [
                    complex(1.2989220069926715, 0.0),
                    complex(-0.44298403029771038, 1.2072042290596245),
                    complex(1.3297992944717407, -0.2793462872505188),
                    complex(1.2880089382154818, 1.8511738881894504),
                    complex(0.42394316847249935, 0.49968688992363497),
                    complex(-1.5945566441278871, -0.56074606567017271),
                    complex(0.11114078523474813, 0.0),
                    complex(-5.7671726681292057e-3, -1.0655906200408936),
                    complex(2.2734347066997658, -1.8530661685119423),
                    complex(0.7861956692836497, 1.3627767901989158),
                    complex(-0.52840667391114637, 1.1352737230419456),
                    complex(-0.29921510815620422, 0.0),
                    complex(-0.41151082515716553, 2.4413645267486572),
                    complex(0.25222346186637878, -0.79533910751342773),
                    complex(-0.89192110300064087, -5.4877474904060364e-2),
                    complex(0.62601166478576098, 0.0),
                    complex(-2.0971806586384922e-2, -1.3141834835731117),
                    complex(0.50057162322073023, 0.45356929690452408),
                    complex(-1.2845993041992188, 0.0),
                    complex(0.21588152686005035, -1.5457194689631288),
                    complex(1.5588787793685963, 0.0),
                ],
            },
        },
        case2: {
            desc: 'trivial(alpha=0) uplo=l',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).packedLower(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).packedLower().toArr(),
            },
        },
    },
    chpr2Errors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                x: [0],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'y has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(1, 0),
                incy: 1,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
                y: [0],
            },
        },
        case3: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [0],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: complex(0, 0),
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                alpha: complex(0, 0),
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                incx: 0,
                incy: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'incy=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                incx: 1,
                incy: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },

        case8: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: 0.1,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
    },

    cher2: {
        case0: {
            desc: 'uplo=u, incxy(1,1), n=6, alpha=(0.2,0.8), x={6},y={6})',
            input: {
                uplo: 'u',
                n: 6,
                incx: 1,
                incy: 1,
                lda: 6,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.6699212789108326, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),

                    complex(-0.759411695954748, -0.95266433991192034),
                    complex(-0.29472044110298157, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),

                    complex(-1.4512749167534824, -0.13248827665237595),
                    complex(-0.44376975470840946, 0.21556829407581407),
                    complex(-5.0599103055578076e-2, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),

                    complex(0.43568331003189087, 0.25014132261276245),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.37739565968513489, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),

                    complex(0.21349581530398931, 0.2921811399220765),
                    complex(0.52621022235371806, -0.73440941757238531),
                    complex(0.95455081082574245, -0.40510939785344929),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-0.87873807786148905, 0.0),
                    complex(0.0, 0.0),

                    complex(-0.29031318383030258, 0.25797904343016909),
                    complex(-0.53358655517175724, 0.18735347000462102),
                    complex(-0.41773073873688016, 0.15855896775211434),
                    complex(-0.64947164058685303, 1.9156390801072121e-2),
                    complex(0.61000005836157611, 0.22623846572916503),
                    complex(1.1878795074504351, 0.0),
                ],
            },
        },
        case1: {
            desc: 'uplo=l, incxy(1,1), n=6, alpha=(0.2,0.8), x={6},y={6})',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0.2, 0.8),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: [
                    complex(1.2989220069926715, 0.0),
                    complex(-0.44298403029771038, 1.2072042290596245),
                    complex(1.3297992944717407, -0.2793462872505188),
                    complex(1.2880089382154818, 1.8511738881894504),
                    complex(0.42394316847249935, 0.49968688992363497),
                    complex(-1.5945566441278871, -0.56074606567017271),
                    complex(0.0, 0.0),
                    complex(0.11114078523474813, 0.0),
                    complex(-5.7671726681292057e-3, -1.0655906200408936),
                    complex(2.2734347066997658, -1.8530661685119423),
                    complex(0.7861956692836497, 1.3627767901989158),
                    complex(-0.52840667391114637, 1.1352737230419456),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-0.29921510815620422, 0.0),
                    complex(-0.41151082515716553, 2.4413645267486572),
                    complex(0.25222346186637878, -0.79533910751342773),
                    complex(-0.89192110300064087, -5.4877474904060364e-2),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.62601166478576098, 0.0),
                    complex(-2.0971806586384922e-2, -1.3141834835731117),
                    complex(0.50057162322073023, 0.45356929690452408),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(-1.2845993041992188, 0.0),
                    complex(0.21588152686005035, -1.5457194689631288),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(0.0, 0.0),
                    complex(1.5588787793685963, 0.0),
                ],
            },
        },
        case2: {
            desc: '(trivial, alpha=(0,0)',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).toArr(),
            },
        },
    },
    cher2Errors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [0],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'y has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                a: [complex(0, 0)],
                y: [0],
            },
        },
        case3: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [0],
            },
        },
        case4a: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                incx: 0,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },

        case7: {
            desc: 'incy=0',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: 0,
                lda: 6,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'lda < max(1, n)',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 4,
                alpha: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    //CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
    chpmv: {
        case0: {
            desc: 'uplo=u, incxy(1,1), n=6, beta(0.3,-0.7), alpha(0.2,-0.8))',
            input: {
                uplo: 'u',
                n: 6,
                incx: 1,
                incy: 1,
                alpha: complex(0.2, -0.8),
                beta: complex(0.3, -0.7),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower().packedUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(1.0565062967743115, 2.1785481727394473),
                    complex(-0.3606433625471539, -0.89535796823641833),
                    complex(0.26097073349390854, -1.1070832376210971),
                    complex(0.18952760417765735, -0.48786643064345248),
                    complex(-0.69670093526894639, -0.67058224989704085),
                    complex(0.2936001246278136, -0.16342854670718146),
                ],
            },
        },
        case1: {
            desc: 'uplo=u, incxy(-1,-1), n=6, beta(0.3,-0.7), alpha(0.2,-0.8))',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0.2, -0.8),
                beta: complex(0.3, -0.7),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(1.079278998718928, 1.7817909884398611),
                    complex(1.1638084092421319, -1.2174375200025294),
                    complex(-9.8419729711608492e-2, -0.9217030074705066),
                    complex(-1.5742564828419501, -0.99686730476185714),
                    complex(1.4777738731392653, 0.47151670288681441),
                    complex(-1.3963398509841487, 1.3006511215030607),
                ],
            },
        },
        case2: {
            desc: 'trivial(alpha=0, beta=1) uplo=l',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                beta: complex(1, 0),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).packedLower(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })().toArr(),
            },
        },
        case3: {
            desc: 'trivial(alpha=0, beta=0) uplo=l',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0).packedLower(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: Array.from({ length: 6 }).fill(complex(0, 0)),
            },
        },
        case4: {
            desc: 'uplo=l, incxy(-1,-1), n=6, beta(1,0), alpha(1,0)',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                alpha: complex(1, 0),
                beta: complex(1, 0),
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                y: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    return v;
                })(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                y: [
                    complex(-1.7225049442869467, 1.9573690212447206),
                    complex(1.7317793933958838, 1.0092722795253204),
                    complex(1.0515376006342987, -0.5220273266242319),
                    complex(0.70976844124876681, -2.1452626845418177),
                    complex(3.9215297000845251e-2, 1.94406543885259),
                    complex(-1.9918115746842791, -1.2672249544431136),
                ],
            },
        },
    },
    chpmvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: 1,
                incy: 1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [0],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'y has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
                y: [0],
            },
        },
        case3: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [0],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                incx: -1,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                incx: 0,
                incy: -1,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },

        case7: {
            desc: 'incy=0',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                incy: 0,
                lda: 6,
                alpha: complex(0, 0),
                beta: complex(0, 0),
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
    },
    // A := alpha*x*x**H + A,
    //CHPR(UPLO,N,ALPHA,X,INCX,AP)
    chpr: {
        case0: {
            desc: 'uplo=u, incx(1), n=6, alpha=0.2, x={6},ap={21})',
            input: {
                uplo: 'u',
                n: 6,
                incx: 1,
                alpha: 0.2,
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: [
                    complex(1.4664377569938685, 0.0),
                    complex(-0.92856705188751221, -0.83204329013824463),
                    complex(-0.29472044110298157, 0.0),
                    complex(-1.2994659767673966, -0.17990848227932066),
                    complex(-0.28946158289909363, 0.26613736152648926),
                    complex(-0.17490710560589112, 0.0),
                    complex(0.43568331003189087, 0.25014132261276245),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.22426788508892059, -0.17262350022792816),
                    complex(0.37739565968513489, 0.0),
                    complex(7.8194520501826237e-2, 0.14056783081885876),
                    complex(0.50360798835754395, -0.94064915180206299),
                    complex(1.0201601128637381, -0.26046736155203748),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-1.0816686910303539, 0.0),
                    complex(-0.26300986834671181, 0.31196009025279664),
                    complex(-0.54288828372955322, 0.24841265380382538),
                    complex(-0.42552053432548659, 0.11192357318741848),
                    complex(-0.64947164058685303, 1.9156390801072121e-2),
                    complex(0.6683753949148884, 0.24178842029114961),
                    complex(1.1698956214925516, 0.0),
                ],
            },
        },
        case1: {
            desc: 'uplo=u, incx(1), n=6, alpha=0.2, x={6},ap={21})',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                alpha: 0.2,
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: [
                    complex(1.2809381210347879, 0.0),
                    complex(-0.98694238844082449, -0.84759324470022923),
                    complex(-0.29472044110298157, -1.1665705442428589),
                    complex(-1.1398672411927047, -0.18069329334156942),
                    complex(-0.28946158289909363, 0.26613736152648926),
                    complex(-0.326518423639795, -0.43068377270992614),
                    complex(0.63861392320075572, 0.0),
                    complex(-1.2375384569168091, 0.61824327707290649),
                    complex(-0.28987718712691601, -0.31726553652934003),
                    complex(0.37739565968513489, -2.223900318145752),
                    complex(7.8194520501826237e-2, 0.14056783081885876),
                    complex(0.50360798835754395, 0.0),
                    complex(1.0857694149017334, -0.11582532525062561),
                    complex(-0.69095385074615479, -0.81496870517730713),
                    complex(-1.2845993041992188, 0.24226348102092743),
                    complex(-0.11139855031280793, 0.0),
                    complex(-0.54288828372955322, 0.24841265380382538),
                    complex(-0.58511926990017871, 0.11270838424966724),
                    complex(-0.64947164058685303, 0.0),
                    complex(0.72675073146820068, 0.25733837485313416),
                    complex(1.3553952574516321, 0.0),
                ],
            },
        },
        case2: {
            desc: 'trivial(alpha=0)',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                alpha: 0,
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                ap: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper().toArr(),
            },
        },
    },
    chprErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                incx: -1,
                alpha: 1,
                x: [0],
                ap: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                n: 6,
                alpha: complex(0, 0),
                incx: 1,
                x: [complex(0, 0)],
                ap: [0],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                n: 6,
                alpha: 0,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0',
            input: {
                uplo: 'l',
                n: -6,
                alpha: 0,
                incx: 1,
                incy: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'l',
                n: 6,
                alpha: 0,
                incx: 0,
                x: [complex(0, 0)],
                y: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
    },
    // x = A*x, or x=A**T*x
    ctbmv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.89893773605253458, 1.1357840738783445),
                    complex(-0.45151985700522917, 0.29850125888377299),
                    complex(-0.71516509118603144, 0.65794798179721126),
                    complex(1.2587631018024554, 1.0538850956187753),
                    complex(0.13566746748087732, -0.6152783090503986),
                    complex(9.4482665585330361e-2, 0.17142208883575871),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.93782339751986754, 0.74413133490913808),
                    complex(-0.45151985700522917, 0.29850125888377299),
                    complex(0.38229682533410747, 0.30850538482755907),
                    complex(1.2587631018024554, 1.0538850956187753),
                    complex(0.26330208478911077, 1.087232850104408),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-1.5857588772442313, 0.33125815615156995),
                    complex(0.54337390686068332, 2.6859788728962997e-2),
                    complex(-2.2279378005181627, 1.263250970435692),
                    complex(-0.89893773605253458, 1.1357840738783445),
                    complex(-0.45151985700522917, 0.29850125888377299),
                    complex(-0.71516509118603144, 0.65794798179721126),
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(0.54337390686068332, 2.6859788728962997e-2),
                    complex(-1.1550642750183755, -0.20168802073263237),
                    complex(-0.89893773605253458, 1.1357840738783445),
                    complex(5.8109940652350645e-2, 1.300847844397051),
                    complex(-0.46766315226263799, 0.14935680643823457),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.61012440075037944, 1.1637949108289476),
                    complex(0.90450126675144427, 0.75266294427094904),
                    complex(-6.9932675283407519e-2, -0.15133458277470524),
                    complex(-0.12669784107884396, 0.87789422422181218),
                    complex(0.3459086754545817, -0.92679784529900511),
                    complex(-0.97741676538200184, -1.0711961221936419),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(0.90450126675144427, 0.75266294427094904),
                    complex(1.0275292412367314, -0.50077717974435743),
                    complex(-0.12669784107884396, 0.87789422422181218),
                    complex(0.47354329276281515, 0.77571331385580145),
                    complex(-1.1896530260810865, -1.5183962364181183),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.17654735879280281, 0.58874771059708308),
                    complex(0.85890551240220603, 1.8168168506139182),
                    complex(-0.77252134008506745, 0.93792187448835485),
                    complex(-0.23040409445892918, -1.2087986124495775),
                    complex(-0.42820606451121151, -0.1429430913374572),
                    complex(-0.36525553403714772, 0.23281314997025904),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.76020145623371604, 1.0296317263052543),
                    complex(0.85890551240220603, 1.8168168506139182),
                    complex(0.30035218541471975, -0.52701711667996953),
                    complex(-0.23040409445892918, -1.2087986124495775),
                    complex(8.1423733146368305e-2, 0.85940349417582085),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=-1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: -1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-4.4081662135385002e-2, -0.47152180004162414),
                    complex(1.4701333254503495, 5.4962968975589277e-2),
                    complex(0.35926976559701229, -1.4420628114452261),
                    complex(0.31909425240906009, -0.28351812300138768),
                    complex(-0.39426565004330849, 2.4148318973788592),
                    complex(-7.279556264438547e-2, -0.38380208237829727),
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.8970153605671185, -0.76156909872986789),
                    complex(0.67646910762479351, -0.79585439150818971),
                    complex(-0.73671647846003907, -1.4532698964078605e-2),
                    complex(0.48292775843658564, 0.79038470663881988),
                    complex(0.51646822672704484, 0.30309190830210242),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.30354170139129621, 0.23178789271847933),
                    complex(-0.55893308970669509, 3.22839923793304),
                    complex(-0.74100210864921268, -0.33250461678806609),
                    complex(-7.262468820471879e-2, -1.1351242539598538),
                    complex(0.29318211076532813, -0.24869434973481597),
                    complex(0.42076612202015307, -0.10280777453071321),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.29188722469987738, -0.61517223133368759),
                    complex(-0.55893308970669509, 3.22839923793304),
                    complex(1.0387957612083905, -0.69225924429680141),
                    complex(-7.262468820471879e-2, -1.1351242539598538),
                    complex(8.7515933527837264e-2, 0.8568022046203525),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case12: {
            desc: 'trivial n=0',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(0.0, 0.0),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(0.0, 0.0),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
    },
    ctbmvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'a has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [0],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'k<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: -5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case9: {
            desc: 'lda<(k+1)',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: 5,
                lda: 4,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case10: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: 5,
                lda: 6,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    // solve: b = A*x, or b=A**T*x, etc
    ctbsv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-8.7140181184168313, 0.53974179877620077),
                    complex(1.8661327422332243, 3.9027499877685838),
                    complex(-3.6580339360462011, 8.8097913188003023),
                    complex(5.6188869330661131, 4.3092860094409744),
                    complex(-0.43765579556553402, -1.3695491241654123),
                    complex(0.0, 0.0),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(3.338660655299988, -0.3048094465802193),
                    complex(-0.45872546160875621, -2.5016404952691627),
                    complex(0.10040131610044789, -1.482719398550286),
                    complex(5.6878452703088889e-2, -1.4873152158496024),
                    complex(2.4240884199307811e-2, 0.90674088189022095),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(0.19681361062735658, 0.14043312545782299),
                    complex(-0.436664155739315, -0.12517747762468392),
                    complex(1.9908862356588675, -0.49306176041025807),
                    complex(-1.1707260559995081, -4.5452002398639015),
                    complex(-0.3624844721825769, 0.31276418173763804),
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.14752502884958729, -0.1638493128589138),
                    complex(1.4363127452141604, -0.57315625854518182),
                    complex(2.2989023266518371, 0.24098851995891635),
                    complex(-1.1879615574537852, -1.4001332192795501),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-7.6101460816946579e-3, -9.4660182789718411e-2),
                    complex(6.8147857987214427e-2, 1.5234806158371664e-3),
                    complex(2.9528998152503227, -0.17560477971933142),
                    complex(3.2401253786020834, -6.2781243245851357),
                    complex(-2.6978707277566238, 7.2418066105756456),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(4)(0, 0);
                    v.s(2)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.90450126675144427, -0.75266294427094904),
                    complex(-0.17021002669464957, -0.33332994457346288),
                    complex(-1.3114892067913817, -1.0467456831262383),
                    complex(2.3201126210298955, 1.3473362542512941),
                    complex(0.56213947293563105, -0.67778084986020015),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.26783294222237042, 5.3412212596831878),
                    complex(-2.4708523091856045, -10.50494649421857),
                    complex(-2.5360982977762716, -1.6617901325948061),
                    complex(7.8824301124023846, -2.4955417126761907),
                    complex(2.7655458811932285, -0.55943960377267044),
                    complex(0.20166478159136628, -4.9273708864892143e-2),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-2.8778446406669276, -3.977772237557331),
                    complex(-3.1812735166230288, -1.3246826397301936),
                    complex(0.87222461178592359, 0.17438318677524328),
                    complex(0.3419771210511211, 1.3634675771268798),
                    complex(0.20611923584205027, 1.1345702378188081),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=-1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: -1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-2.4553590852142104, -5.3900057558999848),
                    complex(3.4605613317430035, 3.7946358385801675),
                    complex(2.3680417148052602, 6.5876435306548492e-2),
                    complex(0.3617001901639596, -8.2795869074406506e-3),
                    complex(0.15114802344991249, 0.41203114326506957),
                    complex(-0.13368054624873069, -0.18738554063649901),
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(0.8970153605671185, 0.76156909872986789),
                    complex(0.71752083895656837, 0.4128360923899177),
                    complex(1.150410910746126, 1.4022465442979499),
                    complex(-0.36504790906346596, -1.4582194308823537),
                    complex(0.95516886987466565, -0.7413058324432491),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(1.2134872328945268, 3.3689875258937509),
                    complex(-0.51454080049378648, -11.238018851608174),
                    complex(1.1080551530027434, -1.5280652728327344),
                    complex(-7.954583528953969, -7.6831586104367161e-2),
                    complex(-2.6639255157882884, 0.22590102372385532),
                    complex(-0.1750596676913814, 0.11158268354623273),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, k=3, n=6, lda=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-2.1629786169498084, -2.8234267947755658),
                    complex(2.6661358253363017, -3.4888111144575142),
                    complex(0.60070653051319522, 0.25805634910918834),
                    complex(0.10717711933379223, 1.3226800014067428),
                    complex(0.20002703546058132, 1.1371715273742764),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case12: {
            desc: 'trivial n=0',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                k: 3,
                lda: 6,
                incx: 1,
                a: matrix_mxn(6, 6, 6),
                x: (() => {
                    const v = vector(6);
                    v.s(2)(0, 0);
                    v.s(4)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(0.0, 0.0),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(0.0, 0.0),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
    },
    ctbsvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'a has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [0],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                k: 5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'k<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: -5,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case9: {
            desc: 'lda<(k+1)',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: 5,
                lda: 4,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case10: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                k: 5,
                lda: 6,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    // CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
    ctpmv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.21480810469340378, 0.69579227075292904),
                    complex(-0.26094969695022141, 2.1899200001689807),
                    complex(-0.33737654149548413, 0.95849359062438277),
                    complex(0.27866515796963282, -3.2357536489400593),
                    complex(-0.44083150014179662, -1.4766224545245905),
                    complex(-0.13564174811293128, -0.31767194384784148),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.21480810469340378, 0.69579227075292904),
                    complex(-0.15915947579584566, 1.766689865635461),
                    complex(0.68550578102442095, 0.65675536979746663),
                    complex(1.8959032428106883, -1.0481669938732936),
                    complex(0.129162241986017, 0.7662624655235204),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(2.9900855818740801e-3, -0.1254483253515466),
                    complex(0.21805658972259057, -0.64067279808502109),
                    complex(0.55454303196612764, 3.3278109254821664e-2),
                    complex(-1.9443462712775996, -0.20935680748047414),
                    complex(-1.4832770655398049, 1.4177745116425551),
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.9473256824028029, -0.22682655495516757),
                    complex(1.7616956826635151, -0.88734821132789865),
                    complex(-1.6881006148507482, 1.4323906112766398),
                    complex(-1.4653889125406279, 1.4596684301016789),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.22095898699729033, 0.20371451287137976),
                    complex(-0.26583054006009044, -9.1245991110092461e-2),
                    complex(-0.45535556765459351, -2.4279571343292985),
                    complex(-1.1230038346504205, -2.3902810640005123),
                    complex(-1.1355567240741828, 1.0540975185653942),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.75705178245981464, -0.39298421193700861),
                    complex(1.161882517186462, -0.24037047926253274),
                    complex(-0.55301009252260691, -0.14739614395240164),
                    complex(-1.1176685710750058, 1.0959914370245181),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-2.4976814036793784, -0.23946010250133121),
                    complex(1.3996314833238857, 0.78835332165500283),
                    complex(2.012987411251268, -1.9013967150345032),
                    complex(-0.36287068806792244, 0.10791794560781565),
                    complex(-0.12708341444083437, -0.87548495323359354),
                    complex(-0.13564174811293128, -0.31767194384784148),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-2.4976814036793784, -0.23946010250133121),
                    complex(1.277472631899097, 0.69428602534440953),
                    complex(2.7422565039314803, -1.4875504719046497),
                    complex(0.84428196262946487, -0.81270837497490467),
                    complex(0.129162241986017, 0.7662624655235204),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=-1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: -1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.99053756426170603, -0.22489208505092062),
                    complex(0.37551077322001747, 1.4700155236647587),
                    complex(1.3389624021526658, 0.18727475284257045),
                    complex(0.24995649269782261, 0.50298068683844832),
                    complex(-0.86662654194144162, 3.1989469298911999e-2),
                    complex(-3.0267127927761095e-2, -0.39887173640357432),
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.64020916573379671, -0.329553690101406),
                    complex(1.0371188810777623, 0.1362709537994804),
                    complex(0.64128318005638185, 1.5767700095488548),
                    complex(-0.78512600062635252, 0.95230004796124257),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.21231009546066959, 1.4470766140474218),
                    complex(-3.5191155812410297, 1.0437484124731138),
                    complex(1.150591221269297, -0.10552584409884025),
                    complex(-0.21019657406961945, -0.75391528212598224),
                    complex(-0.23082261213115096, -0.82038820112026123),
                    complex(-0.13564174811293128, -0.31767194384784148),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-0.21231009546066959, 1.4470766140474218),
                    complex(-3.5314545767923575, 0.89006305075836589),
                    complex(1.0806654353856708, -0.94111695640190396),
                    complex(1.3036001790457676, -0.8687618514814397),
                    complex(-1.2774295761460941e-2, 0.82686750312289536),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case12: {
            desc: 'trivial (n=0)',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
    },
    ctpmvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [0],
                ap: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [0],
            },
        },
        case2: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                incx: 0,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
    },
    ctpsv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-3.3851017924301425, -1.3820753807292496),
                    complex(-1.332498985226668, 0.67573580100024855),
                    complex(-2.0061115606914002, 0.30873159135263567),
                    complex(0.53366394452817878, 0.69694725868851992),
                    complex(3.3264201117576032e-2, -0.76983395647848252),
                    complex(0.0, 0.0),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.15547987558715703, -2.1970472032566155),
                    complex(-0.64651555475774991, -0.84393783126494526),
                    complex(0.41271091050372538, -1.3369077805429912),
                    complex(0.38779544173476843, 0.38706108295362363),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.4232715421031939, -0.26082687970715562),
                    complex(0.39222678598002386, -0.23224612883036766),
                    complex(-1.2997678064087332, -0.94478598552509907),
                    complex(-0.61517805667934433, -0.20792193349925353),
                    complex(6.7048390491203161e-2, -0.61100462807805289),
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.38094570660476545, -0.6227940257180502),
                    complex(0.77065302766175314, 0.23155213977648836),
                    complex(1.9235986450552391, 0.62921553153130372),
                    complex(-0.67065782778448713, -2.0842393309515446),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(0.2011414528077328, -5.1337030960433162e-2),
                    complex(-0.11347507554592454, 1.7911796556818167),
                    complex(0.17118012099267396, 0.4131927253850306),
                    complex(0.11383486554884616, 0.2712765885757436),
                    complex(0.13250679287302938, 0.40677081363532691),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.57121960654775372, -0.45663636873620916),
                    complex(1.0247114696420028, -0.62076677198172869),
                    complex(1.0568841164395277, 1.9635907758261411),
                    complex(0.37156836221836176, -2.7223246591641685),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(12.179271994541113, 14.308616957116691),
                    complex(20.022860582238643, -11.394470245205239),
                    complex(5.3345587897852154, 2.9774530389843976),
                    complex(-3.4438036297245773, 0.90496225731641267),
                    complex(-0.18635839568165291, -1.8489702849007796),
                    complex(-0.10222449470844058, -0.23940899018683348),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-4.6890622533536845, -4.3105166275435485),
                    complex(-4.9095388171859451, 1.5119249268050328),
                    complex(-2.0011605096316609, 0.28698531071705224),
                    complex(1.4229022368946207, 9.6376644440676285e-2),
                    complex(0.15838072700240158, 1.2277112664711085),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=-1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: -1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.39803859497384503, -1.0008794896314415),
                    complex(-0.95461016818598643, -0.10455630604512522),
                    complex(0.1081956677652671, -0.48235598200894436),
                    complex(-2.1147145699491019, -1.3663776876599896),
                    complex(0.8062826920602193, -8.8111846075539443e-2),
                    complex(-0.15013342830572224, -0.16730119413206657),
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.68806222327377164, -0.52006689057181177),
                    complex(1.1537416954878212, -0.99972441484002239),
                    complex(-0.8275770353535008, 7.3583568030202828e-2),
                    complex(1.5551601886103397, -1.4986169107249498),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(-16.223221775482749, 12.881420892412462),
                    complex(-14.48743291807561, -15.055884393595845),
                    complex(-6.4169359859497845, -0.96002124817834023),
                    complex(1.5382948185317014, 4.3925279091258904),
                    complex(-0.48218019458408296, -1.7482432429647807),
                    complex(-0.10222449470844058, -0.23940899018683348),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(2.5181705640049294, -5.4011920508919671),
                    complex(3.8304845810211443, -1.6707772413094739),
                    complex(0.73533020104709457, -0.53087771067682288),
                    complex(0.94106504084051901, 0.16204542951614687),
                    complex(0.30031726474987952, 1.1671062288717335),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
        case12: {
            desc: 'trivial (n=0)',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                incx: 1,
                ap: matrix_mxn(6, 6, 6).slice_used_for_test(1, 6, 1, 6).setLower(0).packedUpper(),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
    },
    ctpsvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [0],
                ap: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'ap has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [0],
            },
        },
        case2: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                incx: 1,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                incx: 0,
                x: [complex(0, 0)],
                ap: [complex(0, 0)],
            },
        },
    },
    //(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    ctrmv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.83150675414885056, 1.9278190152391241),
                    complex(-0.39338349525210226, 2.0694548243325235),
                    complex(-0.40640543562459586, 0.84668404121139496),
                    complex(0.19690462573740852, -3.4126079216408423),
                    complex(-0.42622225763360433, -1.2458980540507965),
                    complex(0.0, 0.0),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.99249198509378145, 1.4460226393835436),
                    complex(-0.29159327409772651, 1.6462246897990038),
                    complex(0.61647688689530922, 0.54494582038447881),
                    complex(1.814142710578464, -1.2250212665740765),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.48802483127278151, 1.2539385477153218),
                    complex(-0.9653771022291957, -0.85185586585022488),
                    complex(-1.2393312224660826, 1.2132717007075515),
                    complex(-2.5654559956894172, -1.3102608123698012),
                    complex(-1.5184106544536302, -3.6775354290217734),
                    complex(3.4744621149400827, -0.57679995107873072),
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.86358688107481996, -1.2750860003837445),
                    complex(-0.2164488999461775, 0.91153347988063538),
                    complex(-0.94821791084836171, 0.87732584269696456),
                    complex(-0.9484169123258166, -1.4346505089736628),
                    complex(3.4744621149400827, -0.57679995107873072),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.48802483127278151, 1.2539385477153218),
                    complex(1.0241460862502108, 2.6733200210011354e-2),
                    complex(0.65454049159191263, -0.82986178215861273),
                    complex(-0.93126308391969648, -2.2538919124303631),
                    complex(-1.0774122839169915, -2.4272069861881223),
                    complex(-1.1294976357766762, 0.95227101262259106),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(1.1259363074045865, -0.3964969343235083),
                    complex(1.6774228141118177, -1.1316000029855289),
                    complex(0.68597500092135899, -6.6305257363597381e-2),
                    complex(-0.50741854178917789, -0.18432206614001156),
                    complex(-1.1294976357766762, 0.95227101262259106),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(2.2251660077830033, 2.323834176613893),
                    complex(0.27152304192522392, -2.3031402183418956),
                    complex(1.0402809982832415, 2.8743260159852113),
                    complex(0.76270793503965795, -2.6553038624241152),
                    complex(-0.42622225763360433, -1.2458980540507965),
                    complex(0.0, 0.0),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(2.0641807768380724, 1.8420378007583125),
                    complex(0.37331326307959944, -2.7263703528754153),
                    complex(2.0631633208031466, 2.5725877951582952),
                    complex(2.3799460198807134, -0.46771720735734945),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-1.1513151820979886, 0.69642190434815632),
                    complex(0.25139557645729682, -1.3313132788431306),
                    complex(0.50669212996350743, -0.56114501912958303),
                    complex(1.1938025271467749, 3.3443289488025765),
                    complex(0.58289064685587677, -0.78703924639420531),
                    complex(-0.23183789661954024, 1.2835785838419986),
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.15897511885275151, -1.4765058693244129),
                    complex(1.2095200637738619, -1.3632466929137586),
                    complex(0.94750069304083251, 0.63502464640267542),
                    complex(0.66981737160638177, 1.5255068342641325),
                    complex(-0.23183789661954024, 1.2835785838419986),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(0.90061994398728085, -1.5995329267787675),
                    complex(5.3055233216601705, 1.9450131306691401),
                    complex(-2.2713161685357997, -1.7723500578825266),
                    complex(0.1066339008687871, 2.604930535910964),
                    complex(5.6844759743704287e-2, -1.3155592146610233),
                    complex(0.0, 0.0),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: -1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-1.2498731094504563, -1.1083395882641662),
                    complex(0.68069188988759555, 0.24910421006826589),
                    complex(0.47157803291335676, -2.7400158851063168),
                    complex(3.6613593566644753, 2.0828881554261853),
                    complex(3.3443933864997035, -3.9680368366347851),
                ],
            },
        },
        case12: {
            desc: '(trivial n=0)',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                incx: -1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
    },
    ctrmvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'a has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                x: [complex(0, 0)],
                a: [0],
            },
        },
        case2: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                lda: 6,
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                incx: 1,
                lda: 6,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                lda: 6,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'lda<n',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                lda: 5,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
    // SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
    ctrsv: {
        case0: {
            desc: 'tr="n", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    // copied from fortran
                    complex(-3.3851017924301425, -1.3820753807292496),
                    complex(-1.332498985226668, 0.67573580100024855),
                    complex(-2.0061115606914002, 0.30873159135263567),
                    complex(0.53366394452817878, 0.69694725868851992),
                    complex(3.3264201117576032e-2, -0.76983395647848252),
                    complex(0.0, 0.0),
                ],
            },
        },
        case1: {
            desc: 'tr="n", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    // copied from fortran
                    complex(0.15547987558715703, -2.1970472032566155),
                    complex(-0.64651555475774991, -0.84393783126494526),
                    complex(0.41271091050372538, -1.3369077805429912),
                    complex(0.38779544173476843, 0.38706108295362363),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
        case2: {
            desc: 'tr="n", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    // copied from fortran
                    complex(0.0, 0.0),
                    complex(0.2011414528077328, -5.1337030960433162e-2),
                    complex(-0.58779672483646239, 1.4444342253366556),
                    complex(-0.59996084429980345, 1.8923968370248347),
                    complex(2.3943121537219199, 1.2601259090017753),
                    complex(-1.1716441125457138, 2.4468331626047548),
                    //
                ],
            },
        },
        case3: {
            desc: 'tr="n", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'n',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    // copied from fortran
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.89736221505014213, -0.55306139391414555),
                    complex(0.74985303977344109, -2.4958541378702988),
                    complex(3.2482401623515491, 3.4359461997861374),
                    complex(-6.1115892109325278, 5.4105978948754974),
                ],
            },
        },
        case4: {
            desc: 'tr="t", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(0.0, 0.0),
                    complex(0.2011414528077328, -5.1337030960433162e-2),
                    complex(-0.11347507554592454, 1.7911796556818167),
                    complex(0.17118012099267396, 0.4131927253850306),
                    complex(0.11383486554884616, 0.2712765885757436),
                    complex(0.13250679287302938, 0.40677081363532691),
                ],
            },
        },
        case5: {
            desc: 'tr="t", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.57121960654775372, -0.45663636873620916),
                    complex(1.0247114696420028, -0.62076677198172869),
                    complex(1.0568841164395277, 1.9635907758261411),
                    complex(0.37156836221836176, -2.7223246591641685),
                ],
            },
        },
        case6: {
            desc: 'tr="t", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(2.5983253692194772, -15.10230045904005),
                    complex(6.3141701681529163, -6.9203837691759773),
                    complex(-3.748870371182194, 6.4167258382816481),
                    complex(0.27337260555960724, 0.88409421189284909),
                    complex(3.3264201117576032e-2, -0.76983395647848252),
                    complex(0.0, 0.0),
                ],
            },
        },
        case7: {
            desc: 'tr="t", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 't',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.42363350716824644, -2.2093867649350898),
                    complex(2.2225892101670914, -1.7538965376185174),
                    complex(-1.1422197242316745, -0.27970358425951164),
                    complex(-0.17800786756748099, -0.37024297626310343),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(0.0, 0.0),
                ],
            },
        },
        case8: {
            desc: 'tr="c", ul="u", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(0.0, 0.0),
                    complex(-0.15262265513448822, 0.1407114065210669),
                    complex(-1.4443838535530182, -0.3990607325181077),
                    complex(-4.6332350725494867e-2, -0.20224298804899468),
                    complex(-1.4838666079245511, -0.94099706199515598),
                    complex(0.39843405102446616, -0.29059808376945545),
                    //
                ],
            },
        },
        case9: {
            desc: 'tr="c", ul="u", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'u',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setLower(0),
                x: (() => {
                    const v = vector(6);
                    v.s(1)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(0.0, 0.0),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.68806222327377164, -0.52006689057181177),
                    complex(1.1537416954878212, -0.99972441484002239),
                    complex(-0.8275770353535008, 7.3583568030202828e-2),
                    complex(1.5551601886103397, -1.4986169107249498),
                ],
            },
        },
        case10: {
            desc: 'tr="c", ul="l", di="n", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'n',
                n: 6,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(1.4997211511167985, 1.6317882003200239),
                    complex(0.66507043863638171, -0.28918592454053232),
                    complex(-0.5731391584186768, -0.55599486430937639),
                    complex(1.2909233609283014e-2, -9.3569865927350199e-2),
                    complex(-0.24941512573253091, -0.72906990246418102),
                    complex(0.0, 0.0),
                ],
            },
        },
        case11: {
            desc: 'tr="c", ul="l", di="u", incx=1, n=6, a={*}, x={*}',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 6,
                incx: -1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: (() => {
                    const v = vector(6);
                    v.s(6)(0, 0);
                    // v.s(3)(0, 0);
                    return v;
                })(),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(1.0115355777646271, 0.66930834493988645),
                    complex(1.6199462693950109, -2.6460116567259719),
                    complex(7.9681750368777013, 2.1980136848408875),
                    complex(-8.2087922632419357, -2.9159109438871607),
                    complex(-7.9036649818434919, -8.2778612243653438),
                ],
            },
        },
        case12: {
            desc: 'trivial (n=0)',
            input: {
                trans: 'c',
                uplo: 'l',
                diag: 'u',
                n: 0,
                incx: 1,
                lda: 6,
                a: matrix_mxn(6, 8).slice_used_for_test(1, 6, 1, 6).setUpper(0),
                x: vector(6),
            },
            expect: {
                x: [
                    //copied from fortran
                    complex(-0.6490100622177124, 0.77214217185974121),
                    complex(-0.11916876584291458, -0.21951562166213989),
                    complex(0.66413569450378418, -0.42481029033660889),
                    complex(1.1009690761566162, -0.41898009181022644),
                    complex(0.14377148449420929, 0.99698686599731445),
                    complex(-0.11775359511375427, -0.27577802538871765),
                ],
            },
        },
    },
    ctrsvErrors: {
        case0: {
            desc: 'x has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                x: [0],
                a: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'a has no imaginary part',
            input: {
                uplo: 'l',
                trans: 't',
                diag: 'u',
                n: 6,
                incx: 1,
                lda: 6,
                x: [complex(0, 0)],
                a: [0],
            },
        },
        case2: {
            desc: 'uplo!="ul"',
            input: {
                uplo: 'x',
                trans: 't',
                diag: 'u',
                n: 6,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'trans!="ntc"',
            input: {
                uplo: 'u',
                trans: 'x',
                diag: 'u',
                lda: 6,
                n: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'diag!="un"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'x',
                n: 6,
                lda: 6,
                incx: 1,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: -6,
                incx: 1,
                lda: 6,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'incx=0',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                lda: 6,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'lda<n',
            input: {
                uplo: 'u',
                trans: 't',
                diag: 'n',
                n: 6,
                lda: 5,
                incx: 0,
                x: [complex(0, 0)],
                a: [complex(0, 0)],
            },
        },
    },
};
