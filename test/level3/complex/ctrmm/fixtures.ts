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
        case1: {
            desc: 'si="l", ul="u",trA="n", diag="u", alpha=(0.2,0.8)',
            input: {
                side: 'l', //A*B
                uplo: 'u',
                transA: 'n',
                diag: 'u',
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
                    complex(-4.5878528270079100E-002, -1.6328225075105753),
                    complex(-1.4266711201397058, 1.3912352824123402),
                    complex(0.47796290633027105, 0.96207423451653451),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-0.76589425391122323, 0.98606678861933172),
                    complex(-1.7952358832719872, -2.3286074808032682),
                    complex(-0.42450625879478310, 0.20280603525803245),
                    complex(1.3334031780571216, 1.7106045501024205),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-0.59457311678842317, -0.75219244957874709),
                    complex(-9.2946084873448687E-003, 0.39429243933295322),
                    complex(-7.1771211921762917E-002, 4.0015163033017165E-002),
                    complex(0.21906001894265437, -0.40454921133144506),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.89429376859854393, 1.1330330926899785),
                    complex(-0.63534501209259031, -1.0787114552230488),
                    complex(-0.75266008085101577, 1.6310068561511820E-002),
                    complex(0.21357793530192870, 0.26739183168697100),
                    complex(0.13333636522293091, -2.2239003181457520),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-1.9492069659996545, -0.61817558648055604),
                    complex(0.18397794917513985, 0.68461246633872708),
                    complex(1.3359479596257482, 0.78149896871615987),
                    complex(-4.5530510627187937E-002, -0.57592815422904708),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148E-002, 0.24226348102092743),
                    complex(1.7695659025695811, -6.3875064777400237E-002),
                    complex(-3.4141288894177890E-002, 0.37099623206041554),
                    complex(6.8604956493863833E-002, -0.29596711801438530),
                    complex(-0.18212487372942077, -0.50651968429266925),
                    complex(0.72675073146820068, 1.9156390801072121E-002),
                    complex(1.1519117355346680, 0.25733837485313416),
                ]
            },
        },
        case2: {
            desc: 'si="l", ul="l",trA="n", diag="n", alpha=(0.2,0.8)',
            input: {
                side: 'l', //A*B
                uplo: 'l',
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
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(-0.37066642795896643, -0.12937536991789439),
                    complex(-0.53847920794539261, -0.90801649455416311),
                    complex(1.4384450317723125, -1.1169017108178481),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(-0.92898483439414603, -0.19018934561877915),
                    complex(0.57848959639076525, -1.2758597080561640),
                    complex(1.3474196179958526, -2.5611815252913859),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.0000000000000000, -0.0000000000000000),
                    complex(-0.29787973916147625, 0.60822070451199339),
                    complex(-0.45228594812363454, -2.2680385476264697),
                    complex(-1.2911708889627647, -1.5158614406578568),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(-0.54795774185350954, 0.48306801661471988),
                    complex(-1.1615060277364126, 1.0118885323693230),
                    complex(-1.4892507275872222, -1.1002656652229008),
                    complex(0.13333636522293091, -2.2239003181457520),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-0.0000000000000000, 0.0000000000000000),
                    complex(0.40963463948849754, -8.9577905665065391E-002),
                    complex(-0.43354575166980625, -0.41052366517593775),
                    complex(5.8696244842979695E-002, 0.10846525402097236),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148E-002, 0.24226348102092743),
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(-0.74215273248827995, 0.12542282209793743),
                    complex(1.7852967366485175, 1.2067822355121136),
                    complex(-0.24207421890203995, -1.2786381411846650),
                    complex(0.72675073146820068, 1.9156390801072121E-002),
                    complex(1.1519117355346680, 0.25733837485313416),
                ]
            },
        },
        case3: {
            desc: 'si="l", ul="l",trA="n", diag="u", alpha=(0.2,0.8)',
            input: {
                side: 'l', //A*B
                uplo: 'l',
                transA: 'n',
                diag: 'u',
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
                    complex(0.0000000000000000, 0.0000000000000000),
                    complex(0.27836381015723655, -0.34688931028379955),
                    complex(-1.1309596860986888, 0.98876983179347988),
                    complex(1.5699501899216979, -0.43540248683290739),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(0.17651376987269174, -0.83341044833377653),
                    complex(0.19114517061722991, -0.20611370435252130),
                    complex(1.7261704251386947, -1.8330466031303638),
                    complex(1.8823116790993515, -1.2659749840239980),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-0.89516910456372933, -0.75171621965939117),
                    complex(9.3132696400312798E-002, 0.35268675223505719),
                    complex(-0.85619069320613450, -2.4372907052016850),
                    complex(-1.0849484695051297, -1.7299206282738928),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.13103864388225439, 0.33757115807490279),
                    complex(-0.34537869936083304, -1.1064123009432514),
                    complex(-1.8772066118506514, 1.0830042752573195),
                    complex(-1.4024342919909158, -0.89691777948657769),
                    complex(0.13333636522293091, -2.2239003181457520),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-0.29840446782714819, 2.6060358078187074E-002),
                    complex(0.21810073837373067, 0.52034416997246713),
                    complex(1.0073710161307641, 0.21550738033713684),
                    complex(0.12976748517256848, -0.25796975062319144),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148E-002, 0.24226348102092743),
                    complex(1.0929374410649750, -0.47358493318535277),
                    complex(-0.96129415423971043, -0.67605414721227164),
                    complex(1.3354772095237124, 0.89691338794251474),
                    complex(-0.26802875492832845, -1.6254385282619381),
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
