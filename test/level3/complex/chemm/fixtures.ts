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

import { complex } from '../../../../src/lib/f_func';

import { matrix_mxn } from '../../../matrices';

export const fixture = {
    // C := alpha*op( A )*op( B ) + beta*C,
    // CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    chemm: {
        case0: {
            desc: 'trivial (n=0)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 0, // c(mxn)
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: matrix_mxn(6, 6).toArr(),
            },
        },
        case1: {
            desc: 'trivial (m=0)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 0, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: matrix_mxn(6, 6).toArr(),
            },
        },
        case2: {
            desc: 'trivial (alpha=0+i0, beta=1+i0)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0, 0),
                beta: complex(1, 0),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: matrix_mxn(6, 6).toArr(),
            },
        },
        case4: {
            desc: '(near trivial), alpha=0+i0, beta=0+i0',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                //
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0, 0),
                beta: complex(0, 0),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: (() => {
                    const m = matrix_mxn(6, 6);
                    for (let i = 1; i <= 6; i++) {
                        m.setCol(i, 1, 4, 0);
                    }
                    return m.toArr();
                })(),
            },
        },
        case5: {
            desc: 'trivial (alpha=0+i0, beta=-0.1+i0.5)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0, 0),
                beta: complex(-0.1, 0.5),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: [
                    complex(-0.62237560749053955, 0.53226113319396973),
                    complex(0.24737988460992777, -0.1201653682745798),
                    complex(-0.75213200052052631, 0.54106923157229403),
                    complex(1.2430207025365014e-2, 0.66414930266044481),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(0.31924869258488275, -0.41900512806656365),
                    complex(0.44549368961858815, -6.4155890297825202e-2),
                    complex(0.58386198939683609, 0.11377346982854686),
                    complex(0.29232997535965843, 1.3088857189798215),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-0.30125784702948177, -0.65703323011405335),
                    complex(0.14261050267437336, -0.12199792232017415),
                    complex(-0.10314716950175895, -0.17622129062732661),
                    complex(0.22950244607256476, -0.1680851394285221),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(-1.6129594200377628e-2, 0.22332940258812528),
                    complex(-1.3168137706243144e-3, -0.6437833610924204),
                    complex(-0.286694849693376, -0.17395827117300522),
                    complex(4.8572183583087236e-2, 0.20596018012258932),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-0.17365376531018797, -6.442627624977959e-2),
                    complex(-4.4838060444009287e-2, 0.25290854202366692),
                    complex(0.36174763279293565, 0.63694962403274946),
                    complex(0.12700804872952975, -0.33389439267542165),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(0.7361198652322688, 2.4656567610946922e-2),
                    complex(-0.12868173935579019, -0.30803825611761382),
                    complex(-8.087529326482068e-2, -0.24149643070759275),
                    complex(3.2303075715112151e-2, -0.33126463825298574),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case6: {
            desc: 'side="l", uplo="u",m=4,n=6,alpha=0.2-i0.8, beta=-0.4+i0.5)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(-0.4, 0.5),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: [
                    complex(-2.9350657340050365e-2, 0.74442972472321212),
                    complex(-0.45724841254510279, 2.4720085014727919),
                    complex(-2.6679770806555987, 0.5913904722523653),
                    complex(-9.3815933352068392e-2, -0.661869673646718),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(1.5935559894779265, -0.62128381917021469),
                    complex(2.1737031164892846, 2.1413157835605925),
                    complex(3.6351520910646191, 0.304199262683018),
                    complex(0.13840214848770671, 1.8006383908576788),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.32836568409053613, 0.32203405755460113),
                    complex(-0.16239848402013413, -1.4478841158354525),
                    complex(0.29888586547997031, -2.3536534270582283),
                    complex(0.82847839806672074, 0.20658982019875485),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(-0.39135653961213285, -1.2753244343960868),
                    complex(0.60691671473854047, -0.46120495921820964),
                    complex(-0.70341197232683361, -0.34856996829018666),
                    complex(0.35560598553536038, -1.6447638016528807),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(1.3449207395328406, 1.3215570006272412),
                    complex(-0.46409757273693786, 0.10578712406097923),
                    complex(-0.23249122621298632, 0.74794125234075148),
                    complex(-0.38751770034385613, 1.013323247499361),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-1.0967932443017383, -0.34955354858902032),
                    complex(1.1465268689873731, -1.8341626752337215),
                    complex(0.81729097243390925, 0.39141011200385634),
                    complex(-0.19520613456232078, -1.0893009966163467),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case7: {
            desc: 'side="l", uplo="u",m=4,n=6,alpha=0.2-i0.8, beta=0+i0)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(0, 0),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: [
                    complex(0.97191123416181879, 0.50981676505314755),
                    complex(-0.80249830577073311, 2.4633199387844824),
                    complex(-1.5169052858488845, 0.42181248767080021),
                    complex(0.27548266942224975, -1.4098228637310937),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(0.99573717717577193, -0.3381138847352283),
                    complex(1.6397932932222989, 1.9558586830974209),
                    complex(3.0495599498415631, -0.15954537563336316),
                    complex(0.56746817711041331, 0.17207548110202797),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.28532641495519784, 1.2286814228388478),
                    complex(-0.39184746285822969, -1.3940848009033966),
                    complex(0.31226850119727223, -2.0975909267832287),
                    complex(0.47552270260740953, 0.26166414017709494),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(-0.24452195045453179, -1.5151170797007516),
                    complex(0.23697198590189406, 0.2576207997762584),
                    complex(-0.48399748916268936, 1.0861288768453325e-2),
                    complex(0.42025250154490373, -1.9025110326155357),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(1.5014424722976405, 1.4936019441758543),
                    complex(-0.26817711353436224, -0.15043506149737251),
                    complex(-0.26850802968163445, -0.17120312143764616),
                    complex(-0.72181190738603673, 1.3124700420818154),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-1.9036250764466338, -0.80173964832748956),
                    complex(1.1123421207973976, -1.416342076357596),
                    complex(0.7681731647874539, 0.70743033996308791),
                    complex(-0.42235070535685332, -0.73844990448468328),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case8: {
            desc: 'side="l", uplo="l",m=4,n=6,alpha=0.2-i0.8, beta=6+i4)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'l',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(0.6, 0.4),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: [
                    complex(2.7403184111072991, -4.1567812522220677),
                    complex(1.5701776230823246, -1.7028121207367017),
                    complex(2.9217292794908212, 2.6645885235228102),
                    complex(0.37030167414141746, -0.67010759673956355),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-2.2630835590522227, -2.1229320066196156),
                    complex(2.4842119823959798, -7.2380204875169287),
                    complex(0.24931288619950953, 1.1393326242309958),
                    complex(-0.15626431993015122, 1.8841563199065459),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-0.46557598774341352, 1.6616017947432198),
                    complex(-1.2238551365692814, -0.46964373631685929),
                    complex(-0.81409954211000601, 1.9643136086659696),
                    complex(0.28587285557125508, 1.4882791703898399),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.89436110606475971, -1.271865303880453),
                    complex(-1.1453435725176404, -0.90558360620368383),
                    complex(1.3523464275293064, -6.98456127100614e-2),
                    complex(1.1498371510409668, 1.7082210529850559),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-2.1478995758838577, 0.37940110303243435),
                    complex(0.35644753206828794, 0.56659829235487935),
                    complex(0.83133135676970071, 0.27374976655088035),
                    complex(-0.46041159244318619, -0.50905345399265156),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-0.95665587081220771, -0.11330529120488797),
                    complex(-1.3978267600056178, 1.7207265051741496),
                    complex(-1.4299226926401771, -2.2074907397804298),
                    complex(-0.895860749554656, 1.0081050580370876),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case9: {
            desc: 'side="l", uplo="l",m=4,n=6,alpha=0.2-i0.8, beta=0+0i)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'l',
                uplo: 'l',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(0, 0),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
            },
            expect: {
                c: [
                    complex(2.3794099978183687, -5.2572592053675464),
                    complex(1.594112403891133, -1.3146109075897425),
                    complex(2.6191713337571962, 1.3896862853742649),
                    complex(-0.50489448098386491, -1.0114715641286953),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-1.8870568972897215, -1.4798347863197581),
                    complex(2.3282269330697893, -6.6209063153986607),
                    complex(-0.21385503471255873, 1.8415818476915469),
                    complex(-2.0252926180947801, 1.5616493786935888),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.55583710858124047, 1.6214363335574582),
                    complex(-1.1411096564460139, -0.21746188326819893),
                    complex(-0.52811552388554139, 1.9243172304508096),
                    complex(0.38209826787649442, 1.8789051474191738),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.61100011936940624, -1.4132121442392671),
                    complex(-0.30276393832627546, -0.56065301559214753),
                    complex(1.7342044784438015, -0.35108443832156466),
                    complex(0.8543503451119846, 1.6608368911139704),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-1.9701439535436795, 0.18700647314980803),
                    complex(4.9864535667231709e-2, 0.37178238124408441),
                    complex(-0.19638998438563049, 0.40383150762655817),
                    complex(-9.2169396312492768e-2, -0.16317671166391745),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-1.3852713096446827, 0.83603641683594687),
                    complex(-0.92571731981307137, 1.7183171309318661),
                    complex(-1.0705714213586131, -2.1832142034369761),
                    complex(-0.48006247787969669, 1.2287208094127564),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case10: {
            desc: 'side="r", uplo="u",m=4,n=6,alpha=0.2-i0.8, beta=0.1+0.5i)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'r',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(0.1, 0.5),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
                //debugC: 'j'
            },
            expect: {
                c: [
                    complex(2.0492507380581135, -4.6850680675356742),
                    complex(-0.25775830863319293, 0.19202397956818207),
                    complex(1.3794382927760527, -0.4405424154070231),
                    complex(1.7237475261652158, 0.76227571044567122),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(-0.81369106407956837, 0.34668457065587932),
                    complex(1.5414087275787951, -2.1579996421612142),
                    complex(-1.4436299452257024, -1.1949543586973261),
                    complex(0.52941132056679197, 2.9833438731343485),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(-0.43239486946219902, 1.7003452332800477),
                    complex(1.8862347506739041, -1.4426502124799074),
                    complex(-0.20421421491147029, -0.14624259107327331),
                    complex(-3.515448911150143e-2, 1.9023132526066009),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(0.80601718477197659, -1.1614707583064985),
                    complex(0.2797683068346497, -0.61066513425175373),
                    complex(1.7846243759415663, 0.12124436054600474),
                    complex(1.1617834020943869, 1.0987454093262699),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(-1.8790441320046147e-2, 1.5761153832527943),
                    complex(-0.25422556413653397, 0.7370997670744015),
                    complex(1.2693742279114388, 2.3969554131871549),
                    complex(0.23838117192449293, -0.77212650684833373),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-2.2885058185422174, -2.2605005049474509),
                    complex(0.78383108046051031, -0.52310022989716165),
                    complex(-1.7265234423765672, -2.5859374695898527),
                    complex(-0.75308659307088432, 1.9179323367452152),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
        case11: {
            desc: 'side="r", uplo="l",m=4,n=6,alpha=0.2-i0.8, beta=0+0i)',
            input: {
                // SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                // SIDE = 'L' or 'l', C := alpha*A*B + beta*C,
                //
                // SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
                // A is mxm or nxn matrix
                // B is mxn matrix
                // C is mxn matrix
                side: 'r',
                uplo: 'l',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, -0.8),
                beta: complex(0, 0),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: matrix_mxn(6, 6),
                b: matrix_mxn(6, 6),
                c: matrix_mxn(6, 6),
                // debugC: 'j'
            },
            expect: {
                c: [
                    complex(2.7620838935196597, 0.52572550560245579),
                    complex(-0.38722128429009905, 1.1122493341536011),
                    complex(3.5326335051329956, -1.6429894187198397),
                    complex(-2.9614312869456008, -2.0975914928727786),
                    complex(0.41464143991470337, 1.7579030990600586),
                    complex(-1.5399500131607056, 0.56074607372283936),
                    complex(2.6066196946963776, -2.5408032465493644),
                    complex(0.18710810744380824, 1.9314235066316638),
                    complex(-0.64836486108789093, 0.78587559107639837),
                    complex(0.4478236020595785, 0.94720972056341235),
                    complex(0.76359343528747559, -1.5637820959091187),
                    complex(-0.79900926351547241, 1.1565370559692383),
                    complex(0.1094756542499874, -1.4785510653310989),
                    complex(1.3309259797098396, -1.0493629461557186),
                    complex(3.4759585804621675, -6.1280154572816539),
                    complex(0.36200306504389429, -2.6446968059067788),
                    complex(0.25222346186637878, 2.4413645267486572),
                    complex(-0.89192110300064087, -0.79533910751342773),
                    complex(-1.804938012002399, 0.65158667891350619),
                    complex(-2.4740051538807633, -0.35407474560864172),
                    complex(-2.0103238075393381, -0.95470719060530307),
                    complex(3.5191883793917236, -6.1574489307043727),
                    complex(0.13333636522293091, -2.223900318145752),
                    complex(0.80418950319290161, -1.2636144161224365),
                    complex(0.74505125745013312, -3.7206741947987183),
                    complex(-2.6813185514918385, 0.96547808598946283),
                    complex(-0.53371050520175256, -2.0393952416895602),
                    complex(2.6098071711445425, -2.5701179676264876),
                    complex(-1.2845993041992188, -0.81496870517730713),
                    complex(4.6726170927286148e-2, 0.24226348102092743),
                    complex(-2.9630782063318857, -0.251201084004402),
                    complex(0.21875545102578733, 1.7054655666082565),
                    complex(-2.4531243858344594, 3.0252959544891533),
                    complex(-2.2812439006104843, 3.5862454770076102),
                    complex(0.72675073146820068, 1.9156390801072121e-2),
                    complex(1.151911735534668, 0.25733837485313416),
                ],
            },
        },
    },
    chemmErrors: {
        case0: {
            desc: 'a has no imaginary part',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                a: [0],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            },
        },
        case1: {
            desc: 'b has no imaginary part',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                b: [0],
                a: [complex(0, 0)],
                c: [complex(0, 0)],
            },
        },
        case2: {
            desc: 'c has no imaginary part',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4, // c(mxn)
                n: 6, // c(mxn)
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6, // ldb >= M
                ldc: 6, // ldc >= M
                c: [0],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case3: {
            desc: 'side!="lr"',
            input: {
                side: 'x',
                uplo: 'u',
                m: 4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case4: {
            desc: 'uplo!="ul"',
            input: {
                side: 'l',
                uplo: 'x',
                m: 4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4,
                n: -6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'm<0',
            input: {
                side: 'l',
                uplo: 'u',
                m: -4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case7: {
            desc: 'lda<max(1, {si = "l"} => m)',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 3,
                ldb: 6,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case8: {
            desc: 'ldb< max(1,m)',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 3,
                ldc: 6,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
        case9: {
            desc: 'ldc< max(1,m)',
            input: {
                side: 'l',
                uplo: 'u',
                m: 4,
                n: 6,
                alpha: complex(0.2, 0.8),
                beta: complex(0.7, -0.3),
                lda: 6,
                ldb: 6,
                ldc: 3,
                c: [complex(0, 0)],
                a: [complex(0, 0)],
                b: [complex(0, 0)],
            },
        },
    },
};
