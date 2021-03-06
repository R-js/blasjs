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

import { complex } from '../../../test-helpers';

import { matrix_mxn } from '../../../matrices';

export const fixture = {
    dgemm: {
        case0: {
            desc: '(trivial) n=0',
            input: {
                cmd: 'xdebug',
                trA: 'n',
                trB: 'n',
                m: 6,
                n: 0,
                k: 4,
                lda: 8,
                ldb: 8,
                ldc: 8,
                beta: 0.2,
                alpha: 0.3,
                a: matrix_mxn(8, 8).real(),
                b: matrix_mxn(8, 8).real(),
                c: matrix_mxn(8, 8).real(),
            },
            expect: {
                c: matrix_mxn(8, 8).real().toArr(),
            },
        },
        case1: {
            desc: '(trivial m=0)',
            input: {
                trA: 'n',
                trB: 'n',
                m: 0,
                n: 8,
                k: 4,
                lda: 8,
                ldb: 8,
                ldc: 8,
                beta: 0,
                alpha: 0,
                a: matrix_mxn(8, 8).real(),
                b: matrix_mxn(8, 8).real(),
                c: matrix_mxn(8, 8).real(),
            },
            expect: {
                c: matrix_mxn(8, 8).real().toArr(),
            },
        },
        case2: {
            desc: '(trivial) alpha=0, beta=1',
            input: {
                trA: 'n',
                trB: 'n',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: 1,
                alpha: 0,
                //bandmatrix_mxn_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: matrix_mxn(8, 8).real(),
                b: matrix_mxn(8, 8).real(),
                c: matrix_mxn(8, 8).real(),
            },
            expect: {
                c: matrix_mxn(8, 8).real().toArr(),
            },
        },
        case3: {
            desc: '(near trivial), alpha=0, beta=0',
            input: {
                cmd: 'xdebug',
                trA: 'n',
                trB: 'c',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: 0,
                alpha: 0,
                //bandmatrix_mxn_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: matrix_mxn(8, 8).real(),
                b: matrix_mxn(8, 8).real(),
                c: matrix_mxn(8, 8).real(),
            },
            expect: {
                c: (() => {
                    const m = matrix_mxn(8, 8).real();
                    for (let i = 1; i <= 8; i++) {
                        m.setCol(i, 1, 6, 0);
                    }
                    return m.toArr();
                })(),
            },
        },
        case4: {
            desc: '( near trivial) alpha=0 beta !=1 and beta != 0',
            input: {
                trA: 'n',
                trB: 'c',
                m: 6, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0.75,
                alpha: 0,
                //bandmatrix_mxn_ku_kl(n = 6, m = 6, lda = m, kl = 3, ku = 2)
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    0.94721567630767822,
                    -0.24467501789331436,
                    0.99734947085380554,
                    0.95432201027870178,
                    0.31098107993602753,
                    -1.1549625098705292,
                    -0.69642528891563416,
                    -0.22104033082723618,
                    -4.3253795010969043e-3,
                    1.8034899830818176,
                    0.57269507646560669,
                    -0.59925694763660431,
                    -0.86074277758598328,
                    -0.21709618717432022,
                    -0.22441133111715317,
                    -0.30863311886787415,
                    0.18916759639978409,
                    -0.66894082725048065,
                    0.32676248252391815,
                    -0.92815384268760681,
                    -0.16820091381669044,
                    0.28304674476385117,
                    0.10000227391719818,
                    0.60314212739467621,
                    -4.2830080725252628e-2,
                    0.37770599126815796,
                    0.81432706117630005,
                    -0.51821538805961609,
                    -0.96344947814941406,
                    3.5044628195464611e-2,
                    -0.17677991464734077,
                    -0.40716621279716492,
                    -0.32498274743556976,
                    -0.48710373044013977,
                    0.54506304860115051,
                    0.86393380165100098,
                ],
            },
        },
        case5: {
            desc: 'trA="n", trB="n",s m=4,n=6, k=3, ld[x]=6, alpha=0.3, beta=-1.2',
            input: {
                trA: 'n',
                trB: 'n',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: -1.2,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    -1.4039963092403844,
                    0.18124124487022381,
                    -1.2107208849353177,
                    -1.444320742259589,
                    0.41464143991470337,
                    -1.5399500131607056,
                    0.84654511690969958,
                    0.4711022660577785,
                    -0.36249415168348781,
                    -3.4519430694935669,
                    0.76359343528747559,
                    -0.79900926351547241,
                    1.1260112382918237,
                    0.51125157946674815,
                    -7.142821162956231e-2,
                    -0.11615811188891931,
                    0.25222346186637878,
                    -0.89192110300064087,
                    6.421047087526327e-2,
                    1.5712993271415543,
                    0.4652053193175591,
                    -1.1516306540304853,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.46722627045877629,
                    -0.73755420551909645,
                    -1.4240404210621318,
                    1.0366043860301484,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    0.49396185302075668,
                    0.76016266322751813,
                    0.46577487040588361,
                    0.35124613688394185,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case6: {
            desc: 'trA="n", trB="n", m=4, n=6, k=3, ld[x]=6, alpha=0.3, beta=0',
            input: {
                trA: 'n',
                trB: 'n',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    0.11154883307425167,
                    -0.21023879931509784,
                    0.38503833184054292,
                    8.2594534860493257e-2,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.2677353896328426,
                    0.11743772268083484,
                    -0.36941475916024308,
                    -0.56635898189985345,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -0.25117726057030149,
                    0.16389766618523197,
                    -0.43048635568469557,
                    -0.60997112169988321,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.5870304636885314,
                    8.6253119830951214e-2,
                    0.19608384651692823,
                    -0.69875584441269589,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.5357544023422437,
                    -0.13322459547614354,
                    -0.12111707140653165,
                    0.20745973218751551,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    0.21111397834564721,
                    0.10869669686512365,
                    -5.419754615287467e-2,
                    -0.42811986278950293,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case7: {
            desc: 'trA="n", trB="n" m=4, n=6, k=3, alpha(0.3), beta(1)',
            input: {
                trA: 'n',
                trB: 'n',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 1,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    1.3745030681511559,
                    -0.53647215650618363,
                    1.7148376263122835,
                    1.3550238818987623,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -1.1963024415203547,
                    -0.17728271842214668,
                    -0.37518193182837228,
                    1.8382943288759033,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -1.3988342973516126,
                    -0.12556391671386166,
                    -0.72970146384089973,
                    -1.0214819468570489,
                    0.25222346186637878,
                    -0.89192110300064087,
                    1.0227137737204224,
                    -1.1512853370858578,
                    -2.8184038571992344e-2,
                    -0.321360184727561,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.59286117664258053,
                    0.37038339288140038,
                    0.96465234349520179,
                    -0.48349411855863927,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -2.4592574517473847e-2,
                    -0.4341915868644296,
                    -0.48750787606696766,
                    -1.0775915033763559,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case8: {
            desc: 'trA="n", trB="n" m=4, n=6, k=3, alpha(0.3), beta(1)',
            input: {
                trA: 'c',
                trB: 'n',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 1,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    2.3039086064868313,
                    -0.65151114030011925,
                    0.80392857434356968,
                    1.4691522943492057,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -1.2538448349965456,
                    -9.9813792277925195e-3,
                    0.34004655561627861,
                    2.3930913417469486,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -1.6735277569094822,
                    5.635214538531419e-2,
                    0.14791521851026912,
                    -0.43391807697828583,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.63240625734282752,
                    -1.2491004259456175,
                    -0.24667513691004092,
                    0.90888091978565244,
                    0.13333636522293091,
                    0.80418950319290161,
                    0.30512477050075026,
                    0.47311061141257815,
                    0.96423497396955082,
                    -0.95843924086395926,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.44474478900959474,
                    -0.42847770252990852,
                    -0.26611774831495649,
                    -0.44957283992730146,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case9: {
            desc: 'trA="n", trB="n" m=4, n=6, k=3, alpha(0.3), beta(1)',
            input: {
                trA: 'c',
                trB: 'n',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    1.040954371409927,
                    -0.32527778310903349,
                    -0.52587072012817104,
                    0.19672294731093659,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.32527778310903349,
                    0.28473906187518905,
                    0.34581372828440782,
                    -1.156196902880839e-2,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -0.52587072012817104,
                    0.34581372828440782,
                    0.44713032666647334,
                    -2.2407251821120321e-2,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.19672294731093659,
                    -1.156196902880839e-2,
                    -2.2407251821120321e-2,
                    0.53148526010051755,
                    0.13333636522293091,
                    0.80418950319290161,
                    0.3622315448010871,
                    -3.0497376944965775e-2,
                    -0.1215344409321826,
                    -0.26748539011780442,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.20903823614647371,
                    0.11441058119964469,
                    0.1671925815991365,
                    0.1998988006595516,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case10: {
            desc: 'trA="n", trB="c" m=4, n=6, k=3, alpha(0.3), beta(0)',
            input: {
                cmd: 'xdebug',
                trA: 'n',
                trB: 'c',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0,
                alpha: 0.3,
                a: (() => {
                    const m = matrix_mxn(6, 6);
                    const m2 = m.real();
                    //console.log(m2);
                    //process.exit(1);
                    return m2;
                })(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    1.1323220981414599,
                    5.8155756369945882e-2,
                    0.60846817867763747,
                    -4.6076554446066703e-2,
                    0.41464143991470337,
                    -1.5399500131607056,
                    5.8155756369945882e-2,
                    8.3122908219772917e-2,
                    -0.10365417583453679,
                    -0.30140785416513954,
                    0.76359343528747559,
                    -0.79900926351547241,
                    0.60846817867763747,
                    -0.10365417583453679,
                    0.5573787535903566,
                    0.5404013774225942,
                    0.25222346186637878,
                    -0.89192110300064087,
                    -4.6076554446066703e-2,
                    -0.30140785416513954,
                    0.5404013774225942,
                    2.2712326344858456,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.14245217765629503,
                    -0.12999764483403706,
                    0.14145512021468692,
                    0.67799604713505268,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -5.3799957506489804e-2,
                    0.29881330774903153,
                    -0.53290206537630114,
                    -1.0541348433374778,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case11: {
            desc: 'trA="n", trB="c" m=4, n=6, k=3, alpha(0.3), beta(1)',
            input: {
                trA: 'n',
                trB: 'c',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 1,
                alpha: 0.3,
                a: matrix_mxn(6, 6).real(),
                b: matrix_mxn(6, 6).real(),
                c: matrix_mxn(6, 6).real(),
            },
            expect: {
                c: [
                    2.395276333218364,
                    -0.26807760082113996,
                    1.9382674731493781,
                    1.2263527925922022,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.87041129551756624,
                    -0.21159753288320862,
                    -0.10942134850266599,
                    2.1032454566106171,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -0.53918885810367356,
                    -0.39311575873363042,
                    0.25816364543415238,
                    0.12889055226542873,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.38960675558582425,
                    -1.5389463110819486,
                    0.31613349233367366,
                    2.6486282941709804,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.19955895195663187,
                    0.37361034352350686,
                    1.2272245351164204,
                    -1.2957803611101958e-2,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.28950651036961084,
                    -0.24407497598052169,
                    -0.96621239529039427,
                    -1.7036064839243306,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case12: {
            //
            desc: 'trA="n", trB="c", m=4, n=6, k=3, alpha(0.3), beta(0.5)',
            input: {
                trA: 'n',
                trB: 'c',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0.5, //
                alpha: 0.3, //
                a: matrix_mxn(6, 6).real(), //
                b: matrix_mxn(6, 6).real(), //
                c: matrix_mxn(6, 6).real(), //
            },
            expect: {
                c: [
                    1.763799215679912,
                    -0.10496092222559705,
                    1.2733678259135077,
                    0.5901381190730679,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.40612776957381014,
                    -6.4237312331717866e-2,
                    -0.10653776216860139,
                    0.90091880122273893,
                    0.76359343528747559,
                    -0.79900926351547241,
                    3.4639660286981871e-2,
                    -0.2483849672840836,
                    0.40777119951225449,
                    0.33464596484401149,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.17176510056987879,
                    -0.92017708262354403,
                    0.42826743487813396,
                    2.459930464328413,
                    0.13333636522293091,
                    0.80418950319290161,
                    -0.17100556480646345,
                    0.12180634934473492,
                    0.68433982766555357,
                    0.3325191217619754,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.17165323393805032,
                    2.736916588425492e-2,
                    -0.74955723033334765,
                    -1.3788706636309043,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case13: {
            //
            desc: 'trA="c", trB="c", m=4, n=6, k=3, alpha(0.3), beta(0.5)',
            input: {
                trA: 'c',
                trB: 'c',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0.5, //
                alpha: 0.3, //
                a: matrix_mxn(6, 6).real(), //
                b: matrix_mxn(6, 6).real(), //
                c: matrix_mxn(6, 6).real(), //
            },
            expect: {
                c: [
                    0.74302595061270371,
                    -0.43085206822838551,
                    0.41372238666556888,
                    1.223245137207666,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.67452232525885392,
                    -2.9922497870655915e-2,
                    0.16101407985116736,
                    1.2885797752188297,
                    0.76359343528747559,
                    -0.79900926351547241,
                    -0.18879018655011265,
                    -0.51414555060978995,
                    -0.58009390976279773,
                    -9.6715660616545351e-3,
                    0.25222346186637878,
                    -0.89192110300064087,
                    0.30043618987643866,
                    -1.1851282103582581,
                    -0.72210506424434351,
                    -0.51005801457012856,
                    0.13333636522293091,
                    0.80418950319290161,
                    0.15443764840970681,
                    6.8346909955205232e-2,
                    0.31117483789830569,
                    -0.59174378784732873,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.9789430016359888,
                    0.22972839886126217,
                    0.46299282935915376,
                    -0.16936478468399807,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
        case14: {
            //
            desc: 'trA="c", trB="c", m=4, n=6, k=3, alpha(0.3), beta(0.5)',
            input: {
                trA: 'c',
                trB: 'c',
                m: 4, // A(M,K), C(M,N)
                n: 6, // B(K,N), A(M,K)
                k: 3,
                lda: 6, // lda >= M
                ldb: 6, // ldb >= K
                ldc: 6, // ldc >= M
                beta: 0, //
                alpha: 0.3, //
                a: matrix_mxn(6, 6).real(), //
                b: matrix_mxn(6, 6).real(), //
                c: matrix_mxn(6, 6).real(), //
            },
            expect: {
                c: [
                    0.11154883307425162,
                    -0.2677353896328426,
                    -0.25117726057030149,
                    0.5870304636885314,
                    0.41464143991470337,
                    -1.5399500131607056,
                    -0.21023879931509784,
                    0.11743772268083487,
                    0.16389766618523197,
                    8.6253119830951214e-2,
                    0.76359343528747559,
                    -0.79900926351547241,
                    0.38503833184054287,
                    -0.36941475916024313,
                    -0.43048635568469557,
                    0.19608384651692823,
                    0.25222346186637878,
                    -0.89192110300064087,
                    8.2594534860493229e-2,
                    -0.56635898189985345,
                    -0.60997112169988321,
                    -0.698755844412696,
                    0.13333636522293091,
                    0.80418950319290161,
                    0.18299103555987523,
                    -0.18345708422356674,
                    -0.23170986955256098,
                    -0.24626686247425131,
                    -1.2845993041992188,
                    4.6726170927286148e-2,
                    -0.86108972520442828,
                    0.50117254072603878,
                    0.67964799431620027,
                    0.15537103560942844,
                    0.72675073146820068,
                    1.151911735534668,
                ],
            },
        },
    },
    dgemmErrors: {
        case3: {
            desc: 'transA!="ntc"',
            input: {
                trA: 'x',
                trB: 't',
                m: 6, // A(M,K), C(M,N)
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
            },
        },
        case4: {
            desc: 'transB!="ntc"',
            input: {
                trA: 'c',
                trB: 'x',
                m: 6, // A(M,K), C(M,N)
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
            },
        },
        case5: {
            desc: 'n<0"',
            input: {
                trA: 'c',
                trB: 'c',
                m: 6, // A(M,K), C(M,N)
                n: -8, // B(K,N), A(M,K)
                k: 4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            },
        },
        case6: {
            desc: 'k<0',
            input: {
                trA: 'c',
                trB: 'c',
                m: 6, // A(M,K), C(M,N)
                n: 8, // B(K,N), A(M,K)
                k: -4,
                lda: 8, // lda >= M
                ldb: 8, // ldb >= K
                ldc: 8, // ldc >= M
                beta: complex(0.2, 0),
                alpha: complex(0, 0),
                a: [complex(0, 0)],
                b: [complex(0, 0)],
                c: [complex(0, 0)],
            },
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
            },
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
            },
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
            },
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
            },
        },
    },
};
