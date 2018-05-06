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

import {
    Complex,
    complex,
    fortranArrComplex64,
    fortranMatrixComplex32,
    fortranMatrixComplex64,
    Matrix,
    muxCmplx
} from '../src/lib/f_func';

const { max, min } = Math;

const data = {
    re36: [
        1.2629542848807933,
        -0.3262333607056494,
        1.3297992629225006,
        1.2724293214294047,
        0.4146414344564082,
        -1.5399500419037095,
        -0.9285670347135381,
        -0.2947204467905602,
        -0.005767172747536955,
        2.404653388857951,
        0.7635934611404596,
        -0.7990092489893682,
        -1.1476570092363514,
        -0.28946157368822334,
        -0.29921511789731614,
        -0.411510832795067,
        0.2522234481561323,
        -0.8919211272845686,
        0.43568329935571865,
        -1.237538421929958,
        -0.22426788527830935,
        0.37739564598170106,
        0.1333363608148414,
        0.8041895097449078,
        -0.057106774383808755,
        0.5036079722337261,
        1.085769362145687,
        -0.6909538396968303,
        -1.2845993538721883,
        0.04672617218835198,
        -0.23570655643950122,
        -0.5428882550102544,
        -0.4333103174567822,
        -0.6494716467962331,
        0.726750747385451,
        1.1519117540872,
        0.9921603654457979,
        -0.42951310949188126,
        1.2383041008533804,
        -0.2793462818542693,
        1.7579030898107073,
        0.5607460908880562,
        -0.4527839725531578,
        -0.8320432961178319,
        -1.166570547084707,
        -1.0655905803882961,
        -1.563782051071005,
        1.1565369971501793,
        0.8320471285723897,
        -0.22732869142475534,
        0.2661373616721048,
        -0.3767027185836281,
        2.4413646288945894,
        -0.7953391172553718,
        -0.054877473711578625,
        0.2501413228541527,
        0.6182432935662469,
        -0.17262350264585732,
        -2.2239002740099374,
        -1.263614384970583,
        0.3587288959713519,
        -0.011045478465663564,
        -0.9406491626186084,
        -0.11582532215695436,
        -0.8149687088699175,
        0.24226348085968588,
        -1.4250983947324998,
        0.36594112304921983,
        0.2484126488725964,
        0.06528818167162072,
        0.01915639166027384,
        0.2573383771555333
    ],
    im36: [
        0.9921603654457979,
        -0.42951310949188126,
        1.2383041008533804,
        -0.2793462818542693,
        1.7579030898107073,
        0.5607460908880562,
        -0.4527839725531578,
        -0.8320432961178319,
        -1.166570547084707,
        -1.0655905803882961,
        -1.563782051071005,
        1.1565369971501793,
        0.8320471285723897,
        -0.22732869142475534,
        0.2661373616721048,
        -0.3767027185836281,
        2.4413646288945894,
        -0.7953391172553718,
        -0.054877473711578625,
        0.2501413228541527,
        0.6182432935662469,
        -0.17262350264585732,
        -2.2239002740099374,
        -1.263614384970583,
        0.3587288959713519,
        -0.011045478465663564,
        -0.9406491626186084,
        -0.11582532215695436,
        -0.8149687088699175,
        0.24226348085968588,
        -1.4250983947324998,
        0.36594112304921983,
        0.2484126488725964,
        0.06528818167162072,
        0.01915639166027384,
        0.2573383771555333,
        1.2629542848807933,
        -0.3262333607056494,
        1.3297992629225006,
        1.2724293214294047,
        0.4146414344564082,
        -1.5399500419037095,
        -0.9285670347135381,
        -0.2947204467905602,
        -0.005767172747536955,
        2.404653388857951,
        0.7635934611404596,
        -0.7990092489893682,
        -1.1476570092363514,
        -0.28946157368822334,
        -0.29921511789731614,
        -0.411510832795067,
        0.2522234481561323,
        -0.8919211272845686,
        0.43568329935571865,
        -1.237538421929958,
        -0.22426788527830935,
        0.37739564598170106,
        0.1333363608148414,
        0.8041895097449078,
        -0.057106774383808755,
        0.5036079722337261,
        1.085769362145687,
        -0.6909538396968303,
        -1.2845993538721883,
        0.04672617218835198,
        -0.23570655643950122,
        -0.5428882550102544,
        -0.4333103174567822,
        -0.6494716467962331,
        0.726750747385451,
        1.1519117540872,
    ],
    rnorm10: [
        -0.6490100777088978,
        -0.11916876241803812,
        0.6641356998941105,
        1.100969102194087,
        0.14377148075806995,
        -0.11775359816595128,
        -0.9120683669483379,
        -1.4375862408299789,
        -0.7970895250719646,
        1.2540831064499711
    ],
    inorm10: [
        0.7721421858045301,
        -0.21951562675343952,
        -0.4248102833772871,
        -0.418980099421959,
        0.9969868609091059,
        -0.27577802908802723,
        1.2560188173061,
        0.6466743904953449,
        1.299312302563431,
        -0.873262111744435
    ]
};

export function vector(n = 6) {
    const nre = new Array(n).fill(0);
    const nim = new Array(n).fill(0);
    let cursor = 0;
    for (let i = 1; i <= n; i++) {
        nre[cursor] = data.rnorm10[cursor];
        nim[cursor] = data.inorm10[cursor];
        cursor++;
    }
    return fortranArrComplex64(muxCmplx(nre, nim))();
}

export function bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4): Matrix {
    if (lda < m) {
        throw new Error(`lda<m, ${lda}<${m}`);
    }
    //console.log({ n, m, lda, kl, ku });
    const nre = new Array(n * lda).fill(0);
    const nim = new Array(n * lda).fill(0);
    //kl
    let cursor = 0;
    for (let j = 1; j <= n; j++) {
        //upper
        for (let i = max(j - ku, 1); i <= min(j - 1, m); i++) {
            nre[(j - 1) * lda + i - 1] = data.re36[cursor];
            nim[(j - 1) * lda + i - 1] = data.im36[cursor++];
        }
        //include diagonal here aswell
        for (let i = j; i <= min(j + kl, m); i++) {
            nre[(j - 1) * lda + i - 1] = data.re36[cursor];
            nim[(j - 1) * lda + i - 1] = data.im36[cursor++];
        }
    }
    return fortranMatrixComplex64(muxCmplx(nre, nim))(lda, n);
}

export function matrix_mxn(lda: number, n: number, m: number = lda, float = 64) {

    if (lda < m) {
        throw new Error(`lda<m, ${lda}<${m}`);
    }

    const nre = new Array(n * lda).fill(0);
    const nim = new Array(n * lda).fill(0);

    let cursor = 0;
    for (let j = 1; j <= n; j++) {
        //upper
        for (let i = 1; i <= m; i++) {
            nre[(j - 1) * lda + i - 1] = data.re36[cursor];
            nim[(j - 1) * lda + i - 1] = data.im36[cursor++];
        }
    }
    //
    const func = float === 64 ? fortranMatrixComplex64 : fortranMatrixComplex32;
    //
    return func(muxCmplx(nre, nim))(lda, n);
}

export function diagonal_nxn(n: number): Matrix {

    const re = new Array(n * n);
    const im = new Array(n * n);

    for (let j = 1; j <= n; j++) {
        for (let i = 1; i <= n; i++) {
            if (i === j) {
                re[(j - 1) * n - 1 + i] = 1;
                im[(j - 1) * n - 1 + i] = 0;
                continue;
            }
            re[(j - 1) * n - 1 + i] = 0;
            im[(j - 1) * n - 1 + i] = 0;
        }
    }
    return fortranMatrixComplex64(muxCmplx(re, im))(n, n);
}
