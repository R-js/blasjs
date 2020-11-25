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
import { fixture } from './fixtures';
import * as blas from '../../../../src/lib';
import { approximatelyWithPrec, each, multiplexer, fortranMatrixComplex64 } from '../../../test-helpers';

const {
    level3: { cgemm },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('cgemm', () => {
        describe('data tests', () => {
            const { cgemm: testData } = fixture;
            each(testData)(
                (
                    {
                        input: {
                            trA,
                            trB,
                            m, // A(M,K), C(M,N)
                            n, // B(K,N), A(M,K)
                            k,
                            lda, // lda >= M
                            ldb, // ldb >= K
                            ldc, // ldc >= M
                            beta,
                            alpha,
                            //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
                            a,
                            b,
                            c,
                        },
                        expect,
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        cgemm(trA, trB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
                        //console.log(`after yr:${y.r}, yi:${y.i}`);
                        //  console.log(c.r.length);
                        /* c.toArr().forEach(
             cpx => console.log(`     +(${cpx.re},\t${cpx.im})`)
           );*/
                        const approx = approximatelyWithPrec(1e-5);
                        multiplexer(c.toArr(), expect.c)(approx);
                    });
                },
            );
        });

        describe('test errors', () => {
            const { cgemmErrors: errors } = fixture;
            each(errors)(
                (
                    {
                        input: {
                            trA,
                            trB,
                            m, // A(M,K), C(M,N)
                            n, // B(K,N), A(M,K)
                            k,
                            lda, // lda >= M
                            ldb, // ldb >= K
                            ldc, // ldc >= M
                            beta,
                            alpha,
                            //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
                            a,
                            b,
                            c,
                        },
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const aM = fortranMatrixComplex64(a)(1, 1);
                        const bM = fortranMatrixComplex64(b)(1, 1);
                        const cM = fortranMatrixComplex64(c)(1, 1);

                        //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);

                        const call = () => cgemm(trA, trB, m, n, k, alpha, aM, lda, bM, ldb, beta, cM, ldc);
                        //call()
                        expect(call).toThrow();
                    });
                },
            );
        });
    });
});
