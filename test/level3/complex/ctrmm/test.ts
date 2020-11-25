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

import * as blas from '../../../../src/lib';
import { approximatelyWithPrec } from '../../../test-helpers';

import { fixture } from './fixtures';

const {
    util: { each, multiplexer, fortranMatrixComplex64 },
    level3: { ctrmm },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('ctrmm', () => {
        describe('data tests', () => {
            const { ctrmm: testData } = fixture;
            each(testData)(
                (
                    {
                        input: {
                            side, //A*B
                            uplo,
                            transA,
                            diag,
                            m, // A = m*m marrix , lda >= m
                            n, // rows in B
                            alpha, // B = alpha * A * B
                            a, // only uses 4x4!!
                            lda, //NBupper 4x4 of A is referenced
                            b, //m*n
                            ldb,
                        },
                        expect,
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        //console.log('a-->');
                        /*a.toArr().forEach(cplx => {
            console.log(`     + (${cplx.re},${cplx.im}),`);
          });*/
                        ctrmm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
                        const approx = approximatelyWithPrec(1e-5);
                        multiplexer(b.toArr(), expect.b)(approx);
                    });
                },
            );
        }); //https://www.youtube.com/watch?v=-XZb3--n4D4
        describe('test errors', () => {
            const { ctrmmErrors: errors } = fixture;
            each(errors)(
                (
                    {
                        input: {
                            side, //A*B
                            uplo,
                            transA,
                            diag,
                            m, // A = m*m marrix , lda >= m
                            n, // rows in B
                            alpha, // B = alpha * A * B
                            lda, //NBupper 4x4 of A is referenced
                            ldb,
                            a, // only uses 4x4!!
                            b,
                        },
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const aM = fortranMatrixComplex64(a)(1, 1);
                        const bM = fortranMatrixComplex64(b)(1, 1);

                        //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);

                        const call = () => ctrmm(side, uplo, transA, diag, m, n, alpha, aM, lda, bM, ldb);
                        //call()
                        expect(call).toThrow();
                    });
                },
            );
        });
    });
});
