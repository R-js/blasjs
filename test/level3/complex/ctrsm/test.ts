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
    level3: { ctrsm },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('ctrsm', () => {
        describe('data tests', () => {
            const { ctrsm: testData } = fixture;
            each(testData)(
                ({ input: { cmd, side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb }, expect, desc }, key) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        ctrsm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
                        if (cmd === 'logDebug') {
                            b.toArr().forEach((cplx) => {
                                console.log(`     + (${cplx.re},${cplx.im}),`);
                            });
                        }
                        const approx = approximatelyWithPrec(1e-5);
                        multiplexer(b.toArr(), expect.b)(approx);
                    });
                },
            );
        }); //https://www.youtube.com/watch?v=-XZb3--n4D4
        describe('test errors', () => {
            const { ctrsmErrors: errors } = fixture;
            each(errors)(({ input: { side, uplo, transA, diag, m, n, alpha, a: aP, lda, b: bP, ldb }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const a = fortranMatrixComplex64(aP)(1, 1);
                    const b = fortranMatrixComplex64(bP)(1, 1);
                    const call = () => ctrsm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
                    //call()
                    expect(call).toThrow();
                });
            });
        });
    });
});
