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
import { real } from '../../../../src/lib/f_func';
import { approximatelyWithPrec } from '../../../test-helpers';
import { fixture } from './fixtures';

const {
    util: { each, multiplexer, fortranMatrixComplex64 },
    level3: { strmm },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('strmm', () => {
        describe('data tests', () => {
            const { strmm: testData } = fixture;
            each(testData)(
                ({ input: { cmd, side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb }, expect, desc }, key) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        // DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
                        strmm(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb);

                        if (cmd === 'debug') {
                            console.log(`b actual:\n ${b.toArr().join(',\n')}\n\n`);
                            console.log(`b expected:\n${(<Array<number>>real(expect.b)).join(',\n')}`);
                            const r = real(expect.b);
                            const d = b.toArr().map((v, i) => v - r[i]);
                            console.log(d);
                        }
                        const approx = approximatelyWithPrec(1e-5);
                        multiplexer(b.toArr(), real(expect.b))(approx);
                    });
                },
            );
        });

        describe('test errors', () => {
            const { strmmErrors: errors } = fixture;
            each(errors)(({ input: { side, uplo, trans, diag, m, n, alpha, lda, ldb }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

                    const bM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

                    const call = () => strmm(side, uplo, trans, diag, m, n, alpha, aM, lda, bM, ldb);
                    //call();
                    expect(call).toThrow();
                });
            });
        });
    });
});
