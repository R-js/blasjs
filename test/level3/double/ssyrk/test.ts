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
import { approximatelyWithPrec, real, each, multiplexer, fortranMatrixComplex64 } from '../../../test-helpers';
const {
    level3: { ssyrk },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('ssyrk', () => {
        describe('data tests', () => {
            const { ssyrk: testData } = fixture;
            each(testData)(({ input: { cmd, trans, uplo, k, n, lda, ldc, beta, alpha, a, c }, expect, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);

                    if (cmd === 'debug') {
                        console.log(`c actual:\n ${c.toArr().join(',\n')}\n\n`);
                        console.log(`c expected:\n${(<Array<number>>real(expect.c)).join(',\n')}`);
                        const r = real(expect.c);
                        const d = c.toArr().map((v, i) => v - r[i]);
                        console.log(d);
                    }
                    const approx = approximatelyWithPrec(1e-5);
                    multiplexer(c.toArr(), real(expect.c))(approx);
                });
            });
        });

        describe('test errors', () => {
            const { ssyrkErrors: errors } = fixture;
            each(errors)(({ input: { uplo, trans, n, k, alpha, beta, lda, ldc }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

                    const cM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

                    const call = () => ssyrk(uplo, trans, n, k, alpha, aM, lda, beta, cM, ldc);
                    //call();
                    expect(call).toThrow();
                });
            });
        });
    });
});
