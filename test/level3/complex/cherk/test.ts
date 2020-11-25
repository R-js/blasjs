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
    level3: { cherk },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('cherk', () => {
        describe('data tests', () => {
            const { cherk: testData } = fixture;
            each(testData)(
                (
                    {
                        input: {
                            cmd,
                            uplo,
                            trans,
                            n, // A = m*m marrix , lda >= m
                            k, // columns
                            alpha, // B = alpha * A * B
                            beta,
                            lda, // physical storage
                            ldc, // physical storage
                            a,
                            c,
                        },
                        expect,
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
                        const approx = approximatelyWithPrec(1e-5);
                        if (cmd === 'logDebug') {
                            c.toArr().forEach((cc) => console.log(`     + (${cc.re},${cc.im}),`));
                        }

                        multiplexer(c.toArr(), expect.c)(approx);
                    });
                },
            );
        });

        describe('test errors', () => {
            const { cherkErrors: errors } = fixture;
            each(errors)(
                (
                    {
                        input: {
                            uplo,
                            trans,
                            n, // A = m*m marrix , lda >= m
                            k, // columns
                            alpha, // B = alpha * A * B
                            beta,
                            lda, //physical storage
                            ldc, // physical storage
                            //
                            a: ap,
                            c: cp,
                        },
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const a = fortranMatrixComplex64(ap)(1, 1);
                        //const b = fortranMatrixComplex64(bp)(1, 1);
                        const c = fortranMatrixComplex64(cp)(1, 1);
                        const call = () => cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
                        //call();
                        expect(call).toThrow();
                    });
                },
            );
        });
    });
});
