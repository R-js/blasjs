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
    level3: { cher2k },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('cher2k', () => {
        describe('data tests', () => {
            const { cher2k: testData } = fixture;
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
                            lda, //physical storage
                            ldb, // physical storage
                            ldc, // physical storage
                            //
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
                        cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
                        const approx = approximatelyWithPrec(1e-5);
                        if (cmd === 'logDebug') {
                            c.toArr().forEach((cc) => console.log(`     + (${cc.re},${cc.im}),`));
                        }

                        multiplexer(c.toArr(), expect.c)(approx);
                    });
                },
            );
        }); //https://www.youtube.com/watch?v=-XZb3--n4D4
        describe('test errors', () => {
            const { cher2kErrors: errors } = fixture;
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
                            ldb, // physical storage
                            ldc, // physical storage
                            //
                            a: ap,
                            b: bp,
                            c: cp,
                        },
                        desc,
                    },
                    key,
                ) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const a = fortranMatrixComplex64(ap)(1, 1);
                        const b = fortranMatrixComplex64(bp)(1, 1);
                        const c = fortranMatrixComplex64(cp)(1, 1);
                        const call = () => cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
                        //call();
                        expect(call).toThrow();
                    });
                },
            );
        });
    });
});
