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
    level3: { csymm },
} = blas;

describe('blas level 3 single/double complex', function n() {
    describe('csymm', () => {
        describe('data tests', () => {
            const { csymm: testData } = fixture;
            each(testData)(
                (
                    {
                        input: {
                            cmd,
                            side,
                            uplo,
                            m,
                            n,
                            alpha,
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
                        csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

                        if (cmd === 'logDebug') {
                            c.toArr().forEach((cc) => console.log(`     + (${cc.re},${cc.im}),`));
                        }

                        const approx = approximatelyWithPrec(1e-5);

                        multiplexer(c.toArr(), expect.c)(approx);
                    });
                },
            );
        }); //https://www.youtube.com/watch?v=-XZb3--n4D4
        describe('test errors', () => {
            const { csymmErrors: errors } = fixture;
            each(errors)(
                (
                    {
                        input: {
                            side,
                            uplo,
                            m,
                            n,
                            alpha,
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
                        const call = () => csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
                        //call();
                        expect(call).toThrow();
                    });
                },
            );
        });
    });
});
