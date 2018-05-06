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

import { assert, expect } from 'chai';
import * as blas from '../../../../src/lib';
import { Matrix, real } from '../../../../src/lib/f_func';
import {
  approximately,
  approximatelyWithPrec
} from '../../../test-helpers';
import { fixture } from './fixtures';

const {
  util: {
    arrayrify,
    numberPrecision,
    each,
    multiplexer,
    fortranArrComplex64,
    fortranMatrixComplex64,
    complex,
    muxCmplx
  },
  level3: {
    ssyr2k
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('ssyr2k', () => {

    describe('data tests', () => {
      const { ssyr2k: testData } = fixture;
      each(testData)(({
        //DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
        input: {
          cmd,
          uplo,
          trans,
          n,
          k,
          alpha,
          beta,
          lda,
          ldb,
          ldc,
          a,
          b,
          c
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

          if (cmd === 'debug') {
            console.log(`c actual:\n ${c.toArr().join(',\n')}\n\n`);
            console.log(`c expected:\n${(<Array<number>>real(expect.c)).join(',\n')}`);
            const r = real(expect.c);
            const d = c.toArr().map((v, i) => v - r[i]);
            console.log(d);

          }
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(c.toArr(), real(expect.c))(approx);
        });
      });
    });

    describe('test errors', () => {
      const { ssyr2kErrors: errors } = fixture;
      each(errors)(({ input: {
        cmd,
        uplo,
        trans,
        n,
        k,
        alpha,
        beta,
        lda,
        ldb,
        ldc
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const bM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);
          const cM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const call = () => ssyr2k(uplo, trans, n, k, alpha, aM, lda, bM, ldb, beta, cM, ldc);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
