import { assert, expect } from 'chai';
import * as blas from '../../../../src/lib';
import { Matrix } from '../../../../src/lib/f_func';
import { matrix_mxn } from '../../../matrices';

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
    ctrsm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('ctrsm', () => {
    describe('data tests', () => {
      const { ctrsm: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
          side,
          uplo,
          transA,
          diag,
          m,
          n,
          alpha,
          a,
          lda,
          b,
          ldb
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          ctrsm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
          if (cmd === 'logDebug') {
            b.toArr().forEach(cplx => {
              console.log(`     + (${cplx.re},${cplx.im}),`);
            });
          }
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(b.toArr(), expect.b)(approx);
        });
      });
    }); //https://www.youtube.com/watch?v=-XZb3--n4D4
    describe('test errors', () => {
      const { ctrsmErrors: errors } = fixture;
      each(errors)(({ input: {
        side,
        uplo,
        transA,
        diag,
        m,
        n,
        alpha,
        a: aP,
        lda,
        b: bP,
        ldb
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const a = fortranMatrixComplex64(aP)(1, 1);
          const b = fortranMatrixComplex64(bP)(1, 1);
          const call = () => ctrsm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
          //call()
          expect(call).to.throw();
        });
      });
    });
  });
});
