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
    csymm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('csymm', () => {
    describe('data tests', () => {
      const { csymm: testData } = fixture;
      each(testData)(({
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
          c
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

          if (cmd === 'logDebug') {
            c.toArr().forEach(cc =>
              console.log(`     + (${cc.re},${cc.im}),`));
          }

          const approx = approximatelyWithPrec(1E-5);

          multiplexer(c.toArr(), expect.c)(approx);
        });
      });
    }); //https://www.youtube.com/watch?v=-XZb3--n4D4
    describe('test errors', () => {
      const { csymmErrors: errors } = fixture;
      each(errors)(({
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
          c: cp
        }, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const a = fortranMatrixComplex64(ap)(1, 1);
          const b = fortranMatrixComplex64(bp)(1, 1);
          const c = fortranMatrixComplex64(cp)(1, 1);
          const call = () => csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
