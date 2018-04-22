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
    cherk
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('cherk', () => {
    describe('data tests', () => {
      const { cherk: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
          uplo,
          trans,
          n, // A = m*m marrix , lda >= m
          k, // columns
          alpha,  // B = alpha * A * B
          beta,
          lda, // physical storage
          ldb, // physical storage
          ldc, // physical storage
          a,
          b,
          c
        }, expect, desc
      }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
          const approx = approximatelyWithPrec(1E-5);
          if (cmd === 'logDebug') {
            c.toArr().forEach(cc =>
              console.log(`     + (${cc.re},${cc.im}),`));
          }

          multiplexer(c.toArr(), expect.c)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cherkErrors: errors } = fixture;
      each(errors)(({ input: {
        uplo,
        trans,
        n, // A = m*m marrix , lda >= m
        k, // columns
        alpha,  // B = alpha * A * B
        beta,
        lda, //physical storage
        ldb, // physical storage
        ldc, // physical storage
        //
        a: ap,
        b: bp,
        c: cp
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const a = fortranMatrixComplex64(ap)(1, 1);
          const b = fortranMatrixComplex64(bp)(1, 1);
          const c = fortranMatrixComplex64(cp)(1, 1);
          const call = () => cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
