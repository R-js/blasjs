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
    csyrk
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('csyrk', () => {
    describe('data tests', () => {
      //  SUBROUTINE csyrk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
      //   C := alpha*A*B**T + alpha*B*A**T + beta*C,
      //
      const { csyrk: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
          uplo,
          trans,
          n,
          k,
          alpha,
          beta,
          lda,
          ldc,
          a,
          c
        }, expect, desc
      }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          // ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
          csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
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
      const { csyrkErrors: errors } = fixture;
      each(errors)(({ input: {
        uplo,
        trans,
        n, // A = m*m marrix , lda >= m
        k, // columns
        alpha,  // B = alpha * A * B
        beta,
        lda, //physical storage
        ldc, //physical storage
        //
        a: ap,
        c: cp
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const a = fortranMatrixComplex64(ap)(1, 1);
          const c = fortranMatrixComplex64(cp)(1, 1);
          const call = () => csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
