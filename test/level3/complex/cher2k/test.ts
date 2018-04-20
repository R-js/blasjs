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
    cher2k
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('cher2k', () => {
    describe('data tests', () => {
      const { cher2k: testData } = fixture;
      each(testData)(({
        input: {
          debug,
          cmd,
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
          a,
          b,
          c
        }, expect, desc
      }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
          const approx = approximatelyWithPrec(1E-5);
          if (debug || cmd === 'debug') {
            c.toArr().forEach(cc =>
              console.log(`     + (${cc.re},${cc.im}),`));
          }

          multiplexer(c.toArr(), expect.c)(approx);
        });
      });
    }); //https://www.youtube.com/watch?v=-XZb3--n4D4
    /* describe('test errors', () => {
       const { cher2kErrors: errors } = fixture;
       each(errors)(({ input: {
         side, //A*B
         uplo,
         transA,
         diag,
         m, // A = m*m marrix , lda >= m
         n, // rows in B
         alpha,  // B = alpha * A * B
         lda, //NBupper 4x4 of A is referenced
         ldb,
         a, // only uses 4x4!!
         b,
       }, desc }, key) => {
         it(`[${key}]/[${desc}]`, function t() {
 
           const aM = fortranMatrixComplex64(a)(1, 1);
           const bM = fortranMatrixComplex64(b)(1, 1);
 
           //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);
 
           const call = () => cher2k(side, uplo, transA, diag, m, n, alpha, aM, lda, bM, ldb);
           //call()
           expect(call).to.throw();
         });
       });
     });*/
  });
});
