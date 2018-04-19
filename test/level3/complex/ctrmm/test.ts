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
    ctrmm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {
  describe('ctrmm', () => {
    describe('data tests', () => {
      const { ctrmm: testData } = fixture;
      each(testData)(({
        input: {
          side, //A*B
          uplo,
          transA,
          diag,
          m, // A = m*m marrix , lda >= m
          n, // rows in B
          alpha,  // B = alpha * A * B
          a, // only uses 4x4!!
          lda, //NBupper 4x4 of A is referenced
          b, //m*n
          ldb
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log('a-->');
          /*a.toArr().forEach(cplx => {
            console.log(`     + (${cplx.re},${cplx.im}),`);
          });*/
          ctrmm(side, uplo, transA, diag, m, n, alpha, a, lda, b, ldb);
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(b.toArr(), expect.b)(approx);
        });
      });
    }); //https://www.youtube.com/watch?v=-XZb3--n4D4
    describe('test errors', () => {
      const { ctrmmErrors: errors } = fixture;
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

          const call = () => ctrmm(side, uplo, transA, diag, m, n, alpha, aM, lda, bM, ldb);
          //call()
          expect(call).to.throw();
        });
      });
    });
  });
});
