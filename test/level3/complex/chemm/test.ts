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
    chemm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('chemm', () => {

    describe('data tests', () => {
      const { chemm: testData } = fixture;
      each(testData)(({
        input: {
          side,
          uplo,
          m, //c(mxn)
          n, //c(mxn)
          alpha,
          beta,
          lda,
          ldb, //ldb >= M
          ldc, //ldc >= M
          a,
          b,
          c,
          debugB,
          debugC
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          if (debugB) {
            console.log('b=');
            b.toArr().forEach(
              cpx => console.log(`    +(${cpx.re},\t${cpx.im})`)
            );
          }
          chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
          if (debugC) {
            console.log('c=');
            c.toArr().forEach(
              cpx => console.log(`    +(${cpx.re},\t${cpx.im})`)
            );
          }
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(c.toArr(), expect.c)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { chemmErrors: errors } = fixture;
      each(errors)(({ input: {
        side,
        uplo,
        m, //c(mxn)
        n, //c(mxn)
        alpha,
        beta,
        lda,
        ldb, //ldb >= M
        ldc, //ldc >= M
        a,
        b,
        c
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const bM = fortranMatrixComplex64(b)(1, 1);
          const cM = fortranMatrixComplex64(c)(1, 1);

          const call = () => chemm(side, uplo, m, n, alpha, aM, lda, bM, ldb, beta, cM, ldc);
          //call()
          expect(call).to.throw();
        });
      });
    });
  });
});
