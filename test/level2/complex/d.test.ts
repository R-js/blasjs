import { assert, expect } from 'chai';
import * as blas from '../../../src/lib';
import { approximately, approximatelyWithPrec } from '../../test-helpers';
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
  level2: {
    cgbmv
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 2 single/double complex', function n() {

  describe('cgbmv', () => {

    describe('data tests', () => {
      const { cgbmv: testData } = fixture;

      each(testData)(({ input: {
        trans,
        m,
        n,
        kl,
        ku,
        lda,
        incx,
        incy,
        beta,
        alpha,
        //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
        a,
        x,
        y
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
          console.log({ y })
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe.skip('test errors', () => {
      const { cgbmvErrors: errors } = fixture;
      each(errors)(({ input: {
        trans,
        m,
        n,
        kl,
        ku,
        lda,
        incx,
        incy,
        beta,
        alpha,
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(muxCmplx([]))(0, 0);
          const sx = fortranArrComplex64(muxCmplx([]))();
          const sy = fortranArrComplex64(muxCmplx([]))();

          const call = () => cgbmv(trans, m, n, kl, ku, alpha, aM, lda, sx, incx, beta, sy, incy);
          expect(call).to.throw();
        });
      });
    });
  });
});
