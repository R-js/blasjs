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

          //console.log('before', { are: a.r, aim: a.i });
          cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
          //console.log(`after yr:${y.r}, yi:${y.i}`);

          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
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
        a,
        x,
        y
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();
          //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);


          const call = () => cgbmv(trans, m, n, kl, ku, alpha, aM, lda, sx, incx, beta, sy, incy);
          expect(call).to.throw();
        });
      });
    });
  });
});
