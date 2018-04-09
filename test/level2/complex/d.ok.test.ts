import { assert, expect } from 'chai';
import * as blas from '../../../src/lib';
import { Matrix } from '../../../src/lib/f_func';
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
    cgbmv,
    cgemv,
    cgerc,
    cgeru
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

          // console.log('before', { are: a.r, aim: a.i });
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

  describe('cgemv', () => {

    describe('data tests', () => {
      const { cgemv: testData } = fixture;

      each(testData)(({ input: {
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          //console.log(`before yr:${y.r}, yi:${y.i}`);
          cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgemvErrors: errors } = fixture;
      each(errors)(({ input: {
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgemv(trans, m, n, alpha, aM, lda, sx, incx, beta, sy, incy);
          expect(call).to.throw();
        });
      });
    });
  });
  describe('cgerc', () => {

    describe('data tests', () => {
      const { cgerc: testData } = fixture;

      each(testData)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log('a:', a.toArr());
          //console.log(`before y:${JSON.stringify(y.toArr())}`);
          //console.log(`before x:${JSON.stringify(x.toArr())}`);
          cgerc(m, n, alpha, x, incx, y, incy, a, lda);
          //a.toArr().forEach(cplx => {
          //  console.log(`(${cplx.re},${cplx.im})`);
          //});
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgercErrors: errors } = fixture;
      each(errors)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgerc(m, n, alpha, sx, incx, sy, incy, aM, lda);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
  describe('cgeru', () => {

    describe('data tests', () => {
      const { cgeru: testData } = fixture;

      each(testData)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log('a:', a.toArr());
          //console.log(`before y:${JSON.stringify(y.toArr())}`);
          //console.log(`before x:${JSON.stringify(x.toArr())}`);
          cgeru(m, n, alpha, x, incx, y, incy, a, lda);
          //a.toArr().forEach(cplx => {
          //  console.log(`(${cplx.re},${cplx.im})`);
          //});
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgeruErrors: errors } = fixture;
      each(errors)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgeru(m, n, alpha, sx, incx, sy, incy, aM, lda);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
