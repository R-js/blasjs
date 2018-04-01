import { assert, expect } from 'chai';
import * as blas from '../../../src/lib';
import { approximitly } from '../../test-helpers';
import { fixture } from './fixtures';


const {
  util: { arrayrify, numberPrecision, each, multiplexer, fortranArrComplex64, complex, muxCmplx },
  level1: { scasum, cswap, csscal, caxpy, ccopy, cdotc, cdotu, crotg, cscal, csrot, icamax }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 1 complex', function n() {
  //

  const precision = numberPrecision(9);

  //
  describe('caxpy', () => {
    //construct tests for dbeta from fixtures

    //abuse as a for-each loop
    //make sure it is an arrow function for `map` but not an arrow function for `it`
    describe('data tests', () => {
      const { caxpy: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const res = fortranArrComplex64(muxCmplx(ore, oim))();

          caxpy(n, ca, cx, incx, cy, incy);
          multiplexer(cy.toArr(), res.toArr())(approximitly);
        });
      });
    });

    describe('error tests', () => {
      const { caxpyErrors: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const call = () => caxpy(n, ca, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });

  });

  describe('ccopy', () => {
    describe('data tests', () => {

      const { ccopy: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const res = fortranArrComplex64(muxCmplx(ore, oim))();

          ccopy(n, cx, incx, cy, incy);
          multiplexer(cy.toArr(), res.toArr())(approximitly);
        });
      });
    });

    describe('error tests', () => {

      const { ccopyErrors: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const call = () => ccopy(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });

  });

  describe('cdotc', () => {

    describe('data tests', () => {
      const { cdotc: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: expected, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const answer = cdotc(n, cx, incx, cy, incy);
          approximitly(answer, expected);
        });
      });
    });

    describe('error tests', () => {

      const { cdotcError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const call = () => cdotc(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });

  });

  describe('cdotu', () => {

    describe('data tests', () => {
      const { cdotu: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: expected, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();

          const answer = cdotu(n, cx, incx, cy, incy);
          approximitly(answer, expected);

        });
      });
    });


    describe('error tests', () => {
      const { cdotuError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const call = () => cdotu(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });
  });

  describe('crotg', () => {
    describe('data tests', () => {
      const { crotg: testData } = fixture;
      each(testData)(({ input: { ca: aIn, cb: bIn }, output: { ca, cb, c: cOut, s: sOut }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          const s = { re: 0, im: 0 };
          const c = { val: 0 };

          crotg(aIn, bIn, c, s); // adjusts aIn
          //console.log({ aIn, ca, bIn, cb, c, s, cOut, sOut });
          approximitly(aIn, ca);
          approximitly(bIn, cb);
          approximitly(c.val, cOut);
          approximitly(s, sOut);
        });
      });
    });

    describe('cscal', () => {
      describe('data tests', () => {
        const { cscal: testData } = fixture;
        each(testData)(({ input: { n, ca, cx: x, incx }, output: { cx: cxOut }, desc }, key) => {
          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
            const result = muxCmplx(cxOut.re, cxOut.im);

            cscal(n, ca, cx, incx); // adjusts aIn

            multiplexer(cx.toArr(), result)(approximitly);
          });
        });
      });

      describe('error tests', () => {
        const { cscalError: testData } = fixture;
        each(testData)(({ input: { n, cx: x, ca, incx }, desc }, key) => {

          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
            const call = () => cscal(n, ca, cx, incx);
            expect(call).to.throw();
          });
        });
      });
    });

    describe('csrot', () => {
      describe('data tests', () => {
        const { csrot: testData } = fixture;

        each(testData)(({ input: { n, cx: x, incx, cy: y, incy, c, s }, output: { cx: cxOut, cy: cyOut }, desc }, key) => {
          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
            const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
            const cxExpected = muxCmplx(cxOut.re, cxOut.im);
            const cyExpected = muxCmplx(cyOut.re, cyOut.im);

            csrot(n, cx, incx, cy, incy, c, s);

            multiplexer(cx.toArr(), cxExpected)(approximitly);
            multiplexer(cy.toArr(), cyExpected)(approximitly);
          });
        });
      });
      describe('error tests', () => {
        const { csrotError: testData } = fixture;
        each(testData)(({ input: { n, cx: x, incx, cy: y, incy, c, s }, desc }, key) => {

          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
            const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();

            const call = () => csrot(n, cx, incx, cy, incy, c, s);

            expect(call).to.throw();
          });
        });
      });
    });
    describe('csscal', () => {
      describe('data tests', () => {
        const { csscal: testData } = fixture;

        each(testData)(({ input: { n, cx: x, incx, sa, }, output: { cx: cxOut }, desc }, key) => {
          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
            const cxExpected = muxCmplx(cxOut.re, cxOut.im);

            csscal(n, sa, cx, incx);

            multiplexer(cx.toArr(), cxExpected)(approximitly);
          });
        });
      });
      describe('error tests', () => {
        const { csscalError: testData } = fixture;
        each(testData)(({ input: { n, cx: x, incx, sa, }, desc }, key) => {

          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

            const call = () => csscal(n, sa, cx, incx);

            expect(call).to.throw();
          });
        });
      });
    });

  });

  describe('cswap', () => {
    describe('data tests', () => {
      const { cswap: testData } = fixture;

      each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, output: { cx: cxOut, cy: cyOut }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const cxExpected = muxCmplx(cxOut.re, cxOut.im);
          const cyExpected = muxCmplx(cyOut.re, cyOut.im);

          cswap(n, cx, incx, cy, incy);

          multiplexer(cx.toArr(), cxExpected)(approximitly);
          multiplexer(cy.toArr(), cyExpected)(approximitly);
        });
      });
    });
    describe('error tests', () => {
      const { cswapError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const call = () => cswap(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });
  });

  describe('icmax', () => {
    describe('data tests', () => {
      const { icamax: testData } = fixture;

      each(testData)(({ input: { n, cx: x, incx }, output: { max }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const result = icamax(n, cx, incx);
          multiplexer(result, max)(approximitly);
        });
      });
    });

    describe('error tests', () => {
      const { icamaxError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, incx }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

          const call = () => icamax(n, cx, incx);

          expect(call).to.throw();
        });
      });
    });
  });

  describe('scasum', () => {
    describe('data tests', () => {
      const { scasum: testData } = fixture;

      each(testData)(({ input: { n, cx: x, incx }, output: { sum }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const result = scasum(n, cx, incx);

          multiplexer(result, sum)(approximitly);
        });
      });
    });
    describe('error tests', () => {
      const { scasumError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, incx }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

          const call = () => scasum(n, cx, incx);

          expect(call).to.throw();
        });
      });
    });
  });
});




