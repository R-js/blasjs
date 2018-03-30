import { assert, expect } from 'chai';
import * as blas from '../src/lib';
import { fixture } from './fixture_l1_complex';
import { approximitly } from './test-helpers';


const {
  util: { arrayrify, numberPrecision, each, multiplexer, fortranArrComplex64, complex, muxCmplx },
  level1: { caxpy, ccopy, cdotc, cdotu, crotg, cscal, csrot }
} = blas;


const { abs } = Math;
const { isNaN, isFinite } = Number;



describe('blas level 1', function n() {
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

  });
});


