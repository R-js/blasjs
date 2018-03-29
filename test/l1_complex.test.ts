import { assert, expect } from 'chai';
import * as blas from '../src/lib';
import { fixture } from './fixture_l1_complex';
import { approximitly } from './test-helpers';


const {
  util: { arrayrify, numberPrecision, each, multiplexer, fortranArrComplex64, complex },
  level1: { caxpy, ccopy, cdotc, cdotu, crotg, cscal }
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

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          const res = fortranArrComplex64(multiplexer(ore, oim)((re, im) => ({ re, im })))();

          caxpy(n, ca, cx, incx, cy, incy);
          multiplexer(Array.from(cy.r), Array.from(res.r))(approximitly);
          multiplexer(Array.from(cy.i), Array.from(res.i))(approximitly);
          //const actuals = dbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.asLog);
          //multiplexer(actuals, expectation)(approximitly);
        });
      });
    });

    describe('error tests', () => {
      const { caxpyErrors: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
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

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          const res = fortranArrComplex64(multiplexer(ore, oim)((re, im) => ({ re, im })))();

          ccopy(n, cx, incx, cy, incy);
          multiplexer(Array.from(cy.r), Array.from(res.r))(approximitly);
          multiplexer(Array.from(cy.i), Array.from(res.i))(approximitly);
          //const actuals = dbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.asLog);
          //multiplexer(actuals, expectation)(approximitly);
        });
      });
    });

    describe('error tests', () => {
      const { ccopyErrors: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          const call = () => ccopy(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });

  });

  describe('cdotc', () => {
    describe('data tests', () => {
      const { cdotc: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          //const res = fortranArrComplex64(multiplexer(ore, oim)((re, im) => ({ re, im })))();

          const answer = cdotc(n, cx, incx, cy, incy);
          //console.log(answer, cx, cy);
          multiplexer(answer.re, ore)(approximitly);
          multiplexer(answer.im, oim)(approximitly);
          //const actuals = dbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.asLog);
          //multiplexer(actuals, expectation)(approximitly);
        });
      });
    });

    describe('error tests', () => {
      const { cdotcError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          const call = () => cdotc(n, cx, incx, cy, incy);
          expect(call).to.throw();
        });
      });
    });

  });

  describe('cdotu', () => {
    describe('data tests', () => {
      const { cdotu: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
          //const res = fortranArrComplex64(multiplexer(ore, oim)((re, im) => ({ re, im })))();

          const answer = cdotu(n, cx, incx, cy, incy);
          //console.log(answer, cx, cy);
          multiplexer(answer.re, ore)(approximitly);
          multiplexer(answer.im, oim)(approximitly);
          //const actuals = dbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.asLog);
          //multiplexer(actuals, expectation)(approximitly);
        });
      });
    });


    describe('error tests', () => {
      const { cdotuError: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, ca, c, output: o, incx, incy }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
          const cy = fortranArrComplex64(multiplexer(y.re, y.im)((re, im) => ({ re, im })))();
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

          //const res = fortranArrComplex64(multiplexer(ore, oim)((re, im) => ({ re, im })))();
          const s = { re: 0, im: 0 };
          const c = { val: 0 };

          crotg(aIn, bIn, c, s); // adjusts aIn

          multiplexer([aIn.re, aIn.im], [ca.re, ca.im])(approximitly);
          multiplexer([bIn.re, bIn.im], [cb.re, cb.im])(approximitly);
          approximitly(c.val, cOut);
          multiplexer([s.re, s.im], [sOut.re, sOut.im])(approximitly);
        });
      });
    });

    describe('cscal', () => {
      describe('data tests', () => {
        const { cscal: testData } = fixture;
        each(testData)(({ input: { n, ca, cx: x, incx }, output: { cx: cxOut }, desc }, key) => {
          it(`[${key}]/[${desc}]`, function t() {
            const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
            cscal(n, ca, cx, incx); // adjusts aIn
            multiplexer([Array.from(cx.r), Array.from(cx.i)], [cxOut.re, cxOut.im])(approximitly);
          });
        });
      });

      describe('error tests', () => {
        const { cscalError: testData } = fixture;
        each(testData)(({ input: { n, cx: x, ca, incx }, desc }, key) => {

          it(`[${key}]/[${desc}]`, function t() {

            const cx = fortranArrComplex64(multiplexer(x.re, x.im)((re, im) => ({ re, im })))();
            const call = () => cscal(n, ca, cx, incx);
            expect(call).to.throw();
          });
        });
      });

    });

  });

});
