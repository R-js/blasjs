import { assert, expect } from 'chai';
import * as blas from '../src/lib';
import { fixture } from './fixture_l1_complex';
import { approximitly } from './test-helpers';


const {
  util: { arrayrify, numberPrecision, each, multiplexer, fortranArrComplex64, complex },
  level1: { caxpy, ccopy }
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


});
