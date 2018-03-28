import { assert } from 'chai';
import * as blas from '../src/lib';
import { fixture } from './fixture_l1_complex';
import { approximitly } from './test-helpers';

const {
  util: { arrayrify, numberPrecision, each, multiplexer }
} = blas;
const { abs } = Math;
const { isNaN, isFinite } = Number;



describe('beta distribution', function n() {
  //

  const precision = numberPrecision(9);

  //


  describe('caxpy', () => {

    //construct tests for dbeta from fixtures
    const { caxpy: testData } = fixture;
    //abuse as a for-each loop
    //make sure it is an arrow function for `map` but not an arrow function for `it`
    each(testData)(({ input: inn, output: expectation, desc }, key) => {

      it(`density test: ${key}/${desc}`, function t() {
        //const actuals = dbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.asLog);
        //multiplexer(actuals, expectation)(approximitly);
      });
    });
  });

  describe('pbeta probability', () => {

    //construct tests for dbeta from fixtures
    const { pbeta: testData } = fixture;
    //abuse as a for-each loop
    //make sure it is an arrow function for `map` but not an arrow function for `it`
    each(testData)(({ input: inn, output: expectation, desc }, key) => {

      it(`propability test: ${key}/${desc}`, function t() {
        //const actuals = pbeta(inn.x, inn.shape1, inn.shape2, inn.ncp, inn.lowerTail, inn.asLog);
        //multiplexer(actuals, expectation)(approximitly);
      });
    });
  });



});
