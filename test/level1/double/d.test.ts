import { assert, expect } from 'chai';
import * as blas from '../../../src/lib';
import { approximately } from '../../test-helpers';
import { fixture } from './fixtures';

const {
  util: { arrayrify, numberPrecision, each, multiplexer, fortranArrComplex64, complex, muxCmplx },
  level1: {
    isamax,
    sasum,
    scnrm2,
    scopy,
    sdot,
    sdsdot,
    snrm2,
    srot,
    srotg,
    srotm,
    srotmg,
    sscal,
    sswap,
    saxpy
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 1 single/double precision', function n() {

  describe('isamax', () => {

    const { isamax: testData } = fixture;

    each(testData)(({ input: { n, cx: x, incx }, output: { max }, desc }, key) => {

      it(`[${key}]/[${desc}]`, function t() {

        const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
        const result = isamax(n, cx, incx);
        multiplexer(result, max)(approximately);

      });
    });
  });

  describe('sasum', () => {
    describe('data tests', () => {
      const { sasum: testData } = fixture;

      each(testData)(({ input: { n, cx: x, incx }, output: { sum }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const result = sasum(n, cx, incx);
          multiplexer(result, sum)(approximately);
        });
      });
    });
  });

  describe('saxpy', () => {

    describe('data tests', () => {
      const { saxpy: testData } = fixture;
      each(testData)(({ input: { n, cx: x, cy: y, sa, c, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const res = fortranArrComplex64(muxCmplx(ore, oim))();

          saxpy(n, sa, cx, incx, cy, incy);
          multiplexer(cy.toArr(), res.toArr())(approximately);
        });
      });
    });
  });
  describe('scnrm2', () => {

    describe('data tests', () => {
      const { scnrm2: testData } = fixture;
      each(testData)(({ input: { n, x, incx }, output, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

          const result = scnrm2(n, cx, incx);
          multiplexer(result, output)(approximately);
        });
      });
    });
  });
  describe('ccopy', () => {
    describe('data tests', () => {

      const { scopy: testData } = fixture;
      each(testData)(({ input: { n, x, y, incx, incy }, output: out, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          const expected = fortranArrComplex64(muxCmplx(out.re, out.im))();

          scopy(n, cx, incx, cy, incy);
          multiplexer(cy.toArr(), expected.toArr())(approximately);
        });
      });
    });
  });
  describe('sdot', () => {
    describe('data tests', () => {

      const { sdot: testData } = fixture;
      each(testData)(({ input: { n, sx, sy, incx, incy }, output: out, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(sx.re, sx.im))();
          const cy = fortranArrComplex64(muxCmplx(sy.re, sy.im))();
          const result = sdot(n, cx, incx, cy, incy);

          approximately(result, out);
        });
      });
    });
  });
  describe('sdsdot', () => {
    describe('data tests', () => {

      const { sdsdot: testData } = fixture;
      each(testData)(({ input: { n, sb, sx, sy, incx, incy }, output: out, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(sx.re, sx.im))();
          const cy = fortranArrComplex64(muxCmplx(sy.re, sy.im))();
          const result = sdsdot(n, sb, cx, incx, cy, incy);

          approximately(result, out);
        });
      });
    });
  });
  describe('snrm2', () => {
    describe('data tests', () => {

      const { snrm2: testData } = fixture;
      each(testData)(({ input: { n, x, incx }, output: out, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const result = snrm2(n, cx, incx);

          approximately(result, out);
        });
      });
    });
  });
  describe('srot', () => {
    describe('data tests', () => {

      const { srot: testData } = fixture;
      each(testData)(({ input, input: { n, x, y, incy, incx, c, s }, output: out, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
          const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
          srot(n, cx, incx, cy, incy, c, s);



          multiplexer(out.x, cx.toArr())(approximately);
          multiplexer(out.y, cy.toArr())(approximately);
        });
      });
    });
  });

  describe('srotg', () => {
    describe('data tests', () => {

      const { srotg: testData } = fixture;
      each(testData)(({ input: { sa, sb }, output: expect, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {

          const p = { sa, sb, c: 0, s: 0 };
          srotg(p);
          const { c, s, sa: asa, sb: asb } = p;
          const { c: ec, s: es, sa: esa, sb: esb } = expect;
          multiplexer([c, s, asa, asb], [ec, es, esa, esb])(approximately);
        });
      });
    });
  });
  describe('srotm', () => {
    describe('data tests', () => {

      const { srotm: testData } = fixture;
      each(testData)(({ input: { n, sx, sy, incx, incy, sparam }, output: expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          const x = fortranArrComplex64(muxCmplx(sx))();
          const y = fortranArrComplex64(muxCmplx(sy))();
          const spar = fortranArrComplex64(muxCmplx(sparam))();
          srotm(n, x, incx, y, incy, spar);
          multiplexer(x.toArr(), expect.x)(approximately);
          multiplexer(y.toArr(), expect.y)(approximately);
          /*console.log({
            x: x.toArr(),
            y: y.toArr(),
            ex: expect.x,
            ey: expect.y
          });*/
        });
      });
    });
  });

  describe('srotmg', () => {
    describe('data tests', () => {

      const { srotmg: testData } = fixture;
      each(testData)(({ input, output: expect, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          input.sparam = fortranArrComplex64(muxCmplx(input.sparam, undefined))();

          srotmg(input);
          const { sd1, sd2, sx1, sy1, sparam } = input;
          const { sd1: esd1, sd2: esd2, sx1: esx1, sy1: esy1, sparam: esparam } = expect;
          multiplexer([sd1, sd2, sx1, sy1], [esd1, esd2, esx1, esy1])(approximately);
          multiplexer(sparam.toArr(), expect)(approximately);
        });
      });
    });
  });

  describe('sscal', () => {
    describe('data tests', () => {

      const { sscal: testData } = fixture;
      each(testData)(({ input: { n, sa, sx, incx }, output: expect, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          const x = fortranArrComplex64(muxCmplx(sx))();
          sscal(n, sa, x, incx);
          //console.log({ x: x.toArr(), x2: expect.x });
          multiplexer(x.toArr(), expect.x)(approximately);
        });
      });
    });
  });
  describe('sswap', () => {
    describe('data tests', () => {

      const { sswap: testData } = fixture;
      each(testData)(({ input: { n, sx, incx, sy, incy }, output: expect, desc }, key) => {

        it(`[${key}]/[${desc}]`, function t() {
          //
          const x = fortranArrComplex64(muxCmplx(sx))();
          const y = fortranArrComplex64(muxCmplx(sy))();
          //
          sswap(n, x, incx, y, incy);
          //console.log({ x: x.toArr(), y: expect.x });
          multiplexer(x.toArr(), expect.x)(approximately);
          multiplexer(y.toArr(), expect.y)(approximately);
        });
      });
    });
  });

});

