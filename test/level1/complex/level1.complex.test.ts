/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//import { expect } from 'chai';
import * as blas from '../../../src/lib';
import { approximately } from '../../test-helpers';
import { fixture } from './fixtures';

const {
    util: { each, multiplexer, fortranArrComplex64, muxCmplx },
    level1: { scasum, cswap, csscal, caxpy, ccopy, cdotc, cdotu, crotg, cscal, csrot, icamax },
} = blas;

describe('blas level 1 complex', function n() {
    describe('caxpy', () => {
        //construct tests for dbeta from fixtures

        //abuse as a for-each loop
        //make sure it is an arrow function for `map` but not an arrow function for `it`
        describe('data tests', () => {
            const { caxpy: testData } = fixture;
            each(testData)(
                ({ input: { n, cx: x, cy: y, ca, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                        const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                        const res = fortranArrComplex64(muxCmplx(ore, oim))();

                        caxpy(n, ca, cx, incx, cy, incy);
                        multiplexer(cy.toArr(), res.toArr())(approximately);
                    });
                },
            );
        });

        describe('error tests', () => {
            const { caxpyErrors: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, ca, incx, incy }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const call = () => caxpy(n, ca, cx, incx, cy, incy);
                    expect(call).toThrow();
                });
            });
        });
    });

    describe('ccopy', () => {
        describe('data tests', () => {
            const { ccopy: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, output: { re: ore, im: oim }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const res = fortranArrComplex64(muxCmplx(ore, oim))();

                    ccopy(n, cx, incx, cy, incy);
                    multiplexer(cy.toArr(), res.toArr())(approximately);
                });
            });
        });

        describe('error tests', () => {
            const { ccopyErrors: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const call = () => ccopy(n, cx, incx, cy, incy);
                    expect(call).toThrow();
                });
            });
        });
    });

    describe('cdotc', () => {
        describe('data tests', () => {
            const { cdotc: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, output: expected, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const answer = cdotc(n, cx, incx, cy, incy);
                    approximately(answer, expected);
                });
            });
        });

        describe('error tests', () => {
            const { cdotcError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const call = () => cdotc(n, cx, incx, cy, incy);
                    expect(call).toThrow();
                });
            });
        });
    });

    describe('cdotu', () => {
        describe('data tests', () => {
            const { cdotu: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, output: expected, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();

                    const answer = cdotu(n, cx, incx, cy, incy);
                    approximately(answer, expected);
                });
            });
        });

        describe('error tests', () => {
            const { cdotuError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const call = () => cdotu(n, cx, incx, cy, incy);
                    expect(call).toThrow();
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
                    approximately(aIn, ca);
                    approximately(bIn, cb);
                    approximately(c.val, cOut);
                    approximately(s, sOut);
                });
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

                    multiplexer(cx.toArr(), result)(approximately);
                });
            });
        });

        describe('error tests', () => {
            const { cscalError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, ca, incx }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const call = () => cscal(n, ca, cx, incx);
                    expect(call).toThrow();
                });
            });
        });
    });

    describe('csrot', () => {
        describe('data tests', () => {
            const { csrot: testData } = fixture;

            each(testData)(
                ({ input: { n, cx: x, incx, cy: y, incy, c, s }, output: { cx: cxOut, cy: cyOut }, desc }, key) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                        const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                        const cxExpected = muxCmplx(cxOut.re, cxOut.im);
                        const cyExpected = muxCmplx(cyOut.re, cyOut.im);

                        csrot(n, cx, incx, cy, incy, c, s);

                        multiplexer(cx.toArr(), cxExpected)(approximately);
                        multiplexer(cy.toArr(), cyExpected)(approximately);
                    });
                },
            );
        });
        describe('error tests', () => {
            const { csrotError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, incx, cy: y, incy, c, s }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();

                    const call = () => csrot(n, cx, incx, cy, incy, c, s);

                    expect(call).toThrow();
                });
            });
        });
    });
    describe('csscal', () => {
        describe('data tests', () => {
            const { csscal: testData } = fixture;

            each(testData)(({ input: { n, cx: x, incx, sa }, output: { cx: cxOut }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cxExpected = muxCmplx(cxOut.re, cxOut.im);

                    csscal(n, sa, cx, incx);

                    multiplexer(cx.toArr(), cxExpected)(approximately);
                });
            });
        });
        describe('error tests', () => {
            const { csscalError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, incx, sa }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

                    const call = () => csscal(n, sa, cx, incx);

                    expect(call).toThrow();
                });
            });
        });
    });

    describe('cswap', () => {
        describe('data tests', () => {
            const { cswap: testData } = fixture;

            each(testData)(
                ({ input: { n, cx: x, cy: y, incx, incy }, output: { cx: cxOut, cy: cyOut }, desc }, key) => {
                    it(`[${key}]/[${desc}]`, function t() {
                        const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                        const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                        const cxExpected = muxCmplx(cxOut.re, cxOut.im);
                        const cyExpected = muxCmplx(cyOut.re, cyOut.im);

                        cswap(n, cx, incx, cy, incy);

                        multiplexer(cx.toArr(), cxExpected)(approximately);
                        multiplexer(cy.toArr(), cyExpected)(approximately);
                    });
                },
            );
        });
        describe('error tests', () => {
            const { cswapError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, cy: y, incx, incy }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();
                    const cy = fortranArrComplex64(muxCmplx(y.re, y.im))();
                    const call = () => cswap(n, cx, incx, cy, incy);
                    expect(call).toThrow();
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
                    multiplexer(result, max)(approximately);
                });
            });
        });

        describe('error tests', () => {
            const { icamaxError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, incx }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

                    const call = () => icamax(n, cx, incx);

                    expect(call).toThrow();
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

                    multiplexer(result, sum)(approximately);
                });
            });
        });
        describe('error tests', () => {
            const { scasumError: testData } = fixture;
            each(testData)(({ input: { n, cx: x, incx }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const cx = fortranArrComplex64(muxCmplx(x.re, x.im))();

                    const call = () => scasum(n, cx, incx);

                    expect(call).toThrow();
                });
            });
        });
    });
});
