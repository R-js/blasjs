import { assert, expect } from 'chai';

import { each, Matrix, mul_cxc, mul_rxc, multiplexer } from '../../src/lib/f_func';
import {
    approximately,
    approximatelyWithPrec,
} from '../test-helpers';

import { fixture } from './fixtures';

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas f_func utility functions', function n() {

    describe('mul_cxc', () => {

        describe('data tests', () => {
            const { mul_cxc: testData } = fixture;
            each(testData)(({
                input: {
                    a,
                    b
                }, expect, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const c = mul_cxc(a, b)
                    const approx = approximatelyWithPrec(1E-5);
                    multiplexer(c, expect.c)(approx);
                });
            });
        });
    });
    describe('mul_rxc', () => {

        describe('data tests', () => {
            const { mul_rxc: testData } = fixture;
            each(testData)(({
                input: {
                    ra,
                    ia,
                    b
                }, expect, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const c = mul_rxc(ra, ia, b)
                    const approx = approximatelyWithPrec(1E-5);
                    multiplexer(c, expect.c)(approx);
                });
            });
        });
    });
    describe('multiplexer', () => {
        describe('data tests', () => {
            const { multiplexer: testData } = fixture;
            each(testData)(({
                input: {
                    i1,
                    i2,
                    i3
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const c = multiplexer(i1, i2, i3)((a, b, c) => `${a}-${b}-${c}`);
                    expect(c).to.be.deep.equal(should.c)
                    expect(should.c).to.be.deep.equal(c)
                });
            });
        });
        describe('test errors', () => {
            const { multiplexerErr: errors } = fixture;
            each(errors)(({ input: {
                i1,
                i2
            }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const call = () => multiplexer(i1, i2)((a, b, c) => `${a}-${b}-${c}`);
                    expect(call).to.throw();
                });
            });
        });
        /*
        describe('test errors', () => {
            const { someFunctionErr: errors } = fixture;
            each(errors)(({ input: {
            }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    //call()
                    //expect(call).to.throw();
                });
            });
        });*/
    });
});
