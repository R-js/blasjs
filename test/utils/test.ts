import { assert, expect } from 'chai';

import {
    arrayrify,
    coerceToArray,
    complex,
    each,
    errorMsg,
    fortranArrComplex32,
    fortranMatrixComplex32,
    imaginary,
    lowerChar,
    map,
    Matrix,
    mimicFMatrix,
    mul_cxc,
    mul_rxc,
    multiplexer,
    muxCmplx,
    numberPrecision,
    possibleScalar,
    real
} from '../../src/lib/f_func';

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
    });
    describe('coerceToArray', () => {
        describe('data tests', () => {
            const { coerceToArray: testData } = fixture;
            each(testData)(({
                input: {
                    i
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {

                    const ans = coerceToArray(i);
                    expect(ans).to.be.deep.equal(should.o);
                });
            });
        });
        describe('test errors', () => {
            const { coerceToArrayErr: errors } = fixture;
            each(errors)(({ input: {
                o
            }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const call = () => coerceToArray(o);
                    expect(call).to.throw();
                });
            });
        });
    });
    describe('arrayrify', () => {
        describe('data tests', () => {
            const { arrayrify: testData } = fixture;
            each(testData)(({
                input: {
                    fn,
                    data
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {

                    const arrFn = arrayrify(fn);
                    const ans = arrFn(data);
                    expect(ans).to.be.deep.equal(should.o);
                });
            });
        });
    });
    describe('map', () => {
        describe('data tests', () => {
            const { map: testData } = fixture;
            each(testData)(({
                input: {
                    fn,
                    data
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {

                    const ans = map(data)(fn);
                    expect(ans).to.be.deep.equal(should.o);
                });
            });
        });
    });
    describe('numberPrecision', () => {
        describe('data tests', () => {
            const { numberPrecision: testData } = fixture;
            each(testData)(({
                input: {
                    prec,
                    data
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const fn = numberPrecision(prec);
                    const ans = fn(data);
                    expect(ans).to.be.deep.equal(should.o);
                });
            });
        });
    });
    describe('lowerChar', () => {
        describe('data tests', () => {
            const { lowerChar: testData } = fixture;
            each(testData)(({
                input: {
                    i
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const ans = i ? map(i)(lowerChar).join('') : lowerChar(i);
                    expect(ans).to.be.equal(should.o);
                });
            });
        });
    });
    describe('fortranArrComplex32', () => {
        describe('data tests', () => {
            const { fortranArrComplex32: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                    index,
                    re,
                    im
                }, expect: should, desc
            }, key) => {
                //
                if (/^case[0-9]{1,}$/.test(`${key}`)) {
                    it(`[${key}]/[${desc}]`, function t() {
                        const f32 = fortranArrComplex32(data)();
                        const typeName = f32.r[Symbol.toStringTag];
                        expect(f32.toArr()).to.be.deep.equal(should.o);
                        expect(typeName).to.be.equal(should.type)
                    });
                    return;
                }
                //
                if (/^case[0-9]\/s\(\)/.test(`${key}`)) {
                    it(`[${key}]/[${desc}]`, function t() {
                        const f32 = fortranArrComplex32(data)();
                        const ans = f32.s(index)(re, im)
                        console.log(f32.r, f32.i);
                        const exp = should.o;
                        console.log(ans, exp);

                        expect(ans).to.be.deep.equal(exp);
                    });
                    return;
                }
                //
            });
        });
        describe('error tests', () => {
            const { fortranArrComplex32Err: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                    index,
                    re,
                    im
                }, desc
            }, key) => {
                if (/^case[0-9]\/s\(\)/.test(`${key}`)) {
                    it(`[${key}]/[${desc}]`, function t() {
                        const f32 = fortranArrComplex32(data)();
                        const call = () => f32.s(index)(re, im)
                        expect(call).to.throw();
                    });
                    return;
                }
            });
        });
    });
    describe('fortranMatrixComplex32', () => {
        describe('data tests', () => {
            const { fortranMatrixComplex32: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                    dim1,
                    dim2
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    //
                    const fm32 = fortranMatrixComplex32(data)(dim1, dim2);
                    const _5: (number) => number = numberPrecision(5) as any;

                    const typeNameReal = fm32.r[Symbol.toStringTag];
                    const typeNameImagin = fm32.i ? fm32.i[Symbol.toStringTag] : undefined;
                    //
                    //console.log(_5(fm32.toArr()));
                    if (fm32.i) {
                        expect(typeNameImagin).to.be.equal(should.type);
                    }
                    expect(typeNameReal).to.be.equal(should.type);

                    expect(_5(fm32.toArr())).to.be.deep.equal(should.m32);
                });
            });
        });
    });
    describe('imaginary', () => {
        describe('data tests', () => {
            const { imaginary: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    //
                    const im = imaginary(data);
                    //console.log(im);
                    const _5: (number) => number = numberPrecision(5) as any;
                    expect(_5(im)).to.be.deep.equal(should.o);
                });
            });
        });
    });
    describe('real', () => {
        describe('data tests', () => {
            const { real: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                }, expect: should, desc
            }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const re = real(data);
                    const _5: (number) => number = numberPrecision(5) as any;
                    expect(_5(re)).to.be.deep.equal(should.o);
                });
            });
        });
    });
    describe('mimicMatrix', () => {
        describe('data tests', () => {
            const { mimicMatrix: testData } = fixture;
            each(testData)(({
                input: {
                    data,
                    nrCol,
                    nrRow
                }, expect: should, desc
            }, key) => {
                if (key === 'case0' || key === 'case1') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const fortArr = fortranArrComplex32(data)();
                        const imCopy = mimicFMatrix(fortArr.r, fortArr.i)(nrRow, nrCol).imaginary();
                        const _5: (number) => number = numberPrecision(5) as any;
                        expect(_5(imCopy.toArr())).to.be.deep.equal(_5(should.o));
                    });
                }
                if (key === 'setLower0') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const fortArr = fortranArrComplex32(data)();
                        const imCopy = mimicFMatrix(fortArr.r, fortArr.i)(nrRow, nrCol).real().setLower();
                        const _6: (number) => number = numberPrecision(6) as any;
                        expect(_6(imCopy.toArr())).to.be.deep.equal(_6(should.o));
                    });
                }
                if (key === 'setUpper0') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const fortArr = fortranArrComplex32(data)();
                        const imCopy = mimicFMatrix(fortArr.r, fortArr.i)(nrRow, nrCol).real().setUpper();
                        const _6: (number) => number = numberPrecision(6) as any;
                        expect(_6(imCopy.toArr())).to.be.deep.equal(_6(should.o));
                    });
                }
                if (key === 'upperBand0') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const imCopy = data.upperBand().toArr();
                        const _6: (number) => number = numberPrecision(6) as any;
                        expect(_6(imCopy)).to.be.deep.equal(_6(muxCmplx(should.o.re, should.o.im)));
                    });
                }
                if (key === 'lowerBand0') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const imCopy = data.lowerBand().toArr();
                        const _6: (number) => number = numberPrecision(6) as any;
                        expect(_6(imCopy)).to.be.deep.equal(_6(should.o));
                    });
                }
                if (key === 'packedUpperBand0') {
                    it(`[${key}]/[${desc}]`, function t() {
                        const imCopy = data.packedUpper().toArr();
                        //console.log(imCopy);
                        const _6: (number) => number = numberPrecision(6) as any;
                        expect(_6(imCopy)).to.be.deep.equal(_6(should.o));
                    });
                }
            });
        });
        describe('test errors', () => {
            const { mimicMatrixErr: errors } = fixture;
            each(errors)(({ input: {
                lda,
                nrCols,
                rowBase,
                colBase
            }, desc }, key) => {
                it(`[${key}]/[${desc}]`, function t() {
                    const fortArr = fortranArrComplex32([complex(0.2, 0.4)])();

                    const call = () => mimicFMatrix(fortArr.r, fortArr.i)(lda, nrCols, 'n', rowBase, colBase);
                    //call();
                    expect(call).to.throw();
                });
            });
        });
    });
    describe('possibleScalar', () => {
        describe('data tests', () => {
            it(`[case0/ i = [2]`, function t() {
                const re = possibleScalar([2]);
                expect(re).to.be.equal(2);
            });
        });

    });
    describe('errMsg', () => {
        describe('data tests', () => {
            it(`default message`, function t() {
                const text = errorMsg(99, '42');
                expect(text).to.be.equal('Unkown Error code used! [42]');
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
