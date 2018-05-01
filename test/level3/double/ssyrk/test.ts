import { assert, expect } from 'chai';
import * as blas from '../../../../src/lib';
import { Matrix, real } from '../../../../src/lib/f_func';
import {
  approximately,
  approximatelyWithPrec
} from '../../../test-helpers';
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
  level3: {
    ssyrk
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('ssyrk', () => {

    describe('data tests', () => {
      const { ssyrk: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
          trans,
          uplo,
          k,
          n,
          lda,
          ldc,
          beta,
          alpha,
          a,
          c
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          ssyrk(
            uplo,
            trans,
            n,
            k, alpha, a, lda, beta, c, ldc);

          if (cmd === 'debug') {
            console.log(`c actual:\n ${c.toArr().join(',\n')}\n\n`);
            console.log(`c expected:\n${(<Array<number>>real(expect.c)).join(',\n')}`);
            const r = real(expect.c);
            const d = c.toArr().map((v, i) => v - r[i]);
            console.log(d);

          }
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(c.toArr(), real(expect.c))(approx);
        });
    });
  });

  describe('test errors', () => {
    const { ssyrkErrors: errors } = fixture;
    each(errors)(({ input: {
      cmd,
      uplo,
      trans,
      n,
      k,
      alpha,
      beta,
      lda,
      ldc
    }, desc }, key) => {
      it(`[${key}]/[${desc}]`, function t() {

        const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

        const cM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

        const call = () => ssyrk(
          uplo,
          trans,
          n,
          k, alpha, aM, lda, beta, cM, ldc);
        //call();
        expect(call).to.throw();
      });
    });
  });
});
});
