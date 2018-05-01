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
    strmm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('strmm', () => {

    describe('data tests', () => {
      const { strmm: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
          side,
          uplo,
          trans,
          diag,
          m,
          n,
          alpha,
          a,
          lda,
          b,
          ldb
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          // DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
          strmm(
            side,
            uplo,
            trans,
            diag,
            m,
            n,
            alpha,
            a,
            lda,
            b,
            ldb);

          if (cmd === 'debug') {
            console.log(`b actual:\n ${b.toArr().join(',\n')}\n\n`);
            console.log(`b expected:\n${(<Array<number>>real(expect.b)).join(',\n')}`);
            const r = real(expect.b);
            const d = b.toArr().map((v, i) => v - r[i]);
            console.log(d);

          }
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(b.toArr(), real(expect.b))(approx);
        });
      });
    });

    describe('test errors', () => {
      const { strmmErrors: errors } = fixture;
      each(errors)(({ input: {
        side,
        uplo,
        trans,
        diag,
        m,
        n,
        alpha,
        a,
        lda,
        b,
        ldb
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const bM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const call = () => strmm(
            side,
            uplo,
            trans,
            diag,
            m,
            n,
            alpha,
            aM,
            lda,
            bM,
            ldb);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
