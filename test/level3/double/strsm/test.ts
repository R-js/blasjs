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
    strsm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {


  describe('strsm', () => {
    describe('data tests', () => {
      const { strsm: testData } = fixture;
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
          strsm(
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
      const { strsmErrors: errors } = fixture;
      each(errors)(({ input: {
        //  DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        side,
        uplo,
        trans,
        diag,
        m,
        n,
        alpha,
        lda,
        ldb,
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const bM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);
          //SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB
          const call = () => strsm(
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
