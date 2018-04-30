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
    ssyr2k
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('ssyr2k', () => {

    describe('data tests', () => {
      const { ssyr2k: testData } = fixture;
      each(testData)(({
        //DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
        input: {
          cmd,
          uplo,
          trans,
          n,
          k,
          alpha,
          beta,
          lda,
          ldb,
          ldc,
          a,
          b,
          c
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

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
      const { ssyr2kErrors: errors } = fixture;
      each(errors)(({ input: {
        cmd,
        uplo,
        trans,
        n,
        k,
        alpha,
        beta,
        lda,
        ldb,
        ldc
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const bM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);
          const cM = fortranMatrixComplex64({ re: 0, im: 0 })(1, 1);

          const call = () => ssyr2k(uplo, trans, n, k, alpha, aM, lda, bM, ldb, beta, cM, ldc);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
