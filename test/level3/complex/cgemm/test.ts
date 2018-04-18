import { assert, expect } from 'chai';
import * as blas from '../../../../src/lib';
import { Matrix } from '../../../../src/lib/f_func';
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
    cgemm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('cgemm', () => {

    describe('data tests', () => {
      const { cgemm: testData } = fixture;
      each(testData)(({
        input: {
          trA,
          trB,
          m, // A(M,K), C(M,N)
          n, // B(K,N), A(M,K)
          k,
          lda, // lda >= M
          ldb, // ldb >= K
          ldc, // ldc >= M
          beta,
          alpha,
          //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
          a,
          b,
          c
        }, expect, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          cgemm(trA, trB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
          //console.log(`after yr:${y.r}, yi:${y.i}`);
          //  console.log(c.r.length);
          /* c.toArr().forEach(
             cpx => console.log(`     +(${cpx.re},\t${cpx.im})`)
           );*/
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(c.toArr(), expect.c)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgemmErrors: errors } = fixture;
      each(errors)(({ input: {
        trA,
        trB,
        m, // A(M,K), C(M,N)
        n, // B(K,N), A(M,K)
        k,
        lda, // lda >= M
        ldb, // ldb >= K
        ldc, // ldc >= M
        beta,
        alpha,
        //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
        a,
        b,
        c
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const bM = fortranMatrixComplex64(b)(1, 1);
          const cM = fortranMatrixComplex64(c)(1, 1);

          //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);


          const call = () => cgemm(trA, trB, m, n, k, alpha, aM, lda, bM, ldb, beta, cM, ldc);
          //call()
          expect(call).to.throw();
        });
      });
    });
  });
});
