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
    dgemm
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 3 single/double complex', function n() {

  describe('dgemm', () => {

    describe('data tests', () => {
      const { dgemm: testData } = fixture;
      each(testData)(({
        input: {
          cmd,
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
          dgemm(trA, trB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
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
      const { dgemmErrors: errors } = fixture;
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

          const aM = fortranMatrixComplex64(a)(1, 1).real();
          const bM = fortranMatrixComplex64(b)(1, 1).real();
          const cM = fortranMatrixComplex64(c)(1, 1).real();

          //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);


          const call = () => dgemm(trA, trB, m, n, k, alpha, aM, lda, bM, ldb, beta, cM, ldc);
          //call()
          expect(call).to.throw();
        });
      });
    });
  });
});
