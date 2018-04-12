import { assert, expect } from 'chai';
import * as blas from '../../../src/lib';
import { Matrix } from '../../../src/lib/f_func';
import { approximately, approximatelyWithPrec } from '../../test-helpers';
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
  level2: {
    cgbmv,
    cgemv,
    cgerc,
    cgeru,
    chbmv,
    chemv,
    cher,
    cher2,
    chpr2
  }
} = blas;

const { abs } = Math;
const { isNaN, isFinite } = Number;

describe('blas level 2 single/double complex', function n() {

  describe('cgbmv', () => {

    describe('data tests', () => {
      const { cgbmv: testData } = fixture;

      each(testData)(({ input: {
        trans,
        m,
        n,
        kl,
        ku,
        lda,
        incx,
        incy,
        beta,
        alpha,
        //bandmatrix_nxm_ku_kl(n = 6, m = 6, lda = m, kl = 4, ku = 4)
        a,
        x,
        y
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          // console.log('before', { are: a.r, aim: a.i });
          cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
          //console.log(`after yr:${y.r}, yi:${y.i}`);

          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgbmvErrors: errors } = fixture;
      each(errors)(({ input: {
        trans,
        m,
        n,
        kl,
        ku,
        lda,
        incx,
        incy,
        beta,
        alpha,
        a,
        x,
        y
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();
          //console.log(` a:${aM.r},${aM.i}, x=${JSON.stringify(sx)}`);


          const call = () => cgbmv(trans, m, n, kl, ku, alpha, aM, lda, sx, incx, beta, sy, incy);
          expect(call).to.throw();
        });
      });
    });
  });

  describe('cgemv', () => {

    describe('data tests', () => {
      const { cgemv: testData } = fixture;

      each(testData)(({ input: {
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          //console.log(`before yr:${y.r}, yi:${y.i}`);
          cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgemvErrors: errors } = fixture;
      each(errors)(({ input: {
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgemv(trans, m, n, alpha, aM, lda, sx, incx, beta, sy, incy);
          expect(call).to.throw();
        });
      });
    });
  });
  describe('cgerc', () => {

    describe('data tests', () => {
      const { cgerc: testData } = fixture;

      each(testData)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log('a:', a.toArr());
          //console.log(`before y:${JSON.stringify(y.toArr())}`);
          //console.log(`before x:${JSON.stringify(x.toArr())}`);
          cgerc(m, n, alpha, x, incx, y, incy, a, lda);
          //a.toArr().forEach(cplx => {
          //  console.log(`(${cplx.re},${cplx.im})`);
          //});
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgercErrors: errors } = fixture;
      each(errors)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgerc(m, n, alpha, sx, incx, sy, incy, aM, lda);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
  describe('cgeru', () => {

    describe('data tests', () => {
      const { cgeru: testData } = fixture;

      each(testData)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log('a:', a.toArr());
          //console.log(`before y:${JSON.stringify(y.toArr())}`);
          //console.log(`before x:${JSON.stringify(x.toArr())}`);
          cgeru(m, n, alpha, x, incx, y, incy, a, lda);
          //a.toArr().forEach(cplx => {
          //  console.log(`(${cplx.re},${cplx.im})`);
          //});
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cgeruErrors: errors } = fixture;
      each(errors)(({ input: {
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => cgeru(m, n, alpha, sx, incx, sy, incy, aM, lda);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
  describe('chbmv', () => {

    describe('data tests', () => {
      const { chbmv: testData } = fixture;

      each(testData)(({ input: {
        uplo,
        n,
        k,
        alpha,
        beta,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          //console.log(`before y:${JSON.stringify(y.toArr())}`);
          //console.log(`before x:${JSON.stringify(x.toArr())}`);
          chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
          //console.log(`y:${JSON.stringify(y.toArr())}`);
          //});
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { chbmvErrors: errors } = fixture;
      each(errors)(({ input: {
        uplo,
        n,
        k,
        alpha,
        beta,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => chbmv(uplo, n, k, alpha, aM, lda, sx, incx, beta, sy, incy);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
  describe('chemv', () => {
    describe('data tests', () => {

      const { chemv: testData } = fixture;
      each(testData)(({ input: {
        uplo,
        n,
        lda,
        incx,
        incy,
        alpha,
        beta,
        a,
        x,
        y
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {
          /* if (key === 'case2') {
             a.toArr().forEach(cplx => {
               console.log(`     +  (${cplx.re},${cplx.im}),`);
             });
           }*/
          chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);

          /* y.toArr().forEach(cplx => {
             console.log(`+  (${cplx.re},${cplx.im}),`);
           });*/
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(y.toArr(), expect.y)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { chemvErrors: errors } = fixture;
      each(errors)(({ input: {
        uplo,
        n,
        alpha,
        beta,
        a,
        lda,
        x,
        incx,
        y,
        incy
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          const call = () => chemv(uplo, n, alpha, aM, lda, sx, incx, beta, sy, incy);
          //call();

          //TODO: specify the Exception message explicitly
          expect(call).to.throw();
        });
      });
    });
  });

  describe('cher', () => {
    describe('data tests', () => {
      const { cher: testData } = fixture;
      each(testData)(({ input: {
        //ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)
        uplo,
        n,
        alpha,
        x,
        incx,
        a,
        lda
      }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          //console.log('before a=');
          /* a.toArr().forEach(cplx => {
             console.log(`+  (${cplx.re},${cplx.im}),`);
           });*/
          //console.log('x', x.toArr());
          cher(uplo, n, alpha, x, incx, a, lda);

          //console.log('after a=')
          /* a.toArr().forEach(cplx => {
             console.log(`+  (${cplx.re},${cplx.im}),`);
           });*/
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cherErrors: errors } = fixture;
      each(errors)(({ input: {
        uplo,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
      }, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const call = () => cher(uplo, n, alpha, sx, incx, aM, lda);
          expect(call).to.throw();
        });
      });
    });
  });
  describe('chpr2', () => {
    describe('data tests', () => {
      const { chpr2: testData } = fixture;
      each(testData)(({
        input: {
          //CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
          uplo,
          n,
          alpha,
          x,
          incx,
          y,
          incy,
          ap
        }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          //console.log('before a=');
          /* a.toArr().forEach(cplx => {
             console.log(`+  (${cplx.re},${cplx.im}),`);
           });*/
          /*console.log('before ap');
          ap.toArr().forEach(cplx => {
            console.log(`     +  (${cplx.re},${cplx.im}),`);
          });*/
          //console.log('before y', y.toArr());
          //console.log('before x', x.toArr());
          chpr2(uplo, n, alpha, x, incx, y, incy, ap);
          //console.log('after a=')
          /*          ap.toArr().forEach(cplx => {
                      console.log(`+  (${cplx.re},${cplx.im}),`);
                    });*/
          const approx = approximatelyWithPrec(1E-5);
          multiplexer(ap.toArr(), expect.ap)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { chpr2Errors: errors } = fixture;
      each(errors)(({
        input: {
          uplo,
          n,
          alpha,
          ap,
          y,
          x,
          incx,
          incy
        }, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aP = fortranArrComplex64(ap)();
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          //console.log(aP, sx, sy, n, incx, incy, uplo)
          const call = () => chpr2(uplo, n, alpha, sx, incx, sy, incy, aP);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
  describe('cher2', () => {
    describe('data tests', () => {
      const { cher2: testData } = fixture;
      each(testData)(({
        input: {
          //CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
          uplo,
          n,
          alpha,
          incx,
          incy,
          lda,
          x,
          y,
          a
        }, expect, desc }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          cher2(uplo, n, alpha, x, incx, y, incy, a, lda);

          const approx = approximatelyWithPrec(1E-5);
          multiplexer(a.toArr(), expect.a)(approx);
        });
      });
    });

    describe('test errors', () => {
      const { cher2Errors: errors } = fixture;
      each(errors)(({
        input: {
          uplo,
          n,
          alpha,
          a,
          y,
          x,
          lda,
          incx,
          incy
        }, desc
      }, key) => {
        it(`[${key}]/[${desc}]`, function t() {

          const aM = fortranMatrixComplex64(a)(1, 1);
          const sx = fortranArrComplex64(x)();
          const sy = fortranArrComplex64(y)();

          //console.log(aP, sx, sy, n, incx, incy, uplo)
          const call = () => cher2(uplo, n, alpha, sx, incx, sy, incy, aM, lda);
          //call();
          expect(call).to.throw();
        });
      });
    });
  });
});
