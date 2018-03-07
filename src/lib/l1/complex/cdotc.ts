import { Complex, complex, FortranArr } from '../../f_func';

/*
  jack dongarra, linpack,  3/11/78.
  jacob bogers, jkfbogers@gmail.com, 03/2018
*/

export function cdotc(
      n: number,
      cx: FortranArr,
      incx: number,
      cy: FortranArr,
      incy: number): Complex {

      let ctemp = complex();
      const xb = cx.base;
      const yb = cy.base;

      if (!cx.i || !cy.i) {
            const err: string[] = [];
            if (!cx.i) {
                  err.push('cx has no imaginary part');
            }
            if (!cy.i) {
                  err.push('cy has no imaginary part');
            }
            throw new Error(err.join(','));
      }

      if (n <= 0) return ctemp;
      if (incx === 1 && incy === 1) {
            for (let i = 1; i <= n; i++) {
                  //conjugate
                  const Xre = cx.r[i - xb];
                  const Xim = -cx.i[i - xb];

                  const Yre = cy.r[i - yb];
                  const Yim = cy.i[i - yb];

                  //multiply
                  //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                  ctemp.re += Xre * Yre - Xim * Yim;
                  ctemp.im += Xre * Yim + Xim * Yre;
            }
            return ctemp;
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;

            for (let i = 1; i <= n; i++) {
                  //conjugate
                  const Xre = cx.r[ix - xb];
                  const Xim = -cx.i[ix - xb];

                  const Yre = cy.r[iy - yb];
                  const Yim = cy.i[iy - yb];

                  //multiply
                  //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                  ctemp.re += Xre * Yre - Xim * Yim;
                  ctemp.im += Xre * Yim + Xim * Yre;
                  ix += incx;
                  iy += incy;
            }
      }
      return ctemp;
}
