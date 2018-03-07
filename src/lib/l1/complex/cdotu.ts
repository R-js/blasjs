import { cmult, Complex, complex, FortranArr } from '../../f_func';

/**
 * 
 *  jack dongarra, linpack, 3/11/78.
 *  jacob bogers, 03/2018, jkfbogers@gmail.com
 * 
 */

export function cdotu(
      n: number,
      cx: FortranArr,
      incx: number,
      cy: FortranArr,
      incy: number): Complex {

      const ctemp = complex();
      const xb = cx.base;
      const yb = cy.base;

      if (!cx.i || !cy.i) {
            throw new Error('cx or/and cy are missing complex parts');
      }

      if (n <= 0) return ctemp;
      if (incx === 1 && incy === 1) {
            for (let i = 1; i <= n; i++) {
                  const c = cmult(cx.r[i - xb], cx.i[i - xb], cy.r[i - yb], cy.i[i - yb]);
                  ctemp.re += c.re;
                  ctemp.im += c.im;
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1
            for (let i = 1; i <= n; i++) {
                  const c = cmult(cx.r[ix - xb], cx.i[ix - xb], cy.r[iy - yb], cy.i[iy - yb]);
                  ctemp.re += c.re;
                  ctemp.im += c.im;
                  ix += incx;
                  iy += incy;
            }
      }
      return ctemp;
} 
