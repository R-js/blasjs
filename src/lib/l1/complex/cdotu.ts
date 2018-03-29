import { Complex, errMissingIm, FortranArr } from '../../f_func';

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

      if (cx.i === undefined) {
            throw new Error(errMissingIm('cx.i'));
      }
      if (cy.i === undefined) {
            throw new Error(errMissingIm('cy.i'));
      }

      if (n <= 0) return { re: 0, im: 0 };

      let re = 0;
      let im = 0;

      const xb = cx.base;
      const yb = cy.base;

      let ix = 1;
      let iy = 1;
      if (incx < 0) ix = (-n + 1) * incx + 1;
      if (incy < 0) iy = (-n + 1) * incy + 1
      for (let i = 1; i <= n; i++) {
            //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
            re += cx.r[ix - xb] * cy.r[iy - yb] - cx.i[ix - xb] * cy.i[iy - yb];
            im += cx.r[ix - xb] * cy.i[iy - yb] + cx.i[ix - xb] * cy.r[iy - yb];
            ix += incx;
            iy += incy;
      }
      return { re, im };
} 
