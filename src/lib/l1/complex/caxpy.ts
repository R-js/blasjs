/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, javascript 03/2018 (jkfbogers@gmail.com)

   CY(IY) = CY(IY) + CA*CX(IX)
*/

import { Complex, errMissingIm, FortranArr } from '../../f_func';

export function caxpy(
      n: number,
      ca: Complex,
      cx: FortranArr, // sx has dimension ( 1 + ( N - 1 )*abs( INCX ) )
      incx: number,
      cy: FortranArr, // sy has dimension ( 1 + ( N - 1 )*abs( INCY ) )
      incy: number): void {

      if (n <= 0) return;
      const caIsZero = ca.im === 0 && ca.re === 0;
      if (caIsZero) return;

      if (cx.i === undefined) {
            throw new Error(errMissingIm('cx.i'));
      }
      if (cy.i === undefined) {
            throw new Error(errMissingIm('cy.i'));
      }


      const kbx = cy.base;
      const kby = cx.base;

      let ix = 1;
      let iy = 1;
      if (incx < 0) {
            ix = (-n + 1) * incx + 1;
      }
      if (incy < 0) {
            iy = (-n + 1) * incy + 1;
      }

      for (let i = 1; i <= n; i++) {
            //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
            let ra = (ca.re * cx.r[ix - kbx] - ca.im * cx.i[ix - kbx]);
            let ia = (ca.re * cx.i[ix - kbx] + ca.im * cx.r[ix - kbx]);

            cy.r[iy - kby] += ra;
            cy.i[iy - kby] += ia;
            //
            ix += incx;
            iy += incy;
      }
      //}
}
