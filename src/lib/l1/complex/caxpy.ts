/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, javascript 03/2018 (jkfbogers@gmail.com)
*/

import { Complex, FortranArr, scabs1 } from '../../f_func';

export function caxpy(
      n: number,
      ca: Complex,
      cx: FortranArr, // sx has dimension ( 1 + ( N - 1 )*abs( INCX ) )
      incx: number,
      cy: FortranArr, // sy has dimension ( 1 + ( N - 1 )*abs( INCY ) )
      incy: number): void {

      const kbx = cy.base;
      const kby = cx.base;

      if (n <= 0) return;
      if (scabs1(ca) === 0) return;
      if (incx === 1 && incy === 1) {
            for (let i = 1; i <= n; i++) {
                  //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
                  let cxIm = cx.i ? cx.i[i - kbx] : 0;
                  let cxRe = cx.r[i - kbx];

                  let ra = (ca.re * cxRe - ca.im * cxIm);
                  let ia = (ca.re * cxIm + ca.im * cxRe);

                  cy.r[i - kby] += ra;
                  if (cy.i) {
                        cy.i[i - kby] += ia;
                  }
            }
      }
      else {
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
                  let cxIm = cx.i ? cx.i[ix - kbx] : 0;
                  let cxRe = cx.r[ix - kbx];

                  let ra = (ca.re * cxRe - ca.im * cxIm);
                  let ia = (ca.re * cxIm + ca.im * cxRe);

                  cy.r[iy - kby] += ra;
                  if (cy.i) {
                        cy.i[iy - kby] += ia;
                  }
                  //
                  ix += incx;
                  iy += incy;
            }
      }
}
