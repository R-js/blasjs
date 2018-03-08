import { complex, errMissingIm, FortranArr } from '../../f_func';

/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, 03/2018, jkfbogers@gmail.com
*/

export function cswap(
      n: number,
      cx: FortranArr, //  dimension(1 + (N - 1) * abs(INCX))
      incx: number,
      cy: FortranArr, // dimension(1 + (N - 1) * abs(INCY))
      incy: number
): void {

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }

      if (!cy.i) {
            throw new Error(errMissingIm('cy'));
      }

      if (n <= 0) return;

      const bx = cx.base;
      const by = cx.base;

      const ctemp = complex();

      let ix = 1;
      let iy = 1;

      if (incx <= 0) ix = (-n + 1) * incx + 1;
      if (incy <= 0) iy = (-n + 1) * incy + 1;

      for (let i = 1; i <= n; i++) {
            const kx = ix - bx;
            const ky = iy - by;
            ctemp.re = cx.r[kx];
            ctemp.im = cx.i[kx];

            cx.r[kx] = cy.r[ky];
            cx.i[kx] = cy.i[ky];

            cy.r[ky] = ctemp.re;
            cy.i[ky] = ctemp.im;

            ix += incx;
            iy += incy;
      }


}
