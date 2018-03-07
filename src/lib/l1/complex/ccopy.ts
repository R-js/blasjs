import { FortranArr } from '../../f_func';

/*
  jack dongarra, linpack, 3/11/78.
  jacob bogers, jkfbogers@gmail.com, 03/2018
*/

export function ccopy(
      n: number,
      cx: FortranArr,
      incx: number,
      cy: FortranArr,
      incy: number): void {

      const xb = cx.base;
      const yb = cy.base;

      if (!cy.i || !cx.i) {
            throw new Error('cy and cx arguments are missing imaginary parts');
      }

      if (n <= 0) return;
      if (incx === 1 && incy === 1) {
            for (let i = 1; i <= n; i++) {
                  cy.r[i - yb] = cx.r[i - xb];
                  cy.i[i - yb] = cx.i[i - xb];
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;
            for (let i = 1; i <= n; i++) {
                  cy.r[iy - yb] = cx.r[ix - xb];
                  cy.i[iy - yb] = cx.i[ix - xb];
                  ix += incx;
                  iy += incy;
            }
      }

}
