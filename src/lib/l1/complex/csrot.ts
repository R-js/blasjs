import { complex, errMissingIm, FortranArr } from '../../f_func';

/*
  Univ. of Tennessee
  Univ. of California Berkeley
  Univ. of Colorado Denver
  Jacob Bogers, 03/2018, jkfbogers@gmail.com
*/

export function csrot(
      n: number,
      cx: FortranArr, //size:( 1 + ( N - 1 )*abs( INCX ) )
      incx: number,
      cy: FortranArr, //( 1 + ( N - 1 )*abs( INCY ) )
      incy: number,
      c: number,
      s: number
): void {

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }

      if (!cy.i) {
            throw new Error(errMissingIm('cy'));
      }

      if (n <= 0) return;

      const ct = complex();
      const xb = cx.base;
      const yb = cy.base;

      let ix = 1;
      let iy = 1;
      if (incx < 0) ix = (-n + 1) * incx + 1;
      if (incy < 0) iy = (-n + 1) * incy + 1;

      //place it outside the loop, faster slightly??

      let xr = 0;
      let yr = 0;
      for (let i = 1; i <= n; i++) {
            xr = ix - xb;
            yr = iy - yb;
            ct.re = cx.r[xr] * c + s * cy.r[yr];
            ct.im = cx.i[xr] * c + s * cy.i[yr];

            cy.r[yr] = c * cy.r[yr] - s * cx.r[xr];
            cy.i[yr] = c * cy.i[yr] - s * cx.i[xr];

            cx.r[xr] = ct.re;
            cx.r[xr] = ct.im;

            ix += incx;
            iy += incy;
      }
}

