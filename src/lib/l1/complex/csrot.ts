import { errMissingIm, FortranArr } from '../../f_func';

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

      let ix = 1;
      let iy = 1;
      if (incx < 0) ix = (-n + 1) * incx + 1;
      if (incy < 0) iy = (-n + 1) * incy + 1;

      //place it outside the loop, faster slightly??

      for (let i = 1; i <= n; i++) {
            // CTEMP = C*CX( IX ) + S*CY( IY )
            // C and S are REAL
            const re = cx.r[ix - cx.base] * c + s * cy.r[iy - cy.base];
            const im = cx.i[ix - cx.base] * c + s * cy.i[iy - cy.base];
            //CY( IY ) =- S*CX( IX ) + C*CY( IY ) 
            cy.r[iy - cy.base] = - s * cx.r[ix - cx.base] + c * cy.r[iy - cy.base];
            cy.i[iy - cy.base] = - s * cx.i[ix - cx.base] + c * cy.i[iy - cy.base];
            //CX( IX ) = CTEMP
            cx.r[ix - cx.base] = re;
            cx.i[ix - cx.base] = im;

            ix += incx;
            iy += incy;
      }
}

