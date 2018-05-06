/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import { errMissingIm, FortranArr } from '../../f_func';

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
      const by = cy.base;

      let ix = 1;
      let iy = 1;

      if (incx <= 0) ix = (-n + 1) * incx + 1;
      if (incy <= 0) iy = (-n + 1) * incy + 1;

      // case incx = -1 incy=1
      //  x[n] = y[1]
      //  y[1] = x[n]
      //-
      // x[n-1] = y[2]
      // y[2]= x[n-1]
      for (let i = 1; i <= n; i++) {
            const kx = ix - bx;
            const ky = iy - by;

            const re = cx.r[kx];
            const im = cx.i[kx];

            cx.r[kx] = cy.r[ky];
            cx.i[kx] = cy.i[ky];

            cy.r[ky] = re;
            cy.i[ky] = im;

            ix += incx;
            iy += incy;
      }


}
