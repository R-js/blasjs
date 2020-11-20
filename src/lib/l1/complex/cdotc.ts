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
import type { Complex, FortranArr } from '../../f_func';
import { errMissingIm } from '../../f_func';

export function cdotc(
      n: number,
      cx: FortranArr,
      incx: number,
      cy: FortranArr,
      incy: number): Complex {

      let tempRe = 0;
      let tempIm = 0;
      const xb = cx.base;
      const yb = cy.base;

      if (cx.i === undefined) {
            throw new Error(errMissingIm('cx.i'));
      }
      if (cy.i === undefined) {
            throw new Error(errMissingIm('cy.i'));
      }

      if (n <= 0) return { re: tempRe, im: tempIm };

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
            tempRe += Xre * Yre - Xim * Yim;
            tempIm += Xre * Yim + Xim * Yre;
            ix += incx;
            iy += incy;
      }

      return { re: tempRe, im: tempIm };
}
