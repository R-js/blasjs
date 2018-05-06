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

import { Complex, errMissingIm, FortranArr } from '../../f_func';

export function cscal(n: number, ca: Complex, cx: FortranArr, incx: number): void {

      if (cx.i === undefined) {
            throw new Error(errMissingIm('cx.i'));
      }

      const bx = cx.base;

      if (n <= 0) return;
      if (incx <= 0) return;

      let nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            const re = ca.re * cx.r[i - cx.base] - ca.im * cx.i[i - cx.base];
            const im = ca.re * cx.i[i - cx.base] + ca.im * cx.r[i - cx.base];
            cx.r[i - bx] = re;
            cx.i[i - bx] = im;
      }
}
