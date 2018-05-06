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

export function csscal(
      n: number,
      sa: number,
      cx: FortranArr, // dimension(1 + (N - 1) * abs(INCX))
      incx: number): void {

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }
      const xb = cx.base;

      if (n <= 0 || incx <= 0) return;
      else {
            let nincx = n * incx;
            for (let i = 1; i <= nincx; i += incx) {
                  const k = i - xb;
                  cx.r[k] = sa !== 0 ? sa * cx.r[k] : 0;
                  cx.i[k] = sa !== 0 ? sa * cx.i[k] : 0;
            }
      }
}
