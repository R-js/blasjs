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

const { abs } = Math;

export function scasum(
      n: number,
      cx: FortranArr,
      incx: number
): number {

      let stemp = 0;

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }

      const xb = cx.base;

      if (n <= 0 || incx <= 0) return 0;

      const nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            stemp = stemp + abs(cx.r[i - xb]) + abs(cx.i[i - xb]);
      }
      return stemp;
}
