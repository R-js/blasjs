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

import { FortranArr } from '../../f_func';
const { abs } = Math;

/**
 *
 * @description ISAMAX finds the index of the first element having maximum absolute value.

 * @param n number of elements in input
 * @param sx SX is REAL array, dimension(1 + (N - 1) * abs(INCX))
 * @param incx storage spacing between elements of SX
 * 
 */

export function isamax(
      n: number,
      sx: FortranArr,
      incx: number
): number {
      let _isamax = 0;
      let smax: number;

      if (n < 1 || incx <= 0) return 0;
      if (n === 1) return 1;

      _isamax = 1;

      const a = sx.r;
      const b = sx.base;

      smax = a[1 - b];
      let ix = 1 + incx; // starts at '2' if incx=1
      for (let i = 2; i <= n; i++) {
            const v = a[ix - b];
            if (abs(v) > smax) {
                  _isamax = i;
                  smax = abs(v);
            }
            ix = ix + incx;
      }
      return _isamax;
}
