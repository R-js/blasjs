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

import type { FortranArr } from '../../f_func';
const { abs } = Math;

/**
 * 
 * @param n number of elements in input vector(s)
 * @param sx   SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 * @param incx  storage spacing between elements of SX
 */
export function sasum(n: number, sx: FortranArr, incx: number): number {

      let stemp = 0;

      if (n <= 0 || incx <= 0) return stemp;
      const a = sx.r;
      const b = sx.base;

      let nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            const k = i - b;
            stemp += abs(a[k]);
      }
      return stemp;

}
