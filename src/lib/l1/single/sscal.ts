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

export function sscal(
      n: number,
      sa: number,
      sx: FortranArr,
      incx: number): void {


      //alias
      const sb = sx.base;


      if (n <= 0 || incx <= 0) return;
      if (sa === 1) return;
      if (incx === 1) {
            /*code for increment equal to 1
            *
            *
            *        clean-up loop*/
            // Munch in batches of 5
            let m = n % 5;
            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        sx.r[i - sb] = sa * sx.r[i - sb];
                  }
                  if (n < 5) return;
            }
            let mp1 = m + 1;
            for (let i = mp1; i <= n; i += 5) {
                  sx.r[i - sb] = sa * sx.r[i - sb];
                  sx.r[i + 1 - sb] *= sa;
                  sx.r[i + 2 - sb] *= sa;
                  sx.r[i + 3 - sb] *= sa;
                  sx.r[i + 4 - sb] *= sa;
            }
      }
      else {

            let NINCX = n * incx;
            for (let i = 1; i <= NINCX; i += incx) {
                  sx.r[i - sb] *= sa;
            }
      }
}
