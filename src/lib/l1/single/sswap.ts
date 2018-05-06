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

export function sswap(
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number,
): void {

      let stemp = 0;
      //integers

      const xb = sx.base;
      const yb = sy.base;

      if (n <= 0) return;
      if (incx === 1 && incy === 1) {
            // process in batches of 3
            let m = n % 3;

            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        stemp = sx.r[i - xb];
                        sx.r[i - xb] = sy.r[i - yb];
                        sy.r[i - yb] = stemp;
                  }
                  if (n < 3) return;

            }
            let mp1 = m + 1;
            // NOTE: why do we do it like this
            // the fortran code motivation is 
            // simular to the compiler switch for 
            // "loop unrolling".
            // we keep it here as a rememberance
            for (let i = mp1; i <= n;) {
                  stemp = sx.r[i - xb];
                  sx.r[i - xb] = sy.r[i - yb];
                  sy.r[i - yb] = stemp;
                  i++;
                  stemp = sx.r[i - xb];
                  sx.r[i - xb] = sy.r[i - yb];
                  sy.r[i - yb] = stemp;
                  i++;
                  stemp = sx.r[i - xb];
                  sx.r[i - xb] = sy.r[i - yb];
                  sy.r[i - yb] = stemp;
                  i++;
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;
            for (let i = 1; i <= n; i++) {
                  stemp = sx.r[ix - xb];
                  sx.r[ix - xb] = sy.r[iy - yb];
                  sy.r[iy - yb] = stemp;
                  ix += incx;
                  iy += incy;
            }
      }
}
