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

export function saxpy(
      n: number,
      sa: number,
      sx: FortranArr, // sx has dimension ( 1 + ( N - 1 )*abs( INCX ) )
      incx: number,
      sy: FortranArr, // sy has dimension ( 1 + ( N - 1 )*abs( INCY ) )
      incy: number): void {

      const kbx = sy.base;
      const aY = sy.r;
      const kby = sx.base;
      const aX = sx.r;

      if (n <= 0) return;
      if (sa === 0) return;
      if (incx === 1 && incy === 1) {

            //clean up loop
            const m = n % 4;
            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        aY[i - kby] = aY[i - kby] + sa * aX[i - kbx];
                  }
            }
            if (n < 4) return;
            const mp1 = m + 1;
            for (let i = mp1; i <= n;) {
                  aY[i - kby] = aY[i - kby] + sa * aX[i - kbx];
                  i++;
                  aY[i - kby] = aY[i - kby] + sa * aX[i - kbx];
                  i++;
                  aY[i - kby] = aY[i - kby] + sa * aX[i - kbx];
                  i++
                  aY[i - kby] = aY[i - kby] + sa * aX[i - kbx];
                  i++;
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) {
                  ix = (-n + 1) * incx + 1;
            }
            if (incy < 0) {
                  iy = (-n + 1) * incy + 1;
            }

            for (let i = 1; i <= n; i++) {
                  aY[iy - kby] = aY[iy - kby] + sa * aX[ix - kbx];
                  ix += incx;
                  iy += incy;
            }
      }
}

