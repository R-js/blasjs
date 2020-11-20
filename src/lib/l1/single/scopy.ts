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

export function scopy(
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number): void {

      if (n <= 0) return;
      const xb = sx.base;
      const yb = sy.base;

      if (incx === 1 && incy === 1) {
            const m = n % 7;
            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        sy.r[i - yb] = sx.r[i - xb];
                  }
                  if (n < 7) {
                        return;
                  }
            }
            const mp1 = m + 1;
            for (let i = mp1; i <= n; i++) {
                  let kx = i - xb;
                  let ky = i - yb;
                  // prolly this helped the compiler(fortran) unwind for loops
                  sy.r[ky++] = sx.r[kx++]; //1
                  sy.r[ky++] = sx.r[kx++];
                  sy.r[ky++] = sx.r[kx++];
                  sy.r[ky++] = sx.r[kx++];
                  sy.r[ky++] = sx.r[kx++];
                  sy.r[ky++] = sx.r[kx++];
                  sy.r[ky++] = sx.r[kx++]; //7
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;
            for (let i = 1; i <= n; i++) {
                  let kx = ix - xb;
                  let ky = iy - yb;
                  sy.r[ky] = sx.r[kx];
                  ix += incx;
                  iy += incy;
            }
      }
}
