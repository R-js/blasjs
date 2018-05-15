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


export function sdot(
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number): number {

      let stemp = 0;
      const bx = sx.base;
      const by = sy.base;
      if (n <= 0) return 0;
      if (incx === 1 && incy === 1) {
            //process in chunks of 5
            let m = n % 5;
            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        stemp += sx.r[i - bx] * sy.r[i - by];
                  }
                  if (n < 5) {
                        return stemp;
                  }
            }
            const mp1 = m + 1; // m could be > 0 here
            for (let i = mp1; i <= n; i += 5) {
                  stemp +=
                        sx.r[i - bx] * sy.r[i - by] +
                        sx.r[i - bx + 1] * sy.r[i - by + 1] +
                        sx.r[i - bx + 2] * sy.r[i - by + 2] +
                        sx.r[i - bx + 3] * sy.r[i - by + 3] +
                        sx.r[i - bx + 4] * sy.r[i - by + 4];
            }
      }
      else {
            let ix = 1;
            let iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;
            for (let i = 1; i <= n; i++) {
                  stemp += sx.r[ix - sx.base] * sy.r[iy - sy.base];
                  ix += incx;
                  iy += incy;
            }
      }
      return stemp;
}
