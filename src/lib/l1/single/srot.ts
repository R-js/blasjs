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

export function srot(
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number,
      c: number,
      s: number): void {

      let STEMP: number;
      let I: number;
      let IX: number;
      let IY: number;

      const bx = sx.base;
      const by = sy.base;

      if (n <= 0) return;

      IX = 1;
      IY = 1;

      if (incx < 0) IX = (-n + 1) * incx + 1;
      if (incy < 0) IY = (-n + 1) * incy + 1;
      for (I = 1; I <= n; I++) {
            const kix = IX - bx;
            const kiy = IY - by;
            STEMP = c * sx.r[kix] + s * sy.r[kiy];
            sy.r[kiy] = - s * sx.r[kix] + c * sy.r[kiy];
            sx.r[kix] = STEMP;
            IX += incx;
            IY += incy;
      }

}
