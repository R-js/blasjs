/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, blasjs 03/2018 (jkfbogers@gmail.com)
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
