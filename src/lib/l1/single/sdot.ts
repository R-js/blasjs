/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, blasjs 03/2018 (jkfbogers@gmail.com)
*/

import { FortranArr } from '../../f_func';


export function sdot(
      n: number,
      sx: FortranArr, //dimension ( 1 + ( N - 1 )*abs( INCX ) )
      incx: number, //dimension ( 1 + ( N - 1 )*abs( INCY ) )
      sy: FortranArr,
      incy: number): number {

      let stemp = 0;
      const bx = sx.base;
      const by = sy.base;
      if (n < 0) return 0;
      if (incx < 1 && incy < 1) {
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
