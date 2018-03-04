/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, blasjs 03/2018 (jkfbogers@gmail.com)
*/

import { FortranArr } from '../../f_func';

export function scopy(
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number): void {

      if (n < 0) return;
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
            for (let i = mp1; i <= n;) {
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
