/*
  jack dongarra, linpack, 3/11/78.
  jacob bogers 03/2018
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
