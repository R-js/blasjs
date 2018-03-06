/*
jack dongarra, linpack, 3/11/78.
jacob bogers, 03/2018
*/

import { FortranArr } from '../../f_func';

export function sscal(
      n: number,
      sa: number,
      sx: FortranArr,
      incx: number): void {


      //alias
      const sb = sx.base;


      if (n <= 0 || incx <= 0) return;
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
                  sx.r[i] = sa * sx.r[i];
                  sx.r[i + 1] *= sa;
                  sx.r[i + 2] *= sa;
                  sx.r[i + 3] *= sa;
                  sx.r[i + 4] *= sa;
            }
      }
      else {

            let NINCX = n * incx;
            for (let i = 1; i <= NINCX; i += incx) {
                  sx.r[i - sb] *= sa;
            }
      }
}
