/*
 jack dongarra, linpack, 3/11/78.
 jacob bogers, javascript 03/2018 (jkfbogers@gmail.com)
*/

import { FortranArr } from '../../f_func';

export function saxpy(
      n: number,
      sa: number,
      sx: FortranArr, // sx has dimension ( 1 + ( N - 1 )*abs( INCX ) )
      incx: number,
      sy: FortranArr, // sy has dimension ( 1 + ( N - 1 )*abs( INCY ) )
      incy: number): void {

      const kbx = sy.base;
      const aY = sy.arr;
      const kby = sx.base;
      const aX = sx.arr;

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

