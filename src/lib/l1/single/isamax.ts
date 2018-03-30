/*
   jack dongarra, linpack, 3/11/78.
   jacob bogers, javascript 03/2018 (jkfbogers@gmail.com)
*/

import { FortranArr } from '../../f_func';
const { abs } = Math;

/**
 *
 * @description ISAMAX finds the index of the first element having maximum absolute value.

 * @param n number of elements in input
 * @param sx SX is REAL array, dimension(1 + (N - 1) * abs(INCX))
 * @param incx storage spacing between elements of SX
 * 
 */

export function isamax(
      n: number,
      sx: FortranArr,
      incx: number
): number {
      let _isamax = 0;
      let smax: number;

      if (n < 1 || incx <= 0) return 0;
      if (n === 1) return 1;

      _isamax = 1;

      const a = sx.r;
      const b = sx.base;

      smax = a[1 - b];
      let ix = 1 + incx; // starts at '2' if incx=1
      for (let i = 2; i <= n; i++) {
            const v = a[ix - b];
            if (abs(v) > smax) {
                  _isamax = ix;
                  smax = abs(v);
            }
            ix = ix + incx;
      }
      return _isamax;
}
