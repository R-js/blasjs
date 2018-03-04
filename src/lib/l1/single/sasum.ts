/*
  jack dongarra, linpack, 3/11/78.
  jacob bogers, javascript port blas, 3/2018
*/


import { FortranArr } from '../../f_func';
const { abs } = Math;

/**
 * 
 * @param n number of elements in input vector(s)
 * @param sx   SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 * @param incx  storage spacing between elements of SX
 */
export function sasum(n: number, sx: FortranArr, incx: number) {

      let stemp = 0;

      if (n <= 0 || incx <= 0) return stemp;
      const a = sx.r;
      const b = sx.base;

      if (incx === 1) {
            let m = n % 6;
            if (m !== 0) {
                  for (let i = 1; i <= m; i++) {
                        const k = i - b;
                        stemp = stemp + abs(a[k]);
                  }
                  if (n < 6) {
                        return stemp;
                  }
            }
            let mp1 = m + 1;
            for (let i = mp1; i <= n; i += 6) {
                  const k = i - b;
                  stemp = stemp + abs(a[k]) + abs(a[k + 1])
                        + abs(a[k + 2]) + abs(a[k + 3]) + abs(a[k + 4]) + abs(a[k + 5]);
            }
      }
      else {
            let nincx = n * incx;
            for (let i = 1; i <= nincx; i += incx) {
                  const k = i - b;
                  stemp += abs(a[k]);
            }
            return stemp;
      }
}
