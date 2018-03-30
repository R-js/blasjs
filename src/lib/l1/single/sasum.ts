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

      let nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            const k = i - b;
            stemp += abs(a[k]);
      }
      return stemp;

}
