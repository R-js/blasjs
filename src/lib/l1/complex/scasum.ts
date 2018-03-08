/*
  jack dongarra, linpack, 3/11/78.
  jacob bogers, jsport, 03/2018    
*/
import { errMissingIm, FortranArr } from '../../f_func';

const { abs } = Math;

export function scasum(
      n: number,
      cx: FortranArr,
      incx: number
): number {

      let stemp = 0;

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }

      const xb = cx.base;

      if (n <= 0 || incx <= 0) return 0;

      const nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            stemp = stemp + abs(cx.r[i - xb]) + abs(cx.i[i - xb]);
      }
      return stemp;
}
