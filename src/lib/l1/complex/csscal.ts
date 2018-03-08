import { errMissingIm, FortranArr } from '../../f_func';

/*
 jack dongarra, linpack, 3/11/78.
 jkfbogers@gmail.com, 03/2018,
*/

export function csscal(
      n: number,
      sa: number, // dimension(1 + (N - 1) * abs(INCX))
      cx: FortranArr, // dimension(1 + (N - 1) * abs(INCX))
      incx: number): void {

      if (!cx.i) {
            throw new Error(errMissingIm('cx'));
      }
      const xb = cx.base;

      if (n <= 0 || incx <= 0) return;
      else {
            let nincx = n * incx;
            for (let i = 1; i <= nincx; i += incx) {
                  const k = i - xb;
                  cx.r[k] = sa * cx.r[k];
                  cx.i[k] = sa * cx.r[k];
            }
      }



}
