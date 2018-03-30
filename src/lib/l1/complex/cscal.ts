import { Complex, errMissingIm, FortranArr } from '../../f_func';

/*  jack dongarra, linpack,  3/11/78.
    jacob bogers, 03/2018, JS port
*/

export function cscal(n: number, ca: Complex, cx: FortranArr, incx: number): void {

      if (cx.i === undefined) {
            throw new Error(errMissingIm('cx.i'));
      }

      const bx = cx.base;

      if (n <= 0) return;
      if (incx <= 0) return;

      let nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            const re = ca.re * cx.r[i - cx.base] - ca.im * cx.i[i - cx.base];
            const im = ca.re * cx.i[i - cx.base] + ca.im * cx.r[i - cx.base];
            cx.r[i - bx] = re;
            cx.i[i - bx] = im;
      }
}
