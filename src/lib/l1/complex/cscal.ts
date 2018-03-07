import { cmult, Complex, FortranArr } from '../../f_func';

/*
    jack dongarra, linpack,  3/11/78.
    jacob bogers, 03/2018, JS port
*/

export function cscal(n: number, ca: Complex, cx: FortranArr, incx: number): void {

      if (n <= 0 || incx <= 0) return;

      const bx = cx.base;

      if (!cx.i) {
            throw new Error('cx doesnt have an imaginary part defined');
      }

      let nincx = n * incx;
      for (let i = 1; i <= nincx; i += incx) {
            const _cx = cmult(ca.re, ca.im, cx.r[i - bx], cx.i[i - bx]);
            cx.r[i - bx] = _cx.re;
            cx.i[i - bx] = _cx.im;
      }

}
