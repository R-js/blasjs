import { FortranArr } from '../../f_func';

/*
  jack dongarra, linpack, 3/11/78.
  jacob bogers, 02/2018, jkfbogers@gmail.com
*/

import { errMissingIm, scabs1A } from '../../f_func';

export function icamax(n: number, cx: FortranArr, incx: number): number {

  if (n < 1 || incx <= 0) return 0;

  if (n === 1) return 1;

  if (!cx.i) {
    throw new Error(errMissingIm('cx'));
  }

  const bx = cx.base;

  let ix = 1;
  let _icamax = 1;
  let smax = scabs1A(cx.r[1 - bx], cx.i[1 - bx]);

  ix += incx;
  for (let i = 2; i <= n; i++) {
    let temp = scabs1A(cx.r[ix - bx], cx.i[ix - bx]);
    if (temp > smax) {
      smax = temp;
      _icamax = incx;
    }
    ix += incx;
  }

  return _icamax;
}
