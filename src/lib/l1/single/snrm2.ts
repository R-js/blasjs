/**
   This version written on 25-October-1982.
   Modified on 14-October-1993 to inline the call to SLASSQ.
   Sven Hammarling, Nag Ltd.
   
   Bogers, J.K.F. blasjs 03/2018 (jkfbogers@gmail.com)

*/

import { FortranArr } from '../../f_func';

const { abs, sqrt, pow } = Math;

export function snrm2(n: number, x: FortranArr, incx: number): number {

      if (n < 1 || incx < 1) return 0;
      if (n === 1) return abs(x.r[1 - x.base]);

      let scale = 0;
      let ssq = 1;

      /*
      The following loop is equivalent to this call to the LAPACK
       auxiliary routine:
       CALL SLASSQ( N, X, INCX, SCALE, SSQ )
       */

      for (let ix = 1; ix <= 1 + (n - 1) * incx; ix += incx) {
            const kix = ix - x.base;
            if (x.r[kix] !== 0) {
                  let absxi = abs(x.r[kix]);
                  if (scale < absxi) {
                        ssq = 1 + ssq * pow(scale / absxi, 2);
                        scale = absxi;
                  }
                  else {
                        ssq = ssq + pow(absxi / scale, 2);
                  }
            }
      }
      return scale * sqrt(ssq);
}
