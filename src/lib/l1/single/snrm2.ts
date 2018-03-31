/**
   This version written on 25-October-1982.
   Modified on 14-October-1993 to inline the call to SLASSQ.
   Sven Hammarling, Nag Ltd.
   
   Bogers, J.K.F. blasjs 03/2018 (jkfbogers@gmail.com)

*/

import { FortranArr } from '../../f_func';

const { abs, sqrt } = Math;

export function snrm2(n: number, x: FortranArr, incx: number): number {

      if (n < 1 || incx < 1) return 0;
      if (n === 1) return abs(x.r[1 - x.base]);

      let scale = 0;
      let ssq = 1;

      /*
      The following loop is equivalent to this call to the LAPACK
       auxiliary routine:
       CALL SLASSQ( N, X, INCX, SCALE=0, SSQ=1 , )

       * DOCUMENTATION FROM SLASSQ
*  =======
*http://www.netlib.org/lapack/explore-3.1.1-html/slassq.f.html
*
*  SLASSQ  returns the values  scl  and  smsq  such that
*  
*     SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.

       */

      for (let ix = 1; ix <= 1 + (n - 1) * incx; ix += incx) {
            const kix = ix - x.base;
            if (x.r[kix] !== 0) {
                  let absxi = abs(x.r[kix]);
                  const ratio = scale / absxi;
                  const ratioP2 = ratio * ratio;
                  if (scale < absxi) {
                        ssq = 1 + ssq * ratioP2;
                        scale = absxi;
                  }
                  else {
                        ssq = ssq + 1 / ratioP2;
                  }
            }
      }
      return scale * sqrt(ssq);
}
