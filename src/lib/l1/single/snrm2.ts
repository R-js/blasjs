/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
