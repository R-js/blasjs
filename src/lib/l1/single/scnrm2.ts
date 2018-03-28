/*
*>  -- jacob bogers, javascript port, 03/2018 (jkfbogers@gmail.com)
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to CLASSQ.
*>     Sven Hammarling, Nag Ltd.
*/

import { FortranArr } from '../../f_func';

const { abs, sqrt } = Math;

/*
*>
*> SCNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    SCNRM2 := sqrt( x**H*x )
*/
export function scnrm2(n: number, x: FortranArr, incx: number): number {

      if (n < 1 || incx < 1) {
            return 0;
      }

      let scale = 0;
      let ssq = 1;

      /*
       The following loop is equivalent to this call to the LAPACK
       auxiliary routine:
       CALL CLASSQ( N, X, INCX, SCALE, SSQ )
      */

      const bx = x.base;

      for (let ix = 1; ix <= 1 + (n - 1) * incx; ix += incx) {
            const i = ix - bx;
            //real
            if (x.r[i] !== 0) {
                  let temp = abs(x.r[i]);
                  if (scale < temp) {
                        //calc once
                        const t1 = scale / temp
                        ssq = 1 + ssq * t1 * t1;
                        scale = temp;
                  }
                  else {
                        const t1 = temp / scale;
                        ssq = ssq + t1 * t1;
                  }
            }
            //img
            if (x.i && x.i[i] !== 0) {
                  let temp = abs(x.i[i]);
                  if (scale < temp) {
                        //calc once
                        const t1 = scale / temp
                        ssq = 1 + ssq * t1 * t1;
                        scale = temp;
                  }
                  else {
                        const t1 = temp / scale;
                        ssq = ssq + t1 * t1;
                  }
            }
      }
      return scale * sqrt(ssq);
}
