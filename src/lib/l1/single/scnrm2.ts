/*
  Reference BLAS level1 routine (version 3.8.0) --
  Reference BLAS is a software package provided by Univ. of Tennessee,    --
  Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  November 2017
  jacob bogers, javascript 03/2018 (jkfbogers@gmail.com)

  NOTE: the original Fortran code contained COMPLEX*16 X(*) 
  as a decleration of its complex array argument, thats 2x REAL*8 vars
  for real and imaginary parts
*/

import { FortranArr } from '../../f_func';

const { pow, abs, sqrt } = Math;

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
      if (!x.r || !x.i) {
            throw new Error('complex real and/or imaginary parts are missing');
      }

      for (let ix = 1; ix <= 1 + (n - 1) * incx; ix += incx) {
            const k = ix - bx;
            for (let key in { r: 'r', i: 'i' }) {
                  if (x[key][k] !== 0) {
                        let temp = abs(x[key][k]);
                        if (scale < temp) {
                              ssq = 1 + ssq * pow(scale / temp, 2);
                              scale = temp;
                        }
                        else {
                              ssq = ssq + pow(temp / scale, 2);
                        }
                  }
            }
      }
      return scale * sqrt(ssq);
}
