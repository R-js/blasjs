/*
   JS port by Jacob Bogers (jkfbogers@gmail.com)
*/

/**
 * *> \brief \b ISAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ISAMAX(N,SX,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ISAMAX finds the index of the first element having maximum absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] SX
*> \verbatim
*>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of SX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
 * 
 */
import { FortranArr } from '../../f_func';

/**
 *
 * @description ISAMAX finds the index of the first element having maximum absolute value.

 * @param n number of elements in input
 * @param sx SX is REAL array, dimension(1 + (N - 1) * abs(INCX))
 * @param incx storage spacing between elements of SX
 * 
 */
const { abs } = Math;

export function isamax(n: number, sx: FortranArr, incx: number): number {
      let _isamax = 0;
      let smax: number;

      if (n < 1 || incx <= 0) return 0;
      if (n === 1) return 1;

      _isamax = 1;

      const a = sx.arr;
      const b = sx.base;

      smax = a[1 - b];
      let ix = 1 + incx; // starts at '2' if incx=1
      for (let i = 2; i <= n; i++) {
            const v = a[ix - b];
            if (abs(v) > smax) {
                  _isamax = i;
                  smax = abs(v);
            }
            ix = ix + incx;
      }
      return _isamax;
}
