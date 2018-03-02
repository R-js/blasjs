/** 
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL FUNCTION SCABS1(Z)
*
*       .. Scalar Arguments ..
*       COMPLEX Z
*  
*  PURPOSE
*  ===========
*  SCABS1 computes |Re(.)| + |Im(.)| of a complex number
*
*
*  Arguments:
*  ==========
*
*  [in] Z
*          Z is COMPLEX
*
*  Authors:
*  ========
*
*  Univ. of Tennessee
*  Univ. of California Berkeley
*  Univ. of Colorado Denver
*  NAG Ltd.
*
*  November 2017
*
*  LEVEL 1
*
*  REAL FUNCTION SCABS1(Z)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*/

const { abs } = Math;

export const scabs1 = (re: number, im: number) => abs(re) + abs(im);
