/*
   -- Jacob Bogers, 03/2018, jkfbogers@mgail
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*  
      Note: the fortran doc states CA is an input var, but it is also set within the
      subroutine
*/

import { Complex } from '../../f_func';

const { sqrt } = Math;

export function crotg(
      ca: Complex, // in
      cb: Complex, // in
      c: { val: number }, // out argument
      s: Complex // out
): void {

      if (ca.re === 0 && ca.im === 0) {
            c.val = 0;
            s.re = 1;
            s.im = 0;
            //NOTE: violates the documentation, but the source code does this
            ca.re = cb.re;
            ca.im = cb.im;
            return;
      }

      const cabs_ca = sqrt(ca.re * ca.re + ca.im * ca.im);
      const cabs_cb = sqrt(cb.re * cb.re + cb.im * cb.im);
      // numerical stability?
      const scale = cabs_ca + cabs_cb;
      const s1 = cabs_ca / scale;
      const s2 = cabs_cb / scale;
      const norm = scale * sqrt(s1 * s1 + s2 * s2);
      //console.log('norm, |ca|,|cb|', norm, cabs_ca, cabs_cb);
      //console.log('scale', scale);

      // c
      const alphaRe = (ca.re / cabs_ca);
      const alphaIm = (ca.im / cabs_ca);
      c.val = cabs_ca / norm;
      //S = ALPHA*CONJG(CB)/NORM
      //(a+ib)*(c-id) = (ac+bd)+i(-ad+bc)
      s.re = (alphaRe * cb.re + alphaIm * cb.im) / norm;
      s.im = (-alphaRe * cb.im + alphaIm * cb.re) / norm;

      //NOTE: violates the documentation, but the source code does this
      ca.re = alphaRe * norm;
      ca.im = alphaIm * norm;
      //
}
