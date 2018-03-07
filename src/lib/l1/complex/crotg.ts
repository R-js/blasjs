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

import { cabs, cmult, complex, Complex } from '../../f_func';

const { sqrt, pow } = Math;

export function crotg(
      p: {
            ca: Complex, // in
            cb: Complex, // in
            c: number, // out argument
            s: Complex // out
      }): void {

      if (p.ca.re === 0 && p.ca.im === 0) {
            p.c = 0;
            p.s = complex(1);
            p.ca.re = p.cb.re;
            p.ca.im = p.cb.im;
      }
      else {
            const cabs_ca = cabs(p.ca.re, p.ca.im);
            const cabs_cb = cabs(p.cb.re, p.cb.im);
            let scale = cabs_ca + cabs_cb;
            let norm = scale * sqrt(
                  pow(cabs_ca / scale, 2) +
                  pow(cabs_cb / scale, 2)
            );
            // c
            let alpha = complex(p.ca.re / cabs_ca, p.ca.im / cabs_ca);
            p.c = cabs_ca / norm;

            // s
            let _s = cmult(alpha.re, alpha.im, p.cb.re, -p.cb.im);
            _s.re /= norm;
            _s.im /= norm;
            p.s.re = _s.re;
            p.s.im = _s.im;

            // ca, strange because CA (accourding in comments)
            //          is an input var not an output var
            p.ca.re = alpha.re * norm;
            p.ca.im = alpha.im * norm;
            //
      }
}
