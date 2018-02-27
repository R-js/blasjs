/*
*
* Original FORTRAN doc:
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SROTG(SA,SB,C,S)
*
*       .. Scalar Arguments ..
*       REAL C,S,SA,SB
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SROTG construct givens plane rotation.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SA
*> \verbatim
*>          SA is REAL
*> \endverbatim
*>
*> \param[in] SB
*> \verbatim
*>          SB is REAL
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is REAL
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is REAL
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
*> \ingroup single_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*> \endverbatim
*>
*  =====================================================================

    SUBROUTINE SROTG(SA, SB, C, S)

   Reference BLAS level1 routine (version 3.8.0) --
*  Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*/

import { fsign } from './_helpers';

export function srotg(sa: number, sb: number): { sa: number, sb: number, c: number, s: number } {
    let R: number;
    let ROE: number;
    let SCALE: number;
    let Z: number;

    const { abs, sqrt, pow } = Math;

    let c: number;
    let s: number;

    ROE = sb;
    if (abs(sa) > abs(sb)) {
        ROE = sa;
    }
    SCALE = abs(sa) + abs(sb);
    if (SCALE === 0.0) {
        c = 1.0;
        s = 0.0;
        R = 0.0;
        Z = 0.0;
    }
    else {
        R = SCALE * sqrt(pow(sa / SCALE, 2) + pow(sb / SCALE, 2));
        R = fsign(1.0, ROE) * R
        c = sa / R
        s = sb / R
        Z = 1.0
        if (abs(sa) > abs(sb)) Z = s;
        if (abs(sb) >= abs(sa) && c !== 0.0) Z = 1.0 / c;
    }
    sa = R;
    sb = Z; // why are they doing this its an in parameter?
    return { sa, sb, c, s };
}