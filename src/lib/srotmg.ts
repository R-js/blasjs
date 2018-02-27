/**> \brief \b SROTMG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
*
*       .. Scalar Arguments ..
*       REAL SD1,SD2,SX1,SY1
*       ..
*       .. Array Arguments ..
*       REAL SPARAM(5)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*>    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*>    SY2)**T.
*>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*>
*>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
*>
*>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
*>    H=(          )    (          )    (          )    (          )
*>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
*>    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
*>    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
*>    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
*>
*>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*>    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] SD1
*> \verbatim
*>          SD1 is REAL
*> \endverbatim
*>
*> \param[in,out] SD2
*> \verbatim
*>          SD2 is REAL
*> \endverbatim
*>
*> \param[in,out] SX1
*> \verbatim
*>          SX1 is REAL
*> \endverbatim
*>
*> \param[in] SY1
*> \verbatim
*>          SY1 is REAL
*> \endverbatim
*>
*> \param[out] SPARAM
*> \verbatim
*>          SPARAM is REAL array, dimension (5)
*>     SPARAM(1)=SFLAG
*>     SPARAM(2)=SH11
*>     SPARAM(3)=SH21
*>     SPARAM(4)=SH12
*>     SPARAM(5)=SH22
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
*  =====================================================================
*/

/* *
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017

*/

const CONSTANTS = Object.freeze({
    '0.5**-24': Math.pow(1 / 2, -24)
});

const { abs, abs: ABS, pow } = Math;

function srotmg(
    sd1: number,
    sd2: number,
    sx1: number,
    sy1: number): {
        SD1: number,
        SD2: number,
        SX1: number,
        SY1: number,
        SPARAM: number[]
    } {

    let SPARAM = Array(5).fill(0);
    let SD1 = sd1, SD2 = sd2, SX1 = sx1, SY1 = sy1;

    //locals
    let GAM = 2 << 11; // 4096
    let GAMSQ = 2 << 23; //16777216
    let ONE = 1;
    let RGAMSQ = CONSTANTS['0.5**-24'];

    let SFLAG, SH11, SH12, SH21, SH22, SP1, SP2, SQ1,
        SQ2, STEMP, SU, TWO = 2, ZERO = 0;


    if (SD1 < ZERO) {
        // GO ZERO - H - D - AND - SX1..
        SFLAG = -ONE;
        SH11 = ZERO;
        SH12 = ZERO;
        SH21 = ZERO;
        SH22 = ZERO;

        SD1 = ZERO;
        SD2 = ZERO;
        SX1 = ZERO;
    }
    else {
        // CASE - SD1 - NONNEGATIVE
        SP2 = SD2 * SY1;
        if (SP2 === ZERO) {
            SFLAG = -TWO;
            SPARAM[0] = SFLAG;
            return { SD1, SD2, SX1, SY1, SPARAM };
        }
        // REGULAR - CASE..
        SP1 = SD1 * SX1;
        SQ2 = SP2 * SY1;
        SQ1 = SP1 * SX1;

        if (abs(SQ1) > abs(SQ2)) {
            SH21 = -SY1 / SX1;
            SH12 = SP2 / SP1;

            SU = ONE - SH12 * SH21;

            if (SU > ZERO) {
                SFLAG = ZERO
                SD1 = SD1 / SU;
                SD2 = SD2 / SU;
                SX1 = SX1 * SU;
            }
        }
        else {
            if (SQ2 < ZERO) {
                // GO ZERO - H - D - AND - SX1..
                SFLAG = -ONE;
                SH11 = ZERO;
                SH12 = ZERO;
                SH21 = ZERO;
                SH22 = ZERO;

                SD1 = ZERO;
                SD2 = ZERO;
                SX1 = ZERO;
            }
            else {
                SFLAG = ONE;
                SH11 = SP1 / SP2;
                SH22 = SX1 / SY1;
                SU = ONE + SH11 * SH22;
                STEMP = SD2 / SU;
                SD2 = SD1 / SU;
                SD1 = STEMP;
                SX1 = SY1 * SU;
            }
        }

        // PROCESURE..SCALE - CHECK
        if (SD1 !== ZERO) {
            while ((SD1 <= RGAMSQ) || (SD1 >= GAMSQ)) {
                if (SFLAG === ZERO) {
                    SH11 = ONE;
                    SH22 = ONE;
                    SFLAG = -ONE;
                }
                else {
                    SH21 = -ONE;
                    SH12 = ONE;
                    SFLAG = -ONE;
                }
                if (SD1 <= RGAMSQ) {
                    SD1 = SD1 * pow(GAM, 2);
                    SX1 = SX1 / GAM;
                    SH11 = SH11 / GAM;
                    SH12 = SH12 / GAM;
                }
                else {
                    SD1 = SD1 / pow(GAM, 2);
                    SX1 = SX1 * GAM;
                    SH11 = SH11 * GAM;
                    SH12 = SH12 * GAM;
                }
            }//while

        }//if

        if (SD2 !== ZERO) {
            while ((ABS(SD2) <= RGAMSQ) || (ABS(SD2) >= GAMSQ)) {
                if (SFLAG === ZERO) {
                    SH11 = ONE;
                    SH22 = ONE;
                    SFLAG = -ONE;
                }
                else {
                    SH21 = -ONE;
                    SH12 = ONE;
                    SFLAG = -ONE;
                }
                if (ABS(SD2) <= RGAMSQ) {
                    SD2 = SD2 * pow(GAM, 2);
                    SH21 = SH21 / GAM;
                    SH22 = SH22 / GAM;
                }
                else {
                    SD2 = SD2 / pow(GAM, 2);
                    SH21 = SH21 * GAM;
                    SH22 = SH22 * GAM;
                }
            }//while
        }//
    }//if
    if (SFLAG < ZERO) {
        SPARAM[1] = SH11;
        SPARAM[2] = SH21;
        SPARAM[3] = SH12;
        SPARAM[4] = SH22;
    }
    else if (SFLAG === ZERO) {
        SPARAM[2] = SH21;
        SPARAM[3] = SH12;
    }
    else {
        SPARAM[1] = SH11;
        SPARAM[4] = SH22;
    }
    SPARAM[0] = SFLAG
    return { SD1, SD2, SX1, SY1, SPARAM };
}




