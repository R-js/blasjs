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

const { abs, abs: ABS } = Math;
import type { FortranArr } from '../../f_func';

export function srotmg(p: {
    sd1: number;
    sd2: number;
    sx1: number;
    sy1: number; //in
    sparam: FortranArr;
}): void {
    let SD1 = p.sd1;
    let SD2 = p.sd2;
    let SX1 = p.sx1;
    const SY1 = p.sy1;

    //locals
    const GAM = 2 << 11; // 4096
    const GAM2 = GAM * GAM;
    const GAMSQ = 2 << 23; //16777216
    const ONE = 1;
    const RGAMSQ = 1 / (2 << 23); //5.960464477539063e-8

    let SFLAG;
    let SH11;
    let SH12;
    let SH21;
    let SH22;
    let SP1;
    let SP2;
    let SQ1;
    let SQ2;
    let STEMP;
    let SU;
    const TWO = 2;
    const ZERO = 0;

    const pb = p.sparam.base;

    if (SD1 < ZERO) {
        //
        // GO ZERO - H - D - AND - SX1..
        //
        SFLAG = -ONE;
        SH11 = ZERO;
        SH12 = ZERO;
        SH21 = ZERO;
        SH22 = ZERO;

        SD1 = ZERO;
        SD2 = ZERO;
        SX1 = ZERO;
    } else {
        // CASE - SD1 - NONNEGATIVE
        SP2 = SD2 * SY1;
        if (SP2 === ZERO) {
            SFLAG = -TWO;
            p.sparam.r[1 - pb] = SFLAG;
            // cleanup this is not fortran
            p.sd1 = SD1;
            p.sd2 = SD2;
            p.sx1 = SX1;
            p.sy1 = SY1;
            return;
        }
        // REGULAR - CASE..
        SP1 = SD1 * SX1;
        SQ2 = SP2 * SY1;

        SQ1 = SP1 * SX1;

        if (abs(SQ1) > abs(SQ2)) {
            SH21 = -SY1 / SX1;
            SH12 = SP2 / SP1;

            SU = ONE - SH12 * SH21;

            /*  
                In this place in the code:
                    (sd1 >= 0)
                    (sp2 !== 0)
                    |sd1*sx1**2| > |sd2*sy1**2|
                
                therefore:
                    SU = 1 - (-sy1/sx1)*((sd2*sy1)/(sd1*sx1))
                    SU = 1 + (sd2*sy1**2)/(sd1*sx1**2) 
                    SU = 1 + sq2/sq1
                    min(SU) = 1 - max(0,1), 1 is "supremum" (never reached)
                
                Summary:
                    the test "su > 0" is moot,
                    because at this place in the code, it will
                    always be su > 0 that way.
            */

            if (SU > ZERO) {
                SFLAG = ZERO;
                SD1 = SD1 / SU;
                SD2 = SD2 / SU;
                SX1 = SX1 * SU;
            }
        } else {
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
            } else {
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

        // PROCEDURE..SCALE - CHECK
        if (SD1 !== ZERO) {
            while (SD1 <= RGAMSQ || SD1 >= GAMSQ) {
                if (SFLAG === ZERO) {
                    SH11 = ONE;
                    SH22 = ONE;
                    SFLAG = -ONE;
                } else {
                    SH21 = -ONE;
                    SH12 = ONE;
                    SFLAG = -ONE;
                }

                if (SD1 <= RGAMSQ) {
                    SD1 = SD1 * GAM2;
                    SX1 = SX1 / GAM;
                    SH11 = <number>SH11 / GAM;
                    SH12 = <number>SH12 / GAM;
                } else {
                    SD1 = SD1 / GAM2;
                    SX1 = SX1 * GAM;
                    SH11 = <number>SH11 * GAM;
                    SH12 = <number>SH12 * GAM;
                }
            } //while
        } //if

        if (SD2 !== ZERO) {
            while (ABS(SD2) <= RGAMSQ || ABS(SD2) >= GAMSQ) {
                if (SFLAG === ZERO) {
                    SH11 = ONE;
                    SH22 = ONE;
                    SFLAG = -ONE;
                } else {
                    SH21 = -ONE;
                    SH12 = ONE;
                    SFLAG = -ONE;
                }
                if (ABS(SD2) <= RGAMSQ) {
                    SD2 = SD2 * GAM2;
                    SH21 = <number>SH21 / GAM;
                    SH22 = <number>SH22 / GAM;
                } else {
                    SD2 = SD2 / GAM2;
                    SH21 = <number>SH21 * GAM;
                    SH22 = <number>SH22 * GAM;
                }
            } //while
        } //if
    } //if
    // shortcut
    const parr = p.sparam.r;

    if (<number>SFLAG < ZERO) {
        parr[2 - pb] = <number>SH11;
        parr[3 - pb] = <number>SH21;
        parr[4 - pb] = <number>SH12;
        parr[5 - pb] = <number>SH22;
    } else if (SFLAG === ZERO) {
        parr[3 - pb] = <number>SH21;
        parr[4 - pb] = <number>SH12;
    } else {
        parr[2 - pb] = <number>SH11;
        parr[5 - pb] = <number>SH22;
    }
    parr[1 - pb] = <number>SFLAG;
    //wrap it up, => push back
    p.sd1 = SD1;
    p.sd2 = SD2;
    p.sx1 = SX1;
    p.sy1 = SY1;
}
