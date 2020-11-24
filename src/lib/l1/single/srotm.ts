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

import type { FortranArr } from '../../f_func';

export function srotm(N: number, SX: FortranArr, INCX: number, SY: FortranArr, INCY: number, SPARAM: FortranArr): void {
    let SH11, SH12, SH21, SH22, W, Z;
    let I, KX, KY, NSTEPS;

    const bx = SX.base;
    const by = SY.base;
    const bp = SPARAM.base;

    const SFLAG = SPARAM.r[1 - bp];

    if (N <= 0 || SFLAG + 2 === 0) return;
    if (INCX === INCY && INCX > 0) {
        NSTEPS = N * INCX;
        if (SFLAG < 0) {
            SH11 = SPARAM.r[2 - bp];
            SH12 = SPARAM.r[4 - bp];
            SH21 = SPARAM.r[3 - bp];
            SH22 = SPARAM.r[5 - bp];

            for (I = 1; I <= NSTEPS; I += INCX) {
                const W = SX.r[I - bx];
                const Z = SY.r[I - by];

                SX.r[I - bx] = W * SH11 + Z * SH12;
                SY.r[I - by] = W * SH21 + Z * SH22;
            }
        } else if (SFLAG === 0) {
            //checked
            SH12 = SPARAM.r[4 - bp];
            SH21 = SPARAM.r[3 - bp];
            for (I = 1; I <= NSTEPS; I += INCX) {
                W = SX.r[I - bx];
                Z = SY.r[I - by];
                SX.r[I - bx] = W + Z * SH12;
                SY.r[I - by] = W * SH21 + Z;
            }
        } else {
            //checked
            SH11 = SPARAM.r[2 - bp];
            SH22 = SPARAM.r[5 - bp];
            for (I = 1; I <= NSTEPS; I += INCX) {
                W = SX.r[I - bx];
                Z = SY.r[I - by];
                SX.r[I - bx] = W * SH11 + Z;
                SY.r[I - by] = -W + SH22 * Z;
            }
        }
    }
    //INCX !== INCY || INCX <= 0
    else {
        KX = 1;
        KY = 1;
        if (INCX < 0) KX = 1 + (1 - N) * INCX;
        if (INCY < 0) KY = 1 + (1 - N) * INCY;

        if (SFLAG < 0) {
            SH11 = SPARAM.r[2 - bp];
            SH12 = SPARAM.r[4 - bp];
            SH21 = SPARAM.r[3 - bp];
            SH22 = SPARAM.r[5 - bp];
            for (I = 1; I <= N; I++) {
                W = SX.r[KX - bx];
                Z = SY.r[KY - by];
                SX.r[KX - bx] = W * SH11 + Z * SH12;
                SY.r[KY - by] = W * SH21 + Z * SH22;
                KX = KX + INCX;
                KY = KY + INCY;
            }
        } else if (SFLAG === 0) {
            //checked
            SH12 = SPARAM.r[4 - bp];
            SH21 = SPARAM.r[3 - bp];
            for (I = 1; I <= N; I++) {
                W = SX.r[KX - bx];
                Z = SY.r[KY - by];

                SX.r[KX - bx] = W + Z * SH12;
                SY.r[KY - by] = W * SH21 + Z;
                KX += INCX;
                KY += INCY;
            }
        } else {
            //checked
            SH11 = SPARAM.r[2 - bp];
            SH22 = SPARAM.r[5 - bp];
            for (I = 1; I <= N; I++) {
                W = SX.r[KX - bx];
                Z = SY.r[KY - by];
                SX.r[KX - bx] = W * SH11 + Z;
                SY.r[KY - by] = -W + SH22 * Z;
                KX = KX + INCX;
                KY = KY + INCY;
            }
        }
    }
}
