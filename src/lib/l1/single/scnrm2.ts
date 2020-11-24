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
            const temp = abs(x.r[i]);
            if (scale < temp) {
                //calc once
                const t1 = scale / temp;
                ssq = 1 + ssq * t1 * t1;
                scale = temp;
            } else {
                const t1 = temp / scale;
                ssq = ssq + t1 * t1;
            }
        }
        //img
        if (x.i && x.i[i] !== 0) {
            const temp = abs(x.i[i]);
            if (scale < temp) {
                //calc once
                const t1 = scale / temp;
                ssq = 1 + ssq * t1 * t1;
                scale = temp;
            } else {
                const t1 = temp / scale;
                ssq = ssq + t1 * t1;
            }
        }
    }
    return scale * sqrt(ssq);
}
