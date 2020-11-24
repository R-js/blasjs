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
import type { Complex, FortranArr } from '../../f_func';
import { errMissingIm } from '../../f_func';

export function caxpy(
    n: number,
    ca: Complex,
    cx: FortranArr, // sx has dimension ( 1 + ( N - 1 )*abs( INCX ) )
    incx: number,
    cy: FortranArr, // sy has dimension ( 1 + ( N - 1 )*abs( INCY ) )
    incy: number,
): void {
    if (n <= 0) return;
    const caIsZero = ca.im === 0 && ca.re === 0;
    if (caIsZero) return;

    if (cx.i === undefined) {
        throw new Error(errMissingIm('cx.i'));
    }
    if (cy.i === undefined) {
        throw new Error(errMissingIm('cy.i'));
    }

    const kbx = cy.base;
    const kby = cx.base;

    let ix = 1;
    let iy = 1;
    if (incx < 0) {
        ix = (-n + 1) * incx + 1;
    }
    if (incy < 0) {
        iy = (-n + 1) * incy + 1;
    }

    for (let i = 1; i <= n; i++) {
        //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
        const ra = ca.re * cx.r[ix - kbx] - ca.im * cx.i[ix - kbx];
        const ia = ca.re * cx.i[ix - kbx] + ca.im * cx.r[ix - kbx];

        cy.r[iy - kby] += ra;
        cy.i[iy - kby] += ia;
        //
        ix += incx;
        iy += incy;
    }
    //}
}
