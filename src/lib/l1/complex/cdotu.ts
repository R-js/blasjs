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

import { Complex, errMissingIm, FortranArr } from '../../f_func';

export function cdotu(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number): Complex {
    if (cx.i === undefined) {
        throw new Error(errMissingIm('cx.i'));
    }
    if (cy.i === undefined) {
        throw new Error(errMissingIm('cy.i'));
    }

    if (n <= 0) return { re: 0, im: 0 };

    let re = 0;
    let im = 0;

    const xb = cx.base;
    const yb = cy.base;

    let ix = 1;
    let iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;
    for (let i = 1; i <= n; i++) {
        //(a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
        re += cx.r[ix - xb] * cy.r[iy - yb] - cx.i[ix - xb] * cy.i[iy - yb];
        im += cx.r[ix - xb] * cy.i[iy - yb] + cx.i[ix - xb] * cy.r[iy - yb];
        ix += incx;
        iy += incy;
    }
    return { re, im };
}
