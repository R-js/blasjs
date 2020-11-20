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

export function sdsdot(
    n: number,
    sb: number,
    sx: FortranArr,
    incx: number,
    sy: FortranArr,
    incy: number): number {

    let sdot = sb;

    const bx = sx.base;
    const by = sy.base;


    if (n <= 0) {
        return sdot;
    }

    let kx = 1;
    let ky = 1;
    if (incx < 0) {
        kx = 1 + (1 - n) * incx;
    }
    if (incy < 0) {
        ky = 1 + (1 - n) * incy;
    }
    for (let i = 1; i <= n; i++) {
        sdot += sx.r[kx - bx] * sy.r[ky - by];
        kx += incx;
        ky += incy;
    }

    return sdot;
}
