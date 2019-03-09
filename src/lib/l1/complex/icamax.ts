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

import { FortranArr } from '../../f_func';

import { errMissingIm } from '../../f_func';
const { abs } = Math;

export function icamax(n: number, cx: FortranArr, incx: number): number {

  if (n <= 0 || incx <= 0) return 0;

  if (n === 1) return 1;

  if (!cx.i) {
    throw new Error(errMissingIm('cx.i'));
  }

  const bx = cx.base;

  let _icamax = 1;
  let smax = abs(cx.r[1 - bx]) + abs(cx.i[1 - bx]);

  let ix = incx +1;
  for (let i = 2; i <= n; i++) {
    let temp = abs(cx.r[ix - bx]) + abs(cx.i[ix - bx]);
    if (temp > smax) {
      smax = temp;
      _icamax = i;
    }
    ix += incx;
  }

  return _icamax;
}
