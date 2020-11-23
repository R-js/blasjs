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

const { abs, sqrt } = Math;

export function srotg(p: { sa: number; sb: number; c: number; s: number }): void {
    const sa = p.sa;
    const sb = p.sb;
    const absA = abs(sa);
    const absB = abs(sb);

    let Z = 0;

    const roe = absA > absB ? sa : sb;
    const scale = absA + absB;

    if (scale === 0.0) {
        p.c = 1.0;
        p.s = 0.0;
        p.sa = 0; //R
        p.sb = 0; //Z
        return;
    }
    const v1 = sa / scale;
    const v2 = sb / scale;
    // scale > 0 so R > 0
    let R = scale * sqrt(v1 * v1 + v2 * v2);
    //R = sign(1.0, roe) * R;
    R = roe >= 0 ? R : -R;

    p.c = sa / R;
    p.s = sb / R;
    Z = 1.0;
    if (absA > absB) Z = p.s;
    if (absA <= absB && p.c !== 0.0) Z = 1.0 / p.c;
    // keep as is below
    // end original source code
    // repackage for exit
    p.sa = R;
    p.sb = Z;
}
