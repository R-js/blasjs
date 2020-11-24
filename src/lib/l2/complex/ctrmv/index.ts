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

import {
    errMissingIm,
    errWrongArg,
    FortranArr,
    FortranArrEComplex,
    lowerChar,
    Matrix,
    MatrixEComplex,
} from '../../../f_func';

import { normalLower } from './normal-lower';
import { normalUpper } from './normal-upper';
import { transLower } from './trans-lower';
import { transUpper } from './trans-upper';

const { max } = Math;

export function ctrmv(
    uplo: 'u' | 'l',
    trans: 'n' | 't' | 'c',
    diag: 'u' | 'n',
    n: number,
    a: Matrix,
    lda: number,
    x: FortranArr,
    incx: number,
): void {
    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (a.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    // faster then String.toLowerCase()
    const ul = lowerChar(uplo);
    const tr = lowerChar(trans);
    const dg = lowerChar(diag);

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    } else if (!'ntc'.includes(tr)) {
        info = 2;
    } else if (!'un'.includes(dg)) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (lda < max(1, n)) {
        info = 6;
    } else if (incx === 0) {
        info = 8;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('ctpsv', info));
    }

    //*     Quick return if possible.
    if (n === 0) return;

    const noconj = tr === 't';
    const nounit = dg === 'n';

    const kx = incx < 0 ? 1 - (n - 1) * incx : 1;

    let proc: (
        kx: number,
        noconj: boolean,
        nounit: boolean,
        x: FortranArrEComplex,
        incx: number,
        a: MatrixEComplex,
        n: number,
    ) => void;

    if (trans === 'n') {
        if (uplo === 'u') {
            proc = normalUpper;
        } else {
            proc = normalLower;
        }
    } else {
        if (uplo === 'u') {
            proc = transUpper;
        } else {
            proc = transLower;
        }
    }

    return proc(kx, noconj, nounit, <FortranArrEComplex>x, incx, <MatrixEComplex>a, n);
}
