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

import { assert } from 'chai';

import { Complex } from '../src/lib/f_func';


export function isComplex(a): a is Complex {
    return (
        (a !== null) &&
        (typeof a === 'object') &&
        ('re' in a) &&
        ('im' in a) &&
        (typeof a['re'] === 'number' && typeof a['im'] === 'number')
    );
}

export function approximatelyWithPrec(prec: number): (act: number | Complex, exp: number) => void {

    return (act: number | Complex, exp: number) => approximately(act, exp, prec);
}


export function approximately(act: number | Complex, exp: number | Complex, prec = 1E-9) {

    if (isComplex(act)) {
        if (isComplex(exp)) {
            assert.approximately(act.re, exp.re, prec, 'real part numbers are NOT close');
            assert.approximately(act.im, exp.im, prec, 'numbers are NOT close');
            return;
        }
        throw new Error('cannot compare complex type with non complex type');
    }

    switch (true) {
        case isNaN(act):
            assert.isNaN(exp);
            break;
        case isFinite(act):
            //console.log(act, exp, prec);
            assert.approximately(act, <number>exp, prec, 'numbers are NOT close');
            break;
        case !isFinite(act):
            assert.equal(act, exp);
            break;
        default:
            throw new Error(`Icompatible values ${act}, ${exp}`);
    }
}


