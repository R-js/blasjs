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

import { level1 } from './l1';
import { level2 } from './l2';
import { level3 } from './l3';

export { level1, level2, level3 };

import {
    arrayrify,
    complex,
    each,
    fortranArrComplex32,
    fortranArrComplex64,
    fortranMatrixComplex32,
    fortranMatrixComplex64,
    map,
    multiplexer,
    muxCmplx,
    numberPrecision
} from './f_func';

export const util = Object.freeze({
    arrayrify,
    complex,
    each,
    fortranArrComplex32,
    fortranArrComplex64,
    fortranMatrixComplex32,
    fortranMatrixComplex64,
    map,
    multiplexer,
    muxCmplx,
    numberPrecision
});

/*

WORKS!!
WORKS!!
WORKS!!

function firstSuccess<T, K>(data: T[], fn: (v: T) => Promise<K>): Promise<K> {
    function loop(arr: T[]) {
        const elt = arr.pop();
        if (!elt) {
            return Promise.reject('list exhausted');
        }
        return fn(elt)
            .then(d => d)
            .catch(d => { console.log(d); return loop(arr); });

    }
    //kickoff
    return loop(data);
}

firstSuccess([6, 18, 20, 30, 40],
    function testFunction(n: number) {
        if (n > 5) {
            return Promise.reject(n);
        }
        return Promise.resolve(n);
    })
    .then(r => console.log('result', r))
    .catch(err => console.log('all trials failed:', err));

*/
