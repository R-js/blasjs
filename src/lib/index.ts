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



import {
    arrayrify,
    Complex,
    complex,
    each,
    FortranArr,
    fortranArrComplex32,
    fortranArrComplex64,
    fortranMatrixComplex32,
    fortranMatrixComplex64,
    fpArray,
    map,
    Matrix,
    multiplexer,
    muxCmplx,
    numberPrecision
} from './f_func';

const util = Object.freeze({
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

const helper = util;

export { Matrix, Complex, fpArray, FortranArr, util, helper, level1, level2, level3 };




