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

//complex
import { cgemm } from './complex/cgemm';
import { chemm } from './complex/chemm';
import { cher2k } from './complex/cher2k';
import { cherk } from './complex/cherk';
import { csymm } from './complex/csymm';
import { csyr2k } from './complex/csyr2k';
import { csyrk } from './complex/csyrk';
import { ctrmm } from './complex/ctrmm';
import { ctrsm } from './complex/ctrsm';
//single
import { sgemm } from './single/sgemm';
import { ssymm } from './single/ssymm';
import { ssyr2k } from './single/ssyr2k';
import { ssyrk } from './single/ssyrk';
import { strmm } from './single/strmm';
import { strsm } from './single/strsm';
export const level3 = {
    //single
    sgemm,
    ssymm,
    ssyr2k,
    ssyrk,
    strmm,
    strsm,
    //double
    dgemm: sgemm,
    dsymm: ssymm,
    dsyr2k: ssyr2k,
    dsyrk: ssyrk,
    dtrmm: strmm,
    dtrsm: strsm,
    //complex
    cgemm,
    chemm,
    cher2k,
    cherk,
    csymm,
    csyr2k,
    csyrk,
    ctrmm,
    ctrsm,
    //double complex
    zgemm: cgemm,
    zhemm: chemm,
    zher2k: cher2k,
    zherk: cherk,
    zsymm: csymm,
    zsyr2k: csyr2k,
    zsyrk: csyrk,
    ztrmm: ctrmm,
    ztrsm: ctrsm
};

//import './doc';
