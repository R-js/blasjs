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
