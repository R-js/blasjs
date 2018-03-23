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

};
