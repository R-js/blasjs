
import { cgbmv } from './complex/cgbmv';
import { cgemv } from './complex/cgemv';
import { cgerc } from './complex/cgerc';
//
import { sgbmv } from './single/sgbmv';
import { sgemv } from './single/sgemv';
import { sger } from './single/sger';
import { ssbmv } from './single/ssbmv';
import { sspmv } from './single/sspmv';
import { sspr } from './single/sspr';
import { sspr2 } from './single/sspr2';
import { ssymv } from './single/ssymv';
import { ssyr } from './single/ssyr';
import { ssyr2 } from './single/ssyr2';
import { stbmv } from './single/stbmv';
import { stbsv } from './single/stbsv';
import { stpmv } from './single/stpmv';
import { stpsv } from './single/stpsv';
import { strmv } from './single/strmv';
import { strsv } from './single/strsv';

export const level2 = {
    //single
    sgbmv, //1
    sgemv,
    sger,
    ssbmv,
    sspmv, //5
    sspr,
    sspr2,
    ssymv,
    ssyr,
    ssyr2, //10
    stbmv,
    stbsv,
    stpmv,
    stpsv,
    strmv, //15
    strsv,
    // double precision
    dgbmv: sgbmv, //1
    dgemv: sgemv,
    dger: sger,
    dsbmv: ssbmv,
    dspmv: sspmv, //5
    dspr: sspr,
    dspr2: sspr2,
    dsymv: ssymv,
    dsyr: ssyr,
    dsyr2: ssyr2, //10
    dtbmv: stbmv,
    dtbsv: stbsv,
    dtpmv: stpmv,
    dtpsv: stpsv,
    dtrmv: strmv, //15
    dtrsv: strsv,
    //complex single
    cgemv,
    cgbmv,
    cgerc
};

