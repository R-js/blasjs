
import { cgbmv } from './complex/cgbmv';
import { cgemv } from './complex/cgemv';
import { cgerc } from './complex/cgerc';
import { cgeru } from './complex/cgeru';
import { chbmv } from './complex/chbmv';
import { chemv } from './complex/chemv';
import { cher } from './complex/cher';
import { cher2 } from './complex/cher2';
import { chpmv } from './complex/chpmv';
import { chpr } from './complex/chpr';
import { chpr2 } from './complex/chpr2';
import { ctbmv } from './complex/ctbmv';
import { ctbsv } from './complex/ctbsv';
import { ctpmv } from './complex/ctpmv';
import { ctpsv } from './complex/ctpsv';
import { ctrmv } from './complex/ctrmv';
import { ctrsv } from './complex/ctrsv';

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
    cgerc,
    cgeru,
    chbmv,
    chemv,
    cher,
    cher2,
    chpmv,
    chpr,
    chpr2,
    ctbmv,
    ctbsv,
    ctpsv,
    ctrmv,
    ctrsv,
    //complex double
    zgbmv: cgbmv,
    zgemv: cgemv,
    zgerc: cgerc,
    zgeru: cgeru,
    zhbmv: chbmv,
    zhemv: chemv,
    zher: cher,
    zher2: cher2,
    zhpmv: chpmv,
    zhpr: chpr,
    zhpr2: chpr2,
    ztbmv: ctbmv,
    ztbsv: ctbsv,
    ztpmv: ctpmv,
    ztpsv: ctpsv,
    ztrmv: ctrmv,
    ztrsv: ctrsv,
};


