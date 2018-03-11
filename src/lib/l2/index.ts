
import { sgbmv } from './single/sgbmv';
import { sgemv } from './single/sgemv';
import { sger } from './single/sger';
import { ssbmv } from './single/ssbmv';
import { sspmv } from './single/sspmv';
import { sspr } from './single/sspr';
import { ssymv } from './single/ssymv';

export const level2 = {
    sgbmv,
    sgemv,
    sger,
    ssbmv,
    sspmv,
    sspr,
    ssymv
};

