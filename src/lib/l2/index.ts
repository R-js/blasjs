

import { sgbmv } from './single/sgbmv';
import { sgemv } from './single/sgemv';
import { sger } from './single/sger';
import { ssbmv } from './single/ssbmv';

export const level2 = {
    sgbmv,
    sgemv,
    sger,
    ssbmv
};

