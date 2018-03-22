import { sgemm } from './single/sgemm';
import { ssymm } from './single/ssymm';
import { ssyr2k } from './single/ssyr2k';

export const level3 = {
    //single
    sgemm,
    ssymm,
    ssyr2k
};
