import { sgemm } from './single/sgemm';
import { ssymm } from './single/ssymm';

export const level3 = {
    //single
    sgemm,
    ssymm
};
