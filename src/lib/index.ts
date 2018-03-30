import { level1 } from './l1';
import { level2 } from './l2';
import { level3 } from './l3';

export { level1, level2, level3 };

import {
    arrayrify,
    complex,
    each,
    fortranArrComplex32,
    fortranArrComplex64,
    multiplexer,
    muxCmplx,
    numberPrecision

} from './f_func';

export const util = Object.freeze({
    fortranArrComplex32,
    fortranArrComplex64,
    complex,
    each,
    numberPrecision,
    arrayrify,
    multiplexer,
    muxCmplx
});
