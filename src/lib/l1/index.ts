
//complex
import { caxpy } from './complex/caxpy';
import { ccopy } from './complex/ccopy';
import { cdotc } from './complex/cdotc';
// single
import { isamax } from './single/isamax';
import { sasum } from './single/sasum';
import { saxpy } from './single/saxpy';
import { scnrm2 } from './single/scnrm2';
import { scopy } from './single/scopy';
import { sdot } from './single/sdot';
import { sdsdot } from './single/sdsdot';
import { snrm2 } from './single/snrm2';
import { srot } from './single/srot';
import { srotg } from './single/srotg';
import { srotm } from './single/srotm';
import { srotmg } from './single/srotmg';
import { sscal } from './single/sscal';
import { sswap } from './single/sswap';

export const level1 = {
    //single
    isamax,
    sasum,
    saxpy,
    scnrm2,
    scopy,
    sdot,
    sdsdot,
    snrm2,
    srot,
    srotg,
    srotm,
    srotmg,
    sscal,
    sswap,
    // double
    dasum: sasum,
    daxpy: saxpy,
    dcopy: scopy,
    ddot: sdot,
    dnrm2: snrm2,
    drot: srot,
    drotg: srotg,
    drotm: srotm,
    drotmg: srotmg,
    dscal: sscal,
    dsdot: sdsdot,
    dswap: sswap,
    idamax: isamax,
    //complex
    caxpy,
    ccopy,
    cdotc
};
