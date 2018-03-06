
import { isamax } from './single/isamax';
import { sasum } from './single/sasum';
import { saxpy } from './single/saxpy';
import { scnrm2 } from './single/scnrm2';
import { scopy } from './single/scopy';
import { sdot } from './single/sdot';
import { sdotdot } from './single/sdsdot';
import { srotmg } from './single/srotmg';

export const level1 = {
    isamax,
    sasum,
    srotmg,
    saxpy,
    scnrm2,
    scopy,
    sdot,
    sdotdot
};
