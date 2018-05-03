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
    fortranMatrixComplex32,
    fortranMatrixComplex64,
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
    muxCmplx,
    fortranMatrixComplex64,
    fortranMatrixComplex32
});


/*
WORKS!!
WORKS!!
WORKS!!

function firstSuccess<T, K>(data: T[], fn: (v: T) => Promise<K>): Promise<K> {
    function loop(arr: T[]) {
        const elt = arr.pop();
        if (!elt) {
            return Promise.reject('list exhausted');
        }
        return fn(elt)
            .then(d => d)
            .catch(d => { console.log(d); return loop(arr); });

    }
    //kickoff
    return loop(data);
}

firstSuccess([6, 18, 20, 30, 40],
    function testFunction(n: number) {
        if (n > 5) {
            return Promise.reject(n);
        }
        return Promise.resolve(n);
    })
    .then(r => console.log('result', r))
    .catch(err => console.log('all trials failed:', err));
*/
