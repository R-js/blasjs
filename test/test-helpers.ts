
import { assert, expect } from 'chai';

import { Complex, isComplex } from '../src/lib/f_func';

export function approximitly(act: number | Complex, exp: number | Complex) {

    if (isComplex(act)) {
        if (isComplex(exp)) {
            assert.approximately(act.re, exp.re, 1e-9, 'real part numbers are NOT close');
            assert.approximately(act.im, exp.im, 1e-9, 'numbers are NOT close');
            return;
        }
        throw new Error('cannot compare complex type with non complex type');
    }

    switch (true) {
        case isNaN(act):
            assert.isNaN(exp);
            break;
        case isFinite(act):
            assert.approximately(act, <number>exp, 1e-9, 'numbers are NOT close');
            break;
        case !isFinite(act):
            assert.equal(act, exp);
            break;
        default:
            throw new Error(`Icompatible values ${act}, ${exp}`);
    }
}


