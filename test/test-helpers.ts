
import { assert, expect } from 'chai';

import { Complex, isComplex } from '../src/lib/f_func';

export function approximatelyWithPrec(prec: number): (act: number | Complex, exp: number) => void {

    return (act: number | Complex, exp: number) => approximately(act, exp, prec);
}


export function approximately(act: number | Complex, exp: number | Complex, prec = 1E-9) {

    if (isComplex(act)) {
        if (isComplex(exp)) {
            assert.approximately(act.re, exp.re, prec, 'real part numbers are NOT close');
            assert.approximately(act.im, exp.im, prec, 'numbers are NOT close');
            return;
        }
        throw new Error('cannot compare complex type with non complex type');
    }

    switch (true) {
        case isNaN(act):
            assert.isNaN(exp);
            break;
        case isFinite(act):
            assert.approximately(act, <number>exp, prec, 'numbers are NOT close');
            break;
        case !isFinite(act):
            assert.equal(act, exp);
            break;
        default:
            throw new Error(`Icompatible values ${act}, ${exp}`);
    }
}


