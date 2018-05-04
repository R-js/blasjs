import {
    complex,
    fortranArrComplex32 as arr32,
    fortranArrComplex64 as arr64,
    fortranMatrixComplex64,
    mul_cxc
} from '../../src/lib/f_func';



import {
    diagonal_nxn,
    matrix_mxn,
    vector
} from '../matrices';

const pI = Infinity;
const nI = -Infinity;
const { PI, sin, cos, abs, sqrt } = Math;

const cospi = x => cos(PI * x);
const sinpi = x => sin(PI * x);


export const fixture = {
    mul_cxc: {
        case0: {
            desc: 'multiply to complex types,c=mul_cxc(a,b), a=0.25-0.9i, b=0.1+0.5i',
            input: {
                a: complex(0.25, -.9),
                b: complex(0.1, 0.5)
            },
            expect: {
                c: complex(0.47499998845160007, 3.5000001043081319E-002)
            },
        },
    },
    mul_rxc: {
        case0: {
            desc: 'multiply to complex fragments with complextype, c=mul_rxc(ra,ia, b), ra=0.25, ia=-0.9i, b=0.1+0.5i',
            input: {
                ra: 0.25,
                ia: -.9,
                b: complex(0.1, 0.5)
            },
            expect: {
                c: complex(0.47499998845160007, 3.5000001043081319E-002)
            },
        },
    },
    multiplexer: {
        case0: {
            desc: 'multiplex( null, undefined, "just a string" )',
            input: {
                i1: null,
                i2: undefined,
                i3: 'just a string'
            },
            expect: {
                c: [
                    'null-undefined-j',
                    'null-undefined-u',
                    'null-undefined-s',
                    'null-undefined-t',
                    'null-undefined- ',
                    'null-undefined-a',
                    'null-undefined- ',
                    'null-undefined-s',
                    'null-undefined-t',
                    'null-undefined-r',
                    'null-undefined-i',
                    'null-undefined-n',
                    'null-undefined-g'
                ]
            },
        },
    },
    multiplexerErr: {
        case0: {
            desc: 'error check, passing a function as argument',
            input: {
                i1: () => { }, //dud
                i2: undefined
            }
        },
    },
    someFunctionErr: {
        case0: {
            desc: '',
            input: {
            }
        },
    },
}
