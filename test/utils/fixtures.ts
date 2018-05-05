import {
    complex,
    //fortranArrComplex32 as arr32,
    // fortranArrComplex64 as arr64,
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
    coerceToArray: {
        case0: {
            desc: 'o=9',
            input: {
                i: 9,
            },
            expect: {
                o: [{ key: 0, val: 9 }]
            }
        },
        case1: {
            desc: 'o=[1,2,3,9]',
            input: {
                i: [1, 2, 3, 9],
            },
            expect: {
                o: [
                    { key: 0, val: 1 },
                    { key: 1, val: 2 },
                    { key: 2, val: 3 },
                    { key: 3, val: 9 }
                ]
            }
        },
        case2: {
            desc: 'o="hello"',
            input: {
                i: 'hello'
            },
            expect: {
                o: [
                    { key: 0, val: 'h' },
                    { key: 1, val: 'e' },
                    { key: 2, val: 'l' },
                    { key: 3, val: 'l' },
                    { key: 4, val: 'o' }
                ]
            }
        },
    },
    coerceToArrayErr: {
        case0: {
            desc: 'o=null',
            input: {
                o: null, //dud
            }
        },
        case1: {
            desc: 'o={}',
            input: {
                o: {}, //dud
            }
        },
        case2: {
            desc: 'o=Function',
            input: {
                o: Function, //dud
            }
        },
    },
    arrayrify: {
        case0: {
            desc: 'test arrayrify helper',
            input: {
                fn: Math.sin,
                data: [PI / 6, PI / 3, PI / 2]
            },
            expect: {
                o: [
                    Math.sin(PI / 6),
                    Math.sin(PI / 3),
                    Math.sin(PI / 2)
                ]
            }
        },
        case1: {
            desc: 'test arrayrify helper, scalar input',
            input: {
                fn: Math.sin,
                data: PI / 6
            },
            expect: {
                o: Math.sin(PI / 6)
            }
        },
    },
    map: {
        case0: {
            desc: 'test map helper',
            input: {
                fn: (v, idx) => `${v}:${idx}`,
                data: 'hello'
            },
            expect: {
                o: ['h:0', 'e:1', 'l:2', 'l:3', 'o:4']
            }
        }
    },
    numberPrecision: {
        case0: {
            desc: 'prec=4',
            input: {
                prec: 4,
                data: [0.1234567890, 1.987676E+100, 0.2323423E-18, 2.343837462342E+28, NaN]
            },
            expect: {
                o: [0.1235, 1.988e+100, 2.323e-19, 2.344e+28, NaN]
            }
        },
        case1: {
            desc: 'test default prec=6',
            input: {
                prec: undefined,
                data: [0.1234567890, 1.987676E+100, 0.2323423E-18, 2.343837462342E+28, NaN]
            },
            expect: {
                o: [0.123457, 1.98768e+100, 2.32342e-19, 2.34384e+28, NaN]
            }
        },
        case2: {
            desc: 'test numberPrecision on JS object {}',
            input: {
                prec: undefined,
                data: { a: 0.1234567890, b: 1.987676E+100, c: Math.sin }
            },
            expect: {
                o: { a: 0.123457, b: 1.98768e+100, c: Math.sin }
            }
        },
        case3: {
            desc: 'test numberPrecision on a Function',
            input: {
                prec: undefined,
                data: Math.sin
            },
            expect: {
                o: Math.sin
            }
        }
    },
    lowerChar: {
        case0: {
            desc: 'test lowerChar helper',
            input: {
                i: 'aBcDeFç&@$*µù%><²³à}',
            },
            expect: {
                o: 'abcdefç&@$*µù%><²³à}'
            }
        },
        case1: {
            desc: 'test lowerChar helper',
            input: {
                i: undefined,
            },
            expect: {
                o: ''
            }
        }
    },
    fortranArrComplex32: {
        case0: {
            desc: 'complex number input',
            input: {
                data: complex(1, 0.5),
            },
            expect: {
                type: 'Float32Array',
                o: [{ re: 1, im: 0.5 }],
            }
        },
        case1: {
            desc: 'real number input',
            input: {
                data: [1, 2, 3],
            },
            expect: {
                type: 'Float32Array',
                o: [1, 2, 3]
            },
        },
        case2: {
            desc: 'complex with no imaginary',
            input: {
                data: [{ re: 4 }],
            },
            expect: {
                type: 'Float32Array',
                o: [4]
            },
        },
        case3: {
            desc: 'complex with no real',
            input: {
                data: [{ im: 4 }],
            },
            expect: {
                type: 'Float32Array',
                o: [{ re: undefined, im: 4 }]
            },
        },
        'case0/s()': {
            desc: 'get  array element',
            input: {
                data: [
                    complex(1, 0.5),
                    complex(2, 3),
                ],
                index: 2
            },
            expect: {
                o: complex(2, 3)
            },
        },
        'case1/s()': {
            desc: 'get (real) array element',
            input: {
                data: [1, 2],
                index: 2
            },
            expect: {
                o: 2
            },
        },
        'case2/s()': {
            desc: 'set element of real array',
            input: {
                data: [1, 2],
                index: 2,
                re: 99
            },
            expect: {
                o: { re: 2, im: 0 }
            },
        },
    },
    fortranArrComplex32Err: {
        'case0/s()': {
            desc: 'set complex value in a real array',
            input: {
                data: [1, 2],
                index: 2,
                re: 4,
                im: 2
            }
        },
        'case1/s()': {
            desc: 'index out of bounds, error',
            input: {
                data: [1, 2],
                index: 4
            },
            expect: {
                o: 0
            },
        },

    },
    complex: {
        case0: {
            desc: 'creating complex number with defaults',
            input: {
                data: complex(),
            },
            expect: {
                o: [{ re: 0, im: 0 }],
            }
        },
    },
    fortranMatrixComplex32: {
        case0: {
            desc: 'create 2x2 complex matrix',
            input: {
                data: [
                    complex(1, 0.5),
                    complex(0.5, 0.5),
                    complex(0.8, 0.2),
                    complex(0.3, -0.4)
                ],
                dim1: 2,
                dim2: 2
            },
            expect: {
                type: 'Float32Array',
                m32: [
                    complex(1, 0.5),
                    complex(0.5, 0.5),
                    complex(0.8, 0.2),
                    complex(0.3, -0.4)
                ],
            }
        },
        case1: {
            desc: 'create 2x2 real matrix',
            input: {
                data: [1, 0.5, 0.8, 0.3],
                dim1: 2,
                dim2: 2
            },
            expect: {
                type: 'Float32Array',
                m32: [1, 0.5, 0.8, 0.3]
            }
        },
    },
    imaginary: {
        case0: {
            desc: 'filter imaginary part of complex scalar, array',
            input: {
                data: [
                    complex(1, 0.5),
                    complex(0.5, 0.5),
                    complex(0.8, 0.2),
                    complex(0.3, -0.4),
                    Math.sin,
                    0
                ],
            },
            expect: {
                o: [0.5, 0.5, 0.2, -0.4, Math.sin, 0]
            }
        }
    },
    real: {
        case0: {
            desc: 'filter real part of complex scalar, array',
            input: {
                data: [
                    complex(1, 0.5),
                    complex(0.5, 0.5),
                    complex(0.8, 0.2),
                    complex(0.3, -0.4),
                    Math.sin,
                    5
                ],
            },
            expect: {
                o: [1, 0.5, 0.8, 0.3, Math.sin, 5]
            }
        }
    },
    mimicMatrix: {
        case0: {
            desc: 'return imaginary part as a new real-valued matrix',
            input: {
                data: vector(9).toArr(),
                nrCol: 3,
                nrRow: 3,
            },
            expect: {
                o: Array.from(vector(9).i)
            }
        },
        case1: {
            desc: '(trivial) return non-existant imaginary of a source matrix!',
            input: {
                data: Array.from(vector(9).r),
                nrCol: 3,
                nrRow: 3,
            },
            expect: {
                o: new Array(9).fill(0)
            }
        },
        setLower0: {
            desc: 'setLower without Imaginary part',
            input: {
                data: Array.from(vector(9).r),
                nrCol: 3,
                nrRow: 3,
            },
            expect: {
                o: [-0.6490100777088978,
                    0,
                    0,
                    1.100969102194087,
                    0.14377148075806995,
                    0,
                -0.9120683669483379,
                -1.4375862408299789,
                -0.7970895250719646]
            }
        },
        setUpper0: {
            desc: 'setUpper without Imaginary part',
            input: {
                data: Array.from(vector(9).r),
                nrCol: 3,
                nrRow: 3,
            },
            expect: {
                o: [-0.6490100777088978,
                -0.11916876241803812,
                    0.6641356998941105,
                    0,
                    0.14377148075806995,
                -0.11775359816595128,
                    0,
                    0,
                - 0.7970895250719646
                ]
            }
        },
        upperBand0: {
            desc: 'upperBand (default argument)',
            input: {
                data: matrix_mxn(3, 3)
            },
            expect: {
                o: {
                    re: [
                        0,
                        0,
                        1.2629542848807933,
                        0,
                        1.2724293214294047,
                        0.4146414344564082,
                        -0.9285670347135381,
                        -0.2947204467905602,
                        -0.005767172747536955,
                    ],
                    im: [
                        0,
                        0,
                        0.9921603654457979,
                        0,
                        -0.2793462818542693,
                        1.7579030898107073,
                        -0.4527839725531578,
                        -0.8320432961178319,
                        -1.166570547084707,
                    ]
                }
            }
        },
        lowerBand0: {
            desc: 'lowerBand (default argument)',
            input: {
                data: matrix_mxn(3, 3, 3, 32).real()
            },
            expect: {
                o: [
                    1.2629542848807933,
                    -0.3262333607056494,
                    1.3297992629225006,
                    0.4146414344564082,
                    -1.5399500419037095,
                    0,
                    -0.005767172747536955,
                    0,
                    0,
                ],
            }
        },
        packedUpperBand0: {
            desc: 'packedLowerBand (default argument)',
            input: {
                data: matrix_mxn(3, 3, 3, 32).real()
            },
            expect: {
                o: [
                    1.2629542848807933,
                    1.2724293214294047,
                    0.4146414344564082,
                    -0.9285670347135381,
                    -0.2947204467905602,
                    -0.005767172747536955,
                ],
            }
        },
        packedUpperBand1: {
            desc: 'packedLowerBand (default argument), 32 bits',
            input: {
                data: matrix_mxn(3, 3, 3, 32).real()
            },
            expect: {
                o: [
                    1.2629542848807933,
                    -0.3262333607056494,
                    1.3297992629225006,
                    0.4146414344564082,
                    -1.5399500419037095,
                    0,
                    -0.005767172747536955,
                    0,
                    0,
                ],
            }
        }
    },
    mimicMatrixErr: {
        case0: {
            desc: 'lda < 0',
            input: {
                lda: -3,
                nrCols: 3,
                rowBase: 1,
                colBase: 1
            },
        },
        case1: {
            desc: 'lda != IsInteger',
            input: {
                lda: 0.2,
                nrCols: -3,
                rowBase: 1,
                colBase: 1
            },
        },
        case2: {
            desc: 'rowBase != Integer',
            input: {
                lda: 3,
                nrCols: 3,
                rowBase: 1.2,
                colBase: 1
            },
        },
        case3: {
            desc: 'rowBase != Integer',
            input: {
                lda: 3,
                nrCols: 3,
                rowBase: 1.2,
                colBase: 1
            },
        },
        case4: {
            desc: 'nrCols < 0',
            input: {
                lda: 3,
                nrCols: -3,
                rowBase: 1,
                colBase: 1
            },
        },
        case5: {
            desc: 'nrCols != Integer',
            input: {
                lda: 3,
                nrCols: 3.2,
                rowBase: 1,
                colBase: 1
            },
        },
        case6: {
            desc: 'colBase != Integer',
            input: {
                lda: 3,
                nrCols: 3,
                rowBase: 1,
                colBase: 1.3
            },
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
