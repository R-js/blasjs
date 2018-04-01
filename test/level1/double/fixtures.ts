import { complex as _, fortranArrComplex64 as arr64 } from '../../../src/lib/f_func';

const pI = Infinity;
const nI = -Infinity;
const { PI, sin, cos, abs, sqrt } = Math;

const cospi = x => cos(PI * x);
const sinpi = x => sin(PI * x);

export const fixture = {
    isamax: {
        case0: {
            desc: 'isamx n=6, incx=1',
            input: {
                n: 6,
                cx: {
                    re: [1, 0, 3, 4, 5, 6]
                },
                incx: 1,
            },
            output: {
                max: 6
            },
        },
        case1: {
            desc: 'isamx n=3, incx=2',
            input: {
                n: 3,
                cx: {
                    re: [1, 2, 3, 4, 5, 6]
                },
                incx: 2,
            },
            output: {
                max: 5
            },
        },
        case2: {
            desc: 'ismax n=3, incx=-1',
            input: {
                n: 6,
                cx: {
                    re: [1, 2, 3, 4, 5, 6],
                },
                incx: -1,
            },
            output: {
                max: 0
            },
        },
        case3: {
            desc: 'isamx n=0, incx=1',
            input: {
                n: 0,
                cx: {
                    re: [1, 2, 3, 4, 5, 6],

                },
                incx: 1,
            },
            output: {
                max: 0
            },
        },
        case4: {
            desc: 'ismax n=1, incx=1',
            input: {
                n: 1,
                cx: {
                    re: [1, 0, 3, 4, 5, 6]
                },
                incx: 1,
            },
            output: {
                max: 1
            },
        },
    },
    sasum: {
        case0: {
            desc: 'sasum n=6, incx=1',
            input: {
                n: 6,
                cx: {
                    re: [1, 0, 3, 4, 5, 6]
                },
                incx: 1,
            },
            output: {
                sum: 1 + 3 + 4 + 5 + 6
            },
        },
        case1: {
            desc: 'scasum n=3, incx=2',
            input: {
                n: 3,
                cx: {
                    re: [1, 2, 3, 4, 5, 6]
                },
                incx: 2,
            },
            output: {
                sum: 1 + 3 + 5
            },
        },
        case2: {
            desc: 'scasum n=0, incx=2',
            input: {
                n: 0,
                cx: {
                    re: [1, 2, 3, 4, 5, 6],

                },
                incx: 2,
            },
            output: {
                sum: 0
            },
        },
        case3: {
            desc: 'scasum n=3, incx=0',
            input: {
                n: 3,
                cx: {
                    re: [1, 2, 3, 4, 5, 6]
                },
                incx: 0,
            },
            output: {
                sum: 0
            },
        },
    },
    // CY(IY) = CY(IY) + CA*CX(IX)
    saxpy: {
        case0: {
            desc: 'cy = cy + sa*cx, incx=1, incy=1',
            input: {
                //> set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,
                        0.070508391424576003, 0.129287735160946243, 1.715064986883281017,
                        0.460916205989202299, -1.265061234606533969, -0.686852851893526073,
                    -0.445661970099958060]
                },
                sa: 23,
                incx: 1,
                incy: 1,
                n: 10
            },
            output: {
                re: [28.487472905706035,
                    -7.7335447857132165,
                    32.14409136136664,
                    29.33638278430088,
                    9.666040727658336,
                    -33.70378597690203,
                    -20.896125592422173,
                    -8.04363151078942,
                    -0.8194978250868761,
                    54.861365973632914]
            }
        },
        case0a: {
            desc: 'test n < 4 , n=3, cy + sa*cx, incx=1, incy=1',
            input: {
                //> set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,

                    ]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,
                    ]
                },
                sa: 23,
                incx: 1,
                incy: 1,
                n: 3
            },
            output: {
                re: [28.487472905706035,
                    -7.7335447857132165,
                    32.14409136136664
                ]
            }
        },
        case1: {
            desc: 'test delooping, n%4=0,n = 8 cy = cy + sa*cx, incx=1, incy=1',
            input: {
                //> set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977
                    ]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,
                        0.070508391424576003, 0.129287735160946243, 1.715064986883281017,
                        0.460916205989202299, -1.265061234606533969
                    ]
                },
                sa: 23,
                incx: 1,
                incy: 1,
                n: 8
            },
            output: {
                re: [28.487472905706035,
                    -7.7335447857132165,
                    32.14409136136664,
                    29.33638278430088,
                    9.666040727658336,
                    -33.70378597690203,
                    -20.896125592422173,
                    -8.04363151078942,
                ]
            }
        },
        case2: {
            desc: 'cy = cy + sa*cx, incx=2, incy=2',
            input: {
                n: 5,
                //> set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,
                        0.070508391424576003, 0.129287735160946243, 1.715064986883281017,
                        0.460916205989202299, -1.265061234606533969, -0.686852851893526073,
                    -0.445661970099958060]
                },
                sa: 23,
                incx: 2,
                incy: 2
            },
            output: {
                re: [28.487472905706035,
                    -0.23017748948327996,
                    32.14409136136664,
                    0.070508391424576,
                    9.666040727658336,
                    1.715064986883281,
                    -20.896125592422173,
                    -1.265061234606534,
                    -0.8194978250868761,
                    -0.44566197009995806]
            }
        },
        case3: {
            desc: 'cy = cy + ca*cx, ca = 0',
            input: {
                n: 1,
                cy: { re: [1], im: [2] },
                cx: { re: [3], im: [4] },
                sa: 0,
                incx: 1,
                incy: 1,
            },
            output: { re: [1], im: [2] }
        },
        case4: {
            desc: 'cy = cy + ca*cx, n=0',
            input: {
                n: 0,
                cy: { re: [3], im: [4] },
                cx: { re: [9], im: [4] },
                sa: 1,
                incx: 1,
                incy: 1,
            },
            output: { re: [3], im: [4] }
        },
        case5: {
            desc: 'cy = cy + ca*cx, incx=-2, incy=-2',
            input: {
                n: 5,
                //> set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,
                        0.070508391424576003, 0.129287735160946243, 1.715064986883281017,
                        0.460916205989202299, -1.265061234606533969, -0.686852851893526073,
                    -0.445661970099958060]
                },
                sa: 23,
                incx: -2,
                incy: -2
            },
            output: {
                re: [28.487472905706035,
                    -0.23017748948327996,
                    32.14409136136664,
                    0.070508391424576,
                    9.666040727658336,
                    1.715064986883281,
                    -20.896125592422173,
                    -1.265061234606534,
                    -0.8194978250868761,
                    -0.44566197009995806]
            }
        },
    },
    //DONE: scnrm2(n: number, x: FortranArr, incx: number): number
    scnrm2: {
        case0: {
            desc: 'Conj(x)*x, incx=1, n=10',
            input: {
                x: {
                    re: [
                        0, 0, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ],
                    im:
                        [0.76359346114045956, 0, -1.14765700923635139,
                            -0.28946157368822334, -0.29921511789731614, -0.41151083279506701,
                            0.25222344815613229, -0.89192112728456863, 0.43568329935571865,
                            -1.23753842192995811
                        ]
                },
                incx: 1,
                n: 10
            },
            output: 4.1815805452999522
        },
        case1: {
            desc: 'Conj(x)*x, incx=1, n=0',
            input: {
                x: {
                    re: [
                        0, 0, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ],
                    im:
                        [0.76359346114045956, 0, -1.14765700923635139,
                            -0.28946157368822334, -0.29921511789731614, -0.41151083279506701,
                            0.25222344815613229, -0.89192112728456863, 0.43568329935571865,
                            -1.23753842192995811
                        ]
                },
                incx: 1,
                n: 0
            },
            output: 0
        },
        case2: {
            desc: 'Conj(x)*x, incx=0, n=10',
            input: {
                x: {
                    re: [
                        0, 0, 1.3297992629225006134,
                        1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,
                        -0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798
                    ],
                    im:
                        [0.76359346114045956, 0, -1.14765700923635139,
                            -0.28946157368822334, -0.29921511789731614, -0.41151083279506701,
                            0.25222344815613229, -0.89192112728456863, 0.43568329935571865,
                            -1.23753842192995811
                        ]
                },
                incx: 0,
                n: 10
            },
            output: 0
        },
    },
    scopy: {
        case0: {
            desc: 'n=4, cx={4}, cy={4}, incx=1, incy=1',
            input: {
                n: 4,
                x: {
                    re: [1, 2, 3, 4]
                },
                y: {
                    re: [0, 0, 0, 0]
                },
                incx: 1,
                incy: 1,
            },
            output: {
                re: [1, 2, 3, 4]
            },
        },
        case1: {
            desc: 'n=10,(>7), incx=1, incy=1',
            input: {
                n: 10,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5],

                },
                y: {
                    re: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                },
                incx: 1,
                incy: 1,
            },
            output: {
                re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5]
            },
        },
        case2: {
            desc: 'n=7,( n % 7 == 0 ), incx=1, incy=1',
            input: {
                n: 7,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7]
                },
                y: {
                    re: [0, 0, 0, 0, 0, 0, 0]
                },
                incx: 1,
                incy: 1,
            },
            output: {
                re: [1, 2, 3, 4, 5, 6, 7]
            },
        },
        case3: {
            desc: 'n=5, incx=2, incy=2',
            input: {
                n: 5,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5],

                },
                y: {
                    re: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                },
                incx: 2,
                incy: 2,
            },
            output: {
                re: [1, 0, 3, 0, 5, 0, 7, 0, 9, 0]
            },
        },
        case4: {
            desc: 'n=5, incx=2, incy=-1',
            input: {
                n: 5,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5],

                },
                y: {
                    re: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                },
                incx: 2,
                incy: -1,
            },
            output: {
                re: [9, 7, 5, 3, 1, 0, 0, 0, 0, 0]
            },
        },
        case5: {
            desc: 'n=5, incx=-2, incy=1',
            input: {
                n: 5,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5],

                },
                y: {
                    re: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                },
                incx: -2,
                incy: 1,
            },
            output: {
                re: [9, 7, 5, 3, 1, 0, 0, 0, 0, 0]
            },
        },
        case6: {
            desc: 'n=0, incx=-2, incy=1',
            input: {
                n: 0,
                x: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 5],

                },
                y: {
                    re: [0, 0, 0, 0, 5, 0, 9, 0, 0, 0],
                },
                incx: -2,
                incy: 1,
            },
            output: {
                re: [0, 0, 0, 0, 5, 0, 9, 0, 0, 0],
            },
        },
    },
    sdot: {
        /*sdot(
                n: number,
                sx: FortranArr, 
                incx: number, 
                sy: FortranArr,
                incy: number
            ): number */
        case0: {
            desc: 'n=4, sx={4}, sy={4}, incx=1, incy=1',
            input: {
                n: 4,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: 1,
                incy: 1,
            },
            output: 70
        },
        case1: {
            desc: 'n=10, incx=1, incy=1',
            input: {
                n: 10,
                sx: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 2]
                },
                sy: {
                    re: [5, 6, 7, 8, 5, 4, 3, 3, -1, -9]
                },
                incx: 1,
                incy: 1,
            },
            output: 1 * 5 + 2 * 6 + 3 * 7 + 4 * 8 + 5 * 5 + 6 * 4 + 7 * 3 + 8 * 3 + 9 * -1 + 2 * -9
        },
        case2: {
            desc: 'n=0, incx=1, incy=1',
            input: {
                n: 0,
                sx: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 2]
                },
                sy: {
                    re: [5, 6, 7, 8, 5, 4, 3, 3, -1, -9]
                },
                incx: 1,
                incy: 1,
            },
            output: 0
        },
        case3: {
            desc: 'n=5, incx=2, incy=2',
            input: {
                n: 5,
                sx: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 2]
                },
                sy: {
                    re: [5, 6, 7, 8, 5, 4, 3, 3, -1, -9]
                },
                incx: 2,
                incy: 2,
            },
            output: 1 * 5 + 3 * 7 + 5 * 5 + 7 * 3 - 9 * 1
        },
        case4: {
            desc: 'n=5, incx=-2, incy=-2',
            input: {
                n: 5,
                sx: {
                    re: [1, 2, 3, 4, 5, 6, 7, 8, 9, 2]
                },
                sy: {
                    re: [5, 6, 7, 8, 5, 4, 3, 3, -1, -9]
                },
                incx: -2,
                incy: -2,
            },
            output: 1 * 5 + 3 * 7 + 5 * 5 + 7 * 3 - 9 * 1
        },
        case5: {
            desc: 'n=7, incx=1, incy=1',
            input: {
                n: 7,
                sx: {
                    re: [1, 2, 3, 4, 5, 6, 7]
                },
                sy: {
                    re: [5, 6, 7, 8, 5, 4, 3]
                },
                incx: 1,
                incy: 1,
            },
            output: 1 * 5 + 2 * 6 + 3 * 7 + 4 * 8 + 5 * 5 + 6 * 4 + 7 * 3
        }
    },
    sdsdot: {
        /*  n: number,
            sb: number,
            sx: FortranArr,
            incx: number,
            sy: FortranArr,
            incy: number */
        case0: {
            desc: 'n=4, sx={4}, sy={4}, incx=1, incy=1',
            input: {
                n: 4,
                sb: 2,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: 1,
                incy: 1,
            },
            output: (2) + (1 * 5) + (2 * 6) + (3 * 7) + (4 * 8)
        },
        case1: {
            desc: 'n=0, sx={4}, sy={4}, incx=1, incy=1',
            input: {
                n: 0,
                sb: 2,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: 1,
                incy: 1,
            },
            output: 2
        },
        case2: {
            desc: 'n=2, sx={4}, sy={4}, incx=2, incy=2',
            input: {
                n: 2,
                sb: 3,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: 2,
                incy: 2,
            },
            output: (3) + (1 * 5) + (3 * 7)
        },
        case3: {
            desc: 'n=2, sx={4}, sy={4}, incx=-1, incy=-1',
            input: {
                n: 2,
                sb: 3,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: -1,
                incy: -1,
            },
            output: (3) + (1 * 5) + (2 * 6)
        },
        case4: {
            desc: 'n=2, sx={4}, sy={4}, incx=1, incy=-1',
            input: {
                n: 2,
                sb: 3,
                sx: {
                    re: [1, 2, 3, 4]
                },
                sy: {
                    re: [5, 6, 7, 8]
                },
                incx: -1,
                incy: 1,
            },
            output: (3) + (1 * 6) + (2 * 5)
        },
    },
    //export function snrm2(n: number, x: FortranArr, incx: number): number
    snrm2: {
        case0: {
            desc: 'n=4, x={4}, incx=1',
            input: {
                n: 4,
                x: {
                    re: [3, 2, 3, 0]
                },

                incx: 1,
            },
            output: sqrt([3, 2, 3, 0].map(x => x * x).reduce((sum, v) => sum + v, 0))
        },
        case1: {
            desc: 'n=4, x={4}, incx=-1',
            input: {
                n: 4,
                x: {
                    re: [1, 2, 3, 4]
                },

                incx: -1,
            },
            output: 0
        },
        case2: {
            desc: 'n=1, x={4}, incx=1',
            input: {
                n: 1,
                x: {
                    re: [3, 2, 3, 4]
                },

                incx: 1,
            },
            output: 3
        },
    },
    srot: {
        /*
      n: number,
      sx: FortranArr,
      incx: number,
      sy: FortranArr,
      incy: number,
      c: number,
      s: number
      */
        case0: {
            desc: 'n=4, x={2},y={2} incx=1, incy=1',
            input: {
                n: 4,
                x: {
                    re: [1, 0]
                },
                incx: 1,
                y: {
                    re: [0, 1]
                },
                incy: 1,
                c: cospi(1 / 6),
                s: sinpi(1 / 6)
            },
            output: {
                x: [0.86602540378443871, 0.5],
                y: [-0.5, 0.86602540378443871]
            }
        },
        case1: {
            desc: 'n=4,  incx=-1',
            input: {
                n: 4,
                x: {
                    re: [1, 0]
                },
                incx: -1,
                y: {
                    re: [0, 1]
                },
                incy: -1,
                c: cospi(1 / 6),
                s: sinpi(1 / 6)
            },
            output: {
                x: [0.86602540378443871, 0.5],
                y: [-0.5, 0.86602540378443871]
            }
        },
        case2: {
            desc: 'n=4,  incx=-1',
            input: {
                n: 0,
                x: {
                    re: [1, 0]
                },
                incx: -1,
                y: {
                    re: [0, 1]
                },
                incy: -1,
                c: cospi(1 / 6),
                s: sinpi(1 / 6)
            },
            output: {
                x: [1, 0],
                y: [0, 1]
            }
        },
    },
    //srotg(p: { sa: number, sb: number, c: number, s: number })
    srotg: {
        case0: {
            desc: 'sa=4, sb=2, => c= s=',
            input: {
                sa: 4,
                sb: 2
            },
            output: {
                sa: 4.47213595499958,
                sb: 0.4472135954999579,
                c: 0.8944271909999159,
                s: 0.4472135954999579
            }
        },
        case1: {
            desc: 'sa=0, sb=0, => c=1 s=0',
            input: {
                sa: 0,
                sb: 0
            },
            output: {
                sa: 0,
                sb: 0,
                c: 1,
                s: 0
            }
        },
        case3: {
            desc: 'sa=4, sb=2, => c= s=',
            input: {
                sa: -2,
                sb: -4
            },
            output: {
                sa: -4.47213595499958,
                sb: 2.2360679774997898,
                s: 0.8944271909999159,
                c: 0.4472135954999579
            }
        },
    },
    srotm: {
        case0: {
            desc: 'n=4, sx={4}, incx=1, sy={4},incy=1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 1,
                incy: 1,
                sparam: [-2, 0, 0, 0, 0]
            },
            //result values from fortran code
            output: {
                x: [1, 2, 3, 4],
                y: [5, 6, 7, 8]
            }
        },
        case1: {
            desc: 'n=0, sx={4}, incx=1, sy={4},incy=1, sparam={5}',
            input: {
                n: 0,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 1,
                incy: 1,
                sparam: [-2, 0, 0, 0, 0]
            },
            //result values from fortran code
            output: {
                x: [1, 2, 3, 4],
                y: [5, 6, 7, 8]
            }
        },
        case2: {
            desc: 'n=4, sx={4}, incx=1, sy={4},incy=1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 1,
                incy: 1,
                sparam: [-1, 2, 3, 4, 5]
            },
            //result values from fortran code
            output: {
                x: [22, 28, 34, 40],
                y: [28, 36, 44, 52]
            }
        },
        case3: {
            desc: 'n=4, sx={4}, incx=1, sy={4},incy=1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 1,
                incy: 1,
                sparam: [0, NaN, 3, 4, NaN]
            },
            //result values from fortran code
            output: {
                x: [21, 26, 31, 36],
                y: [8, 12, 16, 20]
            }
        },
        case4: {
            desc: 'n=4, sx={4}, incx=1, sy={4},incy=1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 1,
                incy: 1,
                sparam: [1, 2, NaN, NaN, 3]
            },
            //result values from fortran code
            output: {
                x: [7, 10, 13, 16],
                y: [14, 16, 18, 20]
            }
        },
        case5: {
            desc: 'n=4, sx={4}, incx=-1, sy={4}, incy=1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: -1,
                incy: 1,
                sparam: [1, 2, NaN, NaN, 3]
            },
            //result values from fortran code
            output: {
                x: [10, 11, 12, 13],
                y: [11, 15, 19, 23]
            }
        },
        case6: {
            desc: 'n=4, sx={4}, incx=1, sy={4},incy=-1, sparam={5}',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: +1,
                incy: -1,
                sparam: [1, 2, NaN, NaN, 3]
            },
            //result values from fortran code
            output: {
                x: [10, 11, 12, 13],
                y: [11, 15, 19, 23]
            }
        },
        case7: {
            desc: 'n=2, sx={4}, incx=+1, sy={4},incy=-2, sparam={5}',
            input: {
                n: 2,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: +2,
                incy: -2,
                sparam: [0, NaN, 4, 5, NaN]
            },
            //result values from fortran code
            output: {
                x: [36, 2, 28, 4],
                y: [17, 6, 11, 8]
            }
        },
        case8: {
            desc: 'n=2, sx={4}, incx=2, sy={4},incy=2, sparam={5}',
            input: {
                n: 2,
                sx: [1, 2, 3, 4],
                sy: [5, 6, 7, 8],
                incx: 2,
                incy: -2,
                sparam: [-1, 1, 2, 3, 4]
            },
            //result values from fortran code
            output: {
                x: [22, 2, 18, 4],
                y: [26, 6, 30, 8]
            }
        },
    },
    srotmg: {
        case0: {
            desc: 'sd1=-4, sd2=2, sx1=3, sx2=9',
            input: {
                sd1: -4,
                sd2: 2,
                sx1: 3,
                sy1: 9,
                sparam: [0, 0, 0, 0, 0]
            },
            //result values from fortran code
            output: {
                sd1: 0,
                sd2: 0,
                sx1: 0,
                sy1: 9,
                sparam: [-1, 0, 0, 0, 0]
            }
        },
        case1: {
            desc: 'sd1=4, sd2=0, sx1=3, sx2=9',
            input: {
                sd1: 4,
                sd2: 0,
                sx1: 3,
                sy1: 9,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 4,
                sd2: 0,
                sx1: 3,
                sy1: 9.00000000,
                sparam: [-2, 0, 0, 0, 0]
            }
        },
        case2: {
            desc: 'sd1=1, sd2=2, sx1=3, sx2=1',
            input: {
                sd1: 1,
                sd2: 2,
                sx1: 3,
                sy1: 1,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 0.81818181818181812,
                sd2: 1.6363636363636362,
                sx1: 3.6666666666666670,
                sy1: 1,
                sparam: [0, 0, -0.33333333333333331, 0.66666666666666663, 0]
            }
        },
        case3: {
            desc: 'sd1=2, sd2=-1, sx1=3, sx2=8',
            input: {
                sd1: 2,
                sd2: -1,
                sx1: 3,
                sy1: 8,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 0,
                sd2: 0,
                sx1: 0,
                sy1: 8,
                sparam: [-1, 0, 0, 0, 0]
            }
        },
        case4: {
            desc: 'sd1=2, sd2=1, sx1=3, sx2=8',
            input: {
                sd1: 2,
                sd2: 1,
                sx1: 3,
                sy1: 8,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 0.78048780487804881,
                sd2: 1.5609756097560976,
                sx1: 10.250000000000000,
                sy1: 8.0000000000000000,
                sparam: [1, 0.75, 0, 0, 0.375]
            }
        },
        case5: {
            desc: 'sd1= 1 / (2 << 23), sd2=1 / (2 << 24), sx1=3, sx2=8',
            input: {
                sd1: 1 / (2 << 23),
                sd2: 1 / (2 << 24),
                sx1: 3,
                sy1: 8,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 0.39024390243902440,
                sd2: 0.78048780487804881,
                sx1: 2.5024414062500000E-003,
                sy1: 8,
                sparam: [-1, 1.8310546875000000E-004, -2.4414062500000000E-004, 1, 9.1552734375000000E-005]
            }
        },
        case6: {
            desc: 'sd1=1/(2 << 23), sd2=1/(2 << 24), sx1=3, sx2=2',
            input: {
                sd1: 1 / (2 << 23),
                sd2: 1 / (2 << 24),
                sx1: 3,
                sy1: 2,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 0.81818181818181812,
                sd2: 0.40909090909090906,
                sx1: 8.9518229166666674E-004,
                sy1: 2,
                sparam: [-1, 2.4414062500000000E-004, -2.4414062500000000E-004, 1, 0]
            }
        },
        case7: {
            desc: 'sd1=2/(2 << 23), sd2=1/(2 << 23), sx1=3, sx2=2',
            input: {
                sd1: 2 / (2 << 23),
                sd2: 1 / (2 << 23),
                sx1: 3,
                sy1: 2,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 9.7534873268821016E-008,
                sd2: 0.81818181818181812,
                sx1: 3.6666666666666670,
                sy1: 2,
                sparam: [-1, 1, -1.6276041666666666E-004, 0.33333333333333331, 2.4414062500000000E-004]
            },
        },
        case8: {
            desc: 'sd1=(2 << 23), sd2=2 * (2 << 23), sx1=3, sx2=2',
            input: {
                sd1: (2 << 23),
                sd2: 2 * (2 << 23),
                sx1: 3,
                sy1: 2,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 8882055.5294117648,
                sd2: 1.0588235294117647,
                sx1: 5.6666666666666661,
                sy1: 2,
                sparam: [-1, 1, -2730, 1.3333333333333333, 4096]
            }
        },
        case9: {
            desc: 'sd1=1, sd2=-1, sx1=3, sx2=2',
            input: {
                sd1: 2 * (2 << 23),
                sd2: 1 * (2 << 23),
                sx1: 3,
                sy1: 2,
                sparam: [0, 0, 0, 0, 0]
            },
            output: {
                sd1: 1.6363636363636362,
                sd2: 13726813.090909090,
                sx1: 15018.666666666668,
                sy1: 2,
                sparam: [-1, 4096, -0.66666666666666663, 1365.3333333333333, 1]
            }
        }
    },
    sscal: {
        case0: {
            desc: 'n=4, sx={4}, sa=1, incx=1',
            input: {
                n: 4,
                sa: 1,
                sx: [1, 2, 3, 4],
                incx: 1
            },
            //result values from fortran code
            output: {
                x: [1, 2, 3, 4],
            }
        },
        case1: {
            desc: 'n=4, sx={4}, sa=2, incx=1',
            input: {
                n: 4,
                sa: 2,
                sx: [1, 2, 3, 4],
                incx: 1
            },
            //result values from fortran code
            output: {
                x: [2, 4, 6, 8],
            }
        },
        case2: {
            desc: 'n=5, sx={5}, sa=2, incx=1',
            input: {
                n: 5,
                sa: 2,
                sx: [3, 1, 2, 3, 4],
                incx: 1
            },
            //result values from fortran code
            output: {
                x: [6, 2, 4, 6, 8],
            }
        },
        case4: {
            desc: 'n=7, sx={4}, sa=2, incx=1',
            input: {
                n: 7,
                sa: 2,
                sx: [-1, 2, 3, 1, 2, 3, 4],
                incx: 1
            },
            //result values from fortran code
            output: {
                x: [-2, 4, 6, 2, 4, 6, 8],
            }
        },
        case5: {
            desc: 'n=3, sx={7}, sa=2, incx=2',
            input: {
                n: 3,
                sa: 2,
                sx: [-1, 2, 3, 1, 2, 3, 4],
                incx: 2
            },
            //result values from fortran code
            output: {
                x: [-2, 2, 6, 1, 4, 3, 4],
            }
        },
        case6: {
            desc: 'n=6, sx={4}, sa=2, incx=-2',
            input: {
                n: 6,
                sa: 2,
                sx: [-1, 2, 3, 1, 2, 3, 4],
                incx: -2
            },
            //result values from fortran code
            output: {
                x: [-1, 2, 3, 1, 2, 3, 4],
            }
        },
        case7: {
            desc: 'n=0, sx={4}, sa=2, incx=-2',
            input: {
                n: 0,
                sa: 2,
                sx: [-1, 2, 3, 1, 2, 3, 4],
                incx: -2
            },
            //result values from fortran code
            output: {
                x: [-1, 2, 3, 1, 2, 3, 4],
            }
        },
    },
    sswap: {
        case0: {
            desc: 'n=4, sx={4}, incx=1, incy=-1',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                incx: 1,
                sy: [1, 2, 3, 4],
                incy: -1 //flip when copying to x, flip when assigning to y
            },
            //result values from fortran code
            output: {
                x: [4, 3, 2, 1],
                y: [4, 3, 2, 1]
            }
        },
        case1: {
            desc: 'n=4, sx={4}, incx=-1, incy=1',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                incx: -1, //flip when copying to y, flip when assigning to x
                sy: [1, 2, 3, 4],
                incy: 1
            },
            //result values from fortran code
            output: {
                x: [4, 3, 2, 1],
                y: [4, 3, 2, 1]
            }
        },
        case2: {
            desc: 'n=4, sx={4}, incx=1, incy=1',
            input: {
                n: 4,
                sx: [1, 2, 3, 4],
                incx: 1,
                sy: [7, 8, 9, 5],
                incy: 1
            },
            //result values from fortran code
            output: {
                x: [7, 8, 9, 5],
                y: [1, 2, 3, 4]
            }
        },
        case3: {
            desc: 'n=3, sx={4}, incx=1, incy=1',
            input: {
                n: 3, //batchsize
                sx: [1, 2, 3],
                incx: 1,
                sy: [7, 8, 9],
                incy: 1
            },
            //result values from fortran code
            output: {
                x: [7, 8, 9],
                y: [1, 2, 3]
            }
        },
        case4: {
            desc: 'n=3, sx={4}, incx=1, incy=1',
            input: {
                n: 2, //batchsize
                sx: [1, 2],
                incx: 1,
                sy: [7, 8],
                incy: 1
            },
            //result values from fortran code
            output: {
                x: [7, 8],
                y: [1, 2]
            }
        },
        case5: {
            desc: 'n=0, sx={4}, incx=1, incy=1',
            input: {
                n: 0, //batchsize
                sx: [1, 2],
                incx: 1,
                sy: [7, 8],
                incy: 1
            },
            //result values from fortran code
            output: {
                x: [1, 2],
                y: [7, 8]
            }
        },
    }
    /*
    cdotu: {
        case0: {
            desc: 'cdotu, n=4, cx={4}, cy={4}, incx=1, incy=1',
            input: {
                n: 4,
                cx: {
                    re: [1, 2, 3, 4],
                    im: [5, 6, 7, 8]
                },
                cy: {
                    re: [5, 6, 7, 8],
                    im: [9, 10, 11, 12]
                },
                incx: 1,
                incy: 1,
            },
            output: {
                re: -208,
                im: 284
            }
        },
        case1: {
            desc: 'cdotu, n=4, cx={4}, cy={4}, incx=-1, incy=-1',
            input: {
                n: 4,
                cx: {
                    re: [1, 2, 3, 4],
                    im: [5, 6, 7, 8]
                },
                cy: {
                    re: [5, 6, 7, 8],
                    im: [9, 10, 11, 12]
                },
                incx: -1,
                incy: -1,
            },
            output: {
                re: -208,
                im: 284
            }
        },
        case2: {
            desc: 'cdotu, n=0, cx={4}, cy={4}, incx=1, incy=1',
            input: {
                n: 0,
                cx: {
                    re: [1, 2, 3, 4],
                    im: [5, 6, 7, 8]
                },
                cy: {
                    re: [5, 6, 7, 8],
                    im: [9, 10, 11, 12]
                },
                incx: 1,
                incy: 1,
            },
            output: {
                re: 0,
                im: 0
            }
        }
    },
    cdotuError: {
        case0: {
            desc: 'missing cx imaginary part will throw',
            input: {
                n: 4,
                cx: {
                    re: [1, 2, 3, 4],
                },
                cy: {
                    re: [5, 6, 7, 8],
                    im: [9, 10, 11, 12]
                },
                incx: 1,
                incy: 1,
            }
        },
        case1: {
            desc: 'missing cy imaginary part will throw',
            input: {
                n: 4,
                cx: {
                    re: [1, 2, 3, 4],
                    im: [9, 10, 11, 12]
                },
                cy: {
                    re: [5, 6, 7, 8],
                },
                incx: 1,
                incy: 1,
            }
        }
    },
 
    crotg: {
        case0: {
            desc: 'crotg givens rotation ca=(11,19), cb=(34,23), c=>0.471621990, s=>(0.793538332,0.384538263)',
            input: {
                ca: { re: 11, im: 19 },
                cb: { re: 34, im: 23 },
            },
            output: {
                ca: { re: 23.323763103564644, im: 40.286499906157111 },
                cb: { re: 34.000000, im: 23.0000000 },
                c: 0.47162200847078728,
                s: { re: 0.79353827566350310, im: 0.38453827661622281 }
            }
        },
        case1: {
            desc: 'cdotu, n=4, ca=0, should return 0',
            input: {
                ca: { re: 0, im: 0 },
                cb: { re: -1, im: -1 }
            },
            output: {
                ca: { re: -1, im: -1 }, //copy of cb
                cb: { re: -1, im: -1 },
                c: 0,
                s: { re: 1, im: 0 }
            }
        }
    },
    cscal: {
        case0: {
            desc: 'cscal n=6, ca={2,4}, cx={6}',
            input: {
                n: 6,
                ca: { re: 2, im: 4 },
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
                incx: 1
            },
            output: {
                cx: {
                    re: [6.24017670861573848, 0.52641506575094199, 2.68266721683514886,
                        -7.07375491257299416, -2.22509097564902181, 0.11613691215005373],
                    im: [3.1946830700960973, -1.8943743364037180, 5.3076627061949289,
                        9.8990240634335205, 3.1857526601065520, -7.7578186655935744]
                }
            },
        },
        case1: {
            desc: 'cscal n=3, incx=2, ca={2,4}, cx={6}',
            input: {
                n: 3,
                ca: { re: 2, im: 4 },
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
                incx: 2
            },
            output: {
                cx: {
                    re: [6.24017670861573848, -0.32623336070564940, 2.68266721683514886,
                        1.27242932142940468, -2.22509097564902181, -1.53995004190370954],
                    im: [3.1946830700960973, -0.2947204467905601977, 5.3076627061949289,
                        2.4046533888579508798, 3.1857526601065520, -0.7990092489893682037]
                }
            },
        },
        case2: {
            desc: 'cscal n=0, cx unchanged',
            input: {
                n: 0,
                ca: { re: 2, im: 4 },
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
                incx: 2
            },
            output: {
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
            },
        },
        case3: {
            desc: 'cscal incx=0 cx unchanged',
            input: {
                n: 3,
                ca: { re: 2, im: 4 },
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
                incx: 0
            },
            output: {
                cx: {
                    re: [1.26295428488079331, -0.32623336070564940, 1.32979926292250061,
                        1.27242932142940468, 0.41464143445640822, -1.53995004190370954],
                    im: [-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,
                        2.4046533888579508798, 0.7635934611404595618, -0.7990092489893682037]
                },
            },
        }
    },
    cscalError: {
        case0: {
            desc: 'cscal imaginary part of cx will throw',
            input: {
                n: 6,
                ca: { re: 2, im: 4 },
                cx: {
                    re: [1.26295428488079331]
                },
                incx: 1
            }
        },
    },
    */

}
