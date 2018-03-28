const pI = Infinity;
const nI = -Infinity;
import { complex as _, fortranArrComplex64 as arr64 } from '../src/lib/f_func';
export const fixture = {
    caxpy: {
        case0: {
            desc: 'cy = cy + ca*cx, incx=1',
            input: {
                //> > set.seed(0)
                //> complex(r=rnorm(10),i=rnorm(10))
                cx: {
                    re: [
                        1.262954285, -0.326233361, 1.329799263, 1.272429321, 0.414641434,
                        -1.539950042, -0.928567035, -0.294720447, -0.005767173, 2.404653389
                    ],
                    im:
                        [0.7635935, -0.7990092, -1.1476570, -0.2894616, -0.2992151, -0.4115108,
                            0.2522234, -0.8919211, 0.4356833, -1.2375384]
                },
                //set.seed(123) (Mersenne-Twister, Inversion)
                cy: {
                    re: [-0.56047565, -0.23017749, 1.55870831, 0.07050839, 0.12928774, 1.71506499,
                        0.46091621, -1.26506123, -0.68685285, -0.44566197],
                    im: [1.2240818, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.7869131,
                        0.4978505, -1.9666172, 0.7013559, -0.4727914]
                },
                ca: { re: 23, im: 32 },
                incx: 1,
                incy: 1
            },
            output: {
                re: [4.052482, 17.834751, 68.869116, 38.599153, 19.240925, -20.535439,
                    -28.967276, 20.497845, -14.761363, 94.462595],
                im: [59.201269, -28.456866, 16.558237, 34.170805, 5.830737, -56.956237,
                    -23.415155, -31.911857, 10.537522, 48.012733]
            }
        },
        case1: {
            desc: 'cy = cy + ca*cx, incx=3',
            input: {
                cx: [],
                cy: [],
                ca: { re: 23, im: 32 },
                incx: 2,
                incy: 1
            },
            output: []
        },
        case2: {
            desc: 'cy = cy + ca*cx, incx=0',
            input: {
                x: [-1, 0.5, 0.5, 0, 0.5, 0.5, 0, 0.5, 1, 0.5, 0],
                shape1: [2, -2, 1, 0, 0, 8, 8, 8, pI, pI, pI],
                shape2: [2, 2, -1, 0, 0, pI, pI, pI, 4, 4, pI],
                ncp: 0,
                asLog: false
            },
            output: [
                0,
                NaN,
                NaN,
                pI,
                0, //5 
                0,
                pI,
                0,
                pI,
                0, //10 
                0
            ]
        },
        case3: {
            desc: 'edge cases x = (0,1), shape1 = Infinity , shape2 = Infinity',
            input: {
                x: [0.5, 0, 0, 1, 1, 1],
                shape1: [pI, 0.5, 1, 3, 3, 9],
                shape2: [pI, 1, 8, 2, 0.5, 1],
                ncp: 0,
                asLog: false
            },
            output: [
                pI, pI, 8, 0, pI, 9
            ]
        },
        case4: {
            desc: 'shape1 <=2, shape2 <=2, x =(0,..,1), ncp=0, log=false',
            input: {
                x: [0.5, 0.2],
                shape1: [4, 6],
                shape2: [1.5, 6],
                ncp: 0,
                asLog: false
            },
            output: [
                0.8700727972, 0.2906652672
            ]
        }
    },
    ccopy: {
        case0: {
            desc: 's1=2, s2=5, ncp=0, l:false, x=0.5',
            input: {
                //n: number,
                //cx: FortranArr,
                //incx: number,
                //cy: FortranArr,
                //incy: number): 
            },
            output: 0.890625
        },
        case1: {
            desc: 's1=2, s2=5, ncp=undef, l=false, x=[0, 0.2, 0.4, 0.6, 0.8, 1]',
            input: {
            },
            output: []
        },
        case2: {
            desc: 's1=2, s2=5, ncp=undef, lt=false, x=[0, 0.2, 0.4, 0.6, 0.8, 1]',
            input: {

            },
            output: []
        },
        case3: {
            desc: 's1=2, s2=5, ncp=undef, lt=true, asl=true, x=[0, 0.2, 0.4, 0.6, 0.8, 1]',
            input: {
            },
            output: []
        },
    }
}
