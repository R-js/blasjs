import { complex as _, fortranArrComplex64 as arr64 } from '../../../src/lib/f_func';

const pI = Infinity;
const nI = -Infinity;
const { PI, sin, cos, abs, sqrt } = Math;

const cospi = x => cos(PI * x);
const sinpi = x => sin(PI * x);

export const fixture = {
    sgbmv: {
        case0: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", kl=1,ku=1,m=6,n=9, alpha=1.5, beta=2',
            input: {
                trans: 'n',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                alpha: 1.5,
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                lda: 3,
                x: [1, 1, 2, 2, 3, 3, 4, 4, 5],
                incx: 1,
                beta: 2.5,
                y: [1, 1, 1, 1, 1, 1],
                incy: 1
            },
            output: {
                y: [4, 1, 4, 1, 4, 1]
            },
        },
        case1: {
            desc: 'sgbmv, y := alpha*A*x + beta*y,n=0',
            input: {
                trans: 'n',
                m: 6,
                n: 0,
                kl: 1,
                ku: 1,
                lda: 3,
                incy: 1,
                incx: 1,
                //
                beta: 2.5,
                alpha: 1.5,

                a: [
                ],
                x: [],

                y: [1, 1, 1, 1, 1, 1],
            },
            output: {
                y: [1, 1, 1, 1, 1, 1]
            },
        },
        case2: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, alpha=0 && beta=1',
            input: {
                trans: 'n',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                lda: 3,
                incy: 1,
                incx: 1,
                //
                beta: 1,
                alpha: 0,

                a: [
                ],
                x: [],

                y: [1, 1, 1, 1, 1, 1],
            },
            output: {
                y: [1, 1, 1, 1, 1, 1]
            },
        },
        case3: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", kl=1,ku=1,m=3,n=9, incy=2, alpha=1.5, beta=2',
            input: {
                trans: 'n',
                m: 3,
                n: 9,
                kl: 1,
                ku: 1,
                alpha: 1.5,
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                lda: 3,
                x: [1, 1, 2, 2, 3, 3, 4, 4, 5],
                incx: 1,
                beta: 2.5,
                y: [1, 1, 1, 1, 1, 1],
                incy: 2
            },
            output: {
                y: [4, 1, 1, 1, 4, 1]
            },
        },
        case4: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", kl=1,ku=1,m=3,n=9, incy=2, alpha=1.5, beta=0',
            input: {
                trans: 'n',
                m: 3,
                n: 9,
                kl: 1,
                ku: 1,
                alpha: 1.5,
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                lda: 3,
                x: [1, 1, 2, 2, 3, 3, 4, 4, 5],
                incx: 1,
                beta: 0,
                y: [1, 1, 1, 1, 1, 1],
                incy: 2
            },
            output: {
                y: [1.5, 1, -1.5, 1, 1.5, 1]
            },
        },
        case5: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", kl=1,ku=1,m=3,n=9, incy=1, alpha=1.5, beta=0',
            input: {
                trans: 'n',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                incx: 1,
                incy: 1,
                lda: 3,
                //
                alpha: 1.5,
                beta: 0,
                //
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                x: [1, 1, 2, 2, 3, 3, 4, 4, 5],
                y: [1, 1, 1, 1, 1, 1],
            },
            output: {
                y: [1.5, -1.5, 1.5, -1.5, 1.5, -1.5]
            },
        },
        case6: {
            desc: 'sgbmv, y := alpha*(A**T)x + beta*y, trans="t", kl=1,ku=1,m=6,n=9, incy=1, alpha=1.5, beta=1',
            input: {
                trans: 't',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                incx: 1,
                incy: 1,
                lda: 3,
                //
                alpha: 1.5,
                beta: 1,
                //
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                x: [1, 1, 2, 2, 3, 3],
                y: [1, 1, 1, 1, 1, 1, 1, 1, 1],
            },
            output: {
                y: [2.5, -0.5, 2.5, -0.5, 2.5, 5.5, -3.5, 1, 1]
            },
        },
        case7: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", kl=1,ku=1,m=6,n=9, alpha=1.5, beta=2',
            input: {
                trans: 'n',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                alpha: 0,
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                lda: 3,
                x: [1, 1, 2, 2, 3, 3, 4, 4, 5],
                incx: 1,
                beta: 2.5,
                y: [1, 1, 1, 1, 1, 1],
                incy: 1
            },
            output: {
                y: [2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
            },
        },
        case8: {
            desc: 'sgbmv, y := alpha*A*x + beta*y, trans="n", incx=-1,incy=-1 kl=1,ku=1,m=6,n=9, alpha=1, beta=0',
            input: {
                trans: 'n',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                alpha: 1,
                a: [
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                    -1, 2, -1,
                ],
                lda: 3,
                x: [1, 1, 2, 2, 3, 3, 4, 4, 3],
                incx: -1,
                beta: 0,
                y: [9, 0, 0, 0, 0, -9],
                incy: -1
            },
            output: {
                //incy=1, [ 1, -1, 1, -1, 1, -1 ],
                //incy=-1, [ -1, 1, -1, 1, -1, 1 ]
                //incx=-1, incy=1, [ 2, 1, 1, -1, 1, -1 ]
                //incx=-1, incy=-1, [ -1, 1, -1, 1, 1, 2 ]
                y: [-1, 1, -1, 1, 1, 2]
                //2, 1, 1, -1, 1, -1
            },
        },
    },
    sgbmvErrors: {
        case0: {
            desc: 'sgbmv, trans != ("n","t")',
            input: {
                trans: 'x',
                m: 6,
                n: 9,
                kl: 1,
                ku: 1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case1: {
            desc: 'sgbmv, n<0',
            input: {
                trans: 'n',
                m: 6,
                n: -1,
                kl: 1,
                ku: 1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case2: {
            desc: 'sgbmv, m<0',
            input: {
                trans: 'n',
                m: -6,
                n: 1,
                kl: 1,
                ku: 1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case3: {
            desc: 'sgbmv, kl<0',
            input: {
                trans: 'n',
                m: 6,
                n: 1,
                kl: -1,
                ku: 1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case4: {
            desc: 'sgbmv, ku<0',
            input: {
                trans: 'n',
                m: 6,
                n: 1,
                kl: 1,
                ku: -1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case5: {
            desc: 'sgbmv, ku<0',
            input: {
                trans: 'n',
                m: 6,
                n: 1,
                kl: 1,
                ku: 1,
                incy: 1,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 2,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case6: {
            desc: 'sgbmv, incy=0',
            input: {
                trans: 'n',
                m: 6,
                n: 1,
                kl: 1,
                ku: 1,
                incy: 0,
                incx: 1,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
        case7: {
            desc: 'sgbmv, incx=0',
            input: {
                trans: 'n',
                m: 6,
                n: 1,
                kl: 1,
                ku: 1,
                incy: 1,
                incx: 0,
                alpha: 1.5,
                a: [
                ],
                lda: 3,
                x: [],
                beta: 2.5,
                y: [],
            }
        },
    }
}
