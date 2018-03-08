//intrinsic routines of fortran

const { abs } = Math;

export function sign(a: number, b?: number): number {
    if (b === undefined) {
        return a;
    }
    const rc = Math.abs(a);
    return b >= 0 ? rc : -rc;
}

//make JS arrays work like fortran ones
//1. export type
//2. export array wrapper factory

//1
export type Complex = { re: number, im: number };
export type fpArray = Float32Array | Float64Array;
export type FortranSetterGetter = (index: number) => (value?: number) => number | Complex;
export type FortranArr = {
    base: number,
    r: fpArray,
    i?: fpArray,
    s: FortranSetterGetter
};
//2
export function mimicFArray(r: fpArray, i?: fpArray) {
    //Lets make some curry
    let func = function n(startIndex: number = 0): FortranArr {

        return Object.freeze({
            base: startIndex,
            r,
            i,
            s: (index: number) => (re?: number, im?: number) => {
                const pRe = r[index - startIndex];
                const pIm = im && i ? i[index - startIndex] : undefined;
                //check
                if (re !== undefined) {
                    r[index - startIndex] = re;
                }
                if (im !== undefined && i === undefined) {
                    throw new Error('You specified a complex number for a real array');
                }
                if (i !== undefined) {
                    r[index - startIndex] = im || 0;
                    return { re: pRe || 0, im: pIm || 0 };
                }
                return pRe;
            }
        });
    }
    //   TODO: maybe keep reference later, who knows
    //    func['buffer'] = arr;
    return func;
}

export function scabs1(c: Complex) {
    return abs(c.re) + abs(c.im);
}

export function scabs1A(re: number, im: number): number {
    return abs(re) + abs(im);
}

export function complex(re: number = 0, im: number = 0): Complex {
    return { re, im };
}

// dont sue
export const cmult = (reA: number, imA: number, reB: number, imB: number): Complex => {
    //   (a + bi)(c+di)= (a*c-b*d)+i(a*d+b*c)
    return {
        re: (reA * reB - imA * imB),
        im: (reA * imB + imA * reB)
    };
}

export const cabs = (reA: number, imA: number) => {
    return Math.sqrt(reA * reA + imA * imA);
}

export enum ERROR {
    MISSING_IMAGINARY = 1,
    MISSING_REAL = 2,
    INTERNAL_NOTFOUND = 99
}

// private
const errors = new Map<ERROR, string>([
    [ERROR.MISSING_IMAGINARY, 'Missing imaginary part for %s'],
    [ERROR.MISSING_REAL, 'Missing real part for %s'],

]);

const ERROR_UKNOWN = 'Unkown Error code used! [%s]';

export function errorMsg(errNo: ERROR, fmt?: string): string {

    let msg = errors.get(errNo);
    if (!msg) {
        msg = ERROR_UKNOWN;
    }
    return msg.replace('%s', fmt || '');
}

export function errMissingIm(fmt: string): string {
    return errorMsg(ERROR.MISSING_IMAGINARY, fmt);
}

