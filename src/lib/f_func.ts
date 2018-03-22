//intrinsic routines of fortran

const { abs } = Math;

const { isInteger } = Number;

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
    s: FortranSetterGetter,
    assertComplex: (msg: string) => void | never;
};
//2
export function mimicFArray(r: fpArray, i?: fpArray) {
    //Lets make some curry
    let func = function n(startIndex: number = 1): FortranArr {

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
            },
            assertComplex: (msg: string) => {
                if (i === undefined) {
                    throw new Error(errMissingIm(msg))
                }
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

// TODO: REMOVE THIS in next iteration!! explicitly it out  
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
    ERR_MISSING_IMAGINARY = 1,
    ERR_MISSING_REAL = 2,
    ERR_WRONG_ARGUMENT = 3,
}

// private
const errors = new Map<ERROR, string>([
    [ERROR.ERR_MISSING_IMAGINARY, 'Missing imaginary part for %s'],
    [ERROR.ERR_MISSING_REAL, 'Missing real part for %s'],
    [ERROR.ERR_WRONG_ARGUMENT, 'Argument [%s] is has invalid value:[%s]']
]);

const ERROR_UKNOWN = 'Unkown Error code used! [%s]';

export function errorMsg(errNo: ERROR, ...fmt: string[]): string {

    let msg = errors.get(errNo);
    if (!msg) {
        msg = ERROR_UKNOWN;
    }
    if (fmt.length === 0) {
        return msg;
    }
    return fmt.reduce((p, v, k) => p.replace('%s', v), msg);
}

//helpers for errorMsg

export function errMissingIm(fmt: string): string {
    return errorMsg(ERROR.ERR_MISSING_IMAGINARY, fmt);
}

export function errWrongArg(arg: string, value): string {
    return errorMsg(ERROR.ERR_WRONG_ARGUMENT, value);
}


// Matrix
// Matrix
// Matrix
// Matrix

export type FortranMatrixSetterGetter = (col: number) => (row?: number) => number | Complex;

export type Matrix = {
    readonly rowBase: number,
    readonly colBase: number,
    readonly nrCols: number, // inclusive note!!
    readonly nrRows: number,
    readonly r: fpArray, //[(ncols+1)*(nrows+1)]
    readonly i?: fpArray, //imaginary part of matrix [(ncols+1)*(nrows+1)]
    readonly coord: (number) => (number) => number;
    //readonly colOf: (number) => number;
    readonly colOfEx: (number) => number;
    //s: FortranMatrixSetterGetter
    readonly assertComplex: (msg: string) => void | never;
    // zap a row with fa valie
    readonly setCol: (col: number, rowStart: number, rowEnd: number, value: number) => void;
};


export function mimicFMatrix2D(r: fpArray, i?: fpArray) {

    return function c1(nrRows: number, nrCols: number, rowBase: number = 1, colBase = 1): Matrix {

        // check rows
        if (nrRows < 0) {
            throw new Error(errWrongArg('nrRows', 'is Negative'));
        }
        if (!isInteger(nrRows)) {
            throw new Error(errWrongArg('nrRows', 'is a NaN'));
        }
        if (!isInteger(rowBase)) {
            throw new Error(errWrongArg('rowBase', 'is a NaN'));
        }
        // check columns
        if (nrCols < 0) {
            throw new Error(errWrongArg('nrCols', 'is Negative'));
        }
        if (!isInteger(nrCols)) {
            throw new Error(errWrongArg('nrCols', 'is a NaN'));
        }
        if (!isInteger(colBase)) {
            throw new Error(errWrongArg('colBase', 'is a NaN'));
        }
        return Object.freeze({
            rowBase,
            colBase,
            nrCols,
            nrRows,
            r,
            i,
            coord: (col) => {
                const tb = (col - colBase) * nrRows;
                return row => tb + (row - rowBase)
            },
            //  colOf: (col) => (col - colBase) * nrRows,
            colOfEx: (col) => (col - colBase) * nrRows - rowBase,
            assertComplex: (msg: string) => {
                if (i === undefined) {
                    throw new Error(errMissingIm(msg))
                }
            },
            setCol(col: number, rowStart: number, rowEnd: number, value: number): void {
                const coords = this.colOfEx(col);
                this.r.fill(value, coords + rowStart, coords + rowEnd + 1);
            }
        });
    }
}

export function xerbla(fn: string, idx: number) {
    return ` ** On entry to ${fn}, parameter number ${idx}, had an illegal value`;
}


