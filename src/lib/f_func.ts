//intrinsic routines of fortran

const { abs, max } = Math;

const { isInteger } = Number;

const { isArray } = Array;

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

type ArrayElt = { key: string | number, val: any };
export type Complex = { re: number, im: number };
export type fpArray = Float32Array | Float64Array;
export type FortranSetterGetter = (index: number) => (value?: number) => number | Complex;
export type FortranArr = {
    base: number,
    r: fpArray,
    i?: fpArray,
    s: FortranSetterGetter,
    assertComplex: (msg: string) => void | never;
    toArr: () => Complex[] | number[];

};

export function isComplex(a): a is Complex {
    return (
        (a !== null) &&
        (typeof a === 'object') &&
        ('re' in a) &&
        ('im' in a) &&
        (typeof a['re'] === 'number' && typeof a['im'] === 'number')
    );
}
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
            },
            toArr: () => {
                return multiplexer(Array.from(r), i && Array.from(i))((re, im) => i ? { re, im } : re);
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
export function flatten<T>(...rest: (T | T[])[]): T[] {
    let rc: T[] = [];
    for (const itm of rest) {
        if (isArray(itm)) {
            let rc2: T[] = flatten(...itm) as any;
            rc.push(...rc2);
            continue;
        }
        rc.push(itm as any);
    }
    return rc as any;
}

export function muxCmplx(re: number[], im?: number[]): Complex[] {
    return multiplexer(re, im)((re, im) => ({ re, im }))
}

function demuxComplex(...rest: (Complex)[]): { reals: number[], imags?: number[] } {

    const c = flatten(rest);

    const collect: { reals: number[], imags: number[] } = { reals: [], imags: [] };

    c.reduce((prev, v) => {
        if (typeof v === 'number') {
            prev.reals.push(v);
            return prev;
        }
        let re = ('re' in v) ? v.re : v[0];
        let im = ('im' in v) ? v.im : (v['1'] ? v[1] : undefined);
        prev.reals.push(re);
        if (im !== undefined) {
            prev.imags.push(im);
        }
        return prev;
    }, collect);
    //only reals?
    if (collect.reals.length > 0 && collect.imags.length === 0) {
        delete collect.imags;
    }
    return collect;
}

export function fortranArrComplex32(...rest: (Complex | Complex[])[]): (offset?: number) => FortranArr {

    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFArray(new Float32Array(collect.reals), i);
}

export function fortranArrComplex64(...rest: (Complex | Complex[])[]): (offset?: number) => FortranArr {

    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFArray(new Float64Array(collect.reals), i);
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
    [ERROR.ERR_WRONG_ARGUMENT, 'function:[%s], argument [%s] is has invalid value.']
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

export function errWrongArg(arg: string, value: any): string {
    return errorMsg(ERROR.ERR_WRONG_ARGUMENT, arg, JSON.stringify(value));
}


// Matrix
// Matrix
// Matrix
// Matrix

export type FortranMatrixSetterGetter = (colSize: number) => (nrCols: number) => number | Complex;

export type Matrix = {
    readonly rowBase: number,
    readonly colBase: number,
    readonly nrCols: number, // inclusive note!!
    readonly colSize: number,
    readonly r: fpArray, //[(ncols+1)*(nrows+1)]
    readonly i?: fpArray, //imaginary part of matrix [(ncols+1)*(nrows+1)]
    //readonly colOf: (number) => number;
    readonly colOfEx: (number) => number;
    //s: FortranMatrixSetterGetter
    // zap a row with fa valie
    readonly setCol: (col: number, rowStart: number, rowEnd: number, value: number) => void;
    readonly slice: (rowStart: number, rowEnd: number, colStart: number, colEnd: number) => Matrix;
};


function mimicFMatrix(r: fpArray, i?: fpArray) {

    return function c1(lda: number, nrCols: number, rowBase: number = 1, colBase = 1): Matrix {

        // check rows
        if (lda < 0) {
            throw new Error(errWrongArg('nrRows', 'is Negative'));
        }
        if (!isInteger(lda)) {
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
            colSize: lda,
            r,
            i,
            coord: (col) => {
                const tb = (col - colBase) * lda;
                return row => tb + (row - rowBase)
            },
            //  colOf: (col) => (col - colBase) * nrRows,
            colOfEx: (col) => (col - colBase) * lda - rowBase,
            setCol(col: number, rowStart: number, rowEnd: number, value: number): void {
                const coords = this.colOfEx(col);
                this.r.fill(value, coords + rowStart, coords + rowEnd + 1);
                if (this.i) {
                    this.i.fill(value, coords + rowStart, coords + rowEnd + 1);
                }
            },
            slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix {
                const rowSize = (rowEnd - rowStart) + 1;
                const colSize = (colEnd - colStart) + 1;

                let re = new Float64Array(Array.from<number>({ length: rowSize * colSize }));
                let im: Float64Array | undefined;
                if (i) {
                    im = new Float64Array(Array.from<number>({ length: rowSize * colSize }));
                }
                for (let j = colStart; j <= colEnd; j++) {
                    const base = (j - colStart) * rowSize - 1;
                    const coorJ = this.colOfEx(j);
                    for (let i = rowStart; i <= rowEnd; i++) {
                        // console.log(j, coorJ + i, base + i, r[coorJ + i]);
                        re[base + i] = r[coorJ + i];
                        if (im) {
                            im[base + (i - rowStart)] = i[coorJ + i];
                        }
                    }
                }
                return mimicFMatrix(re, im)(rowSize, colSize);
            }
        });
    }
}

export function fortranMatrixComplex32(...rest: (Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {

    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFMatrix(new Float32Array(collect.reals), i);
}

export function fortranMatrixComplex64(...rest: (Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {

    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFMatrix(new Float64Array(collect.reals), i);
}

export function xerbla(fn: string, idx: number) {
    return ` ** On entry to ${fn}, parameter number ${idx}, had an illegal value`;
}

export function lowerChar<T extends string>(c: T): T {
    return String.fromCharCode(c.charCodeAt(0) | 0X20) as any;
}

export const map = iter();
export const each = iter(false);

export function numberPrecision(prec: number = 6) {
    function convert(x: number): number {
        if (isNaN(x)) {
            return NaN;
        }
        return Number.parseFloat(x.toPrecision(prec));
    }

    return arrayrify(convert);
}

function iter<T>(wantMap = true) {
    return function n(xx: T): { (fn: (x: any, idx?: number | string) => any): any | any[] } {
        const fx: ArrayElt[] = coerceToArray(xx) as any;
        return function k(fn: (x: any, idx?: number | string) => any): any | any[] {
            return wantMap ? possibleScalar(fx.map(o => fn(o.val, o.key))) : fx.forEach(o => fn(o.val, o.key));
        };
    }
}

export function arrayrify<T, R>(fn: (x: T, ...rest: any[]) => R) {
    return function n(x: T | T[], ...rest: any[]): R | R[] {
        const fp = Array.isArray(x) ? x : [x];
        const result = fp.map(p => fn(p, ...rest));
        return result.length === 1 ? result[0] : result;
    };
}


function coerceToArray(o: any): { key: string | number, val: any }[] {
    if (o === null || o === undefined) {
        throw new TypeError('Illegal argument excepton: input needs to NOT be "null" or "undefined".');
    }
    if (typeof o === 'number') {
        return [{ key: 0, val: o }] as any;
    }
    if (isArray(o)) {
        return o.map((x, idx) => ({ key: idx, val: x }) as any);
    }
    if (typeof o === 'string') {
        return o.split('').map((x, idx) => ({ key: idx, val: x } as any));
    }
    if (typeof o === 'object') {
        const names = Object.getOwnPropertyNames(o);
        if (names.length === 0) {
            throw new Error('Input argument is an Object with no properties');
        }
        return names.map(name => ({ key: name, val: o[name] })) as any;
    }
    throw new Error('unreachable code');
}


function possibleScalar<T>(x: T[]): T | T[] {
    return x.length === 1 ? x[0] : x;
}


export function multiplexer(...rest: (any | any[])[]) {
    //analyze  
    const analyzed: any[] = [];

    for (let k = 0; k < rest.length; k++) {
        const arg = rest[k];
        // null is special
        if (arg === null) {
            analyzed.push([arg]);
            continue;
        }
        if (['undefined', 'boolean', 'number'].indexOf(typeof arg) >= 0) {
            analyzed.push([arg]);
            continue;
        }
        if (typeof arg === 'string') {
            analyzed.push(arg.split(''));
            continue;
        }
        if (Array.isArray(arg)) {
            analyzed.push(arg);
            continue;
        }
        if (arg instanceof Object) {
            analyzed.push(arg);
            continue;
            //throw new Error('Sorry, looping over properties not yet supported');
        }
        if (arg instanceof Function) {
            throw new Error('Sorry function arguments are not yet supported');
        }
    }//for
    // find the longest array
    const _max = max(...analyzed.map(a => a.length));

    return function k(fn: (...rest: any[]) => any): any | any[] {
        const rc: any[] = [];

        for (let k = 0; k < _max; k++) {
            const result: any[] = [];
            for (let j = 0; j < analyzed.length; j++) {
                const arr: any[] = analyzed[j];
                const idx = k % arr.length;
                result.push(arr[idx]);
            }
            rc.push(fn(...result));
        }
        return possibleScalar(rc);
    };
}

