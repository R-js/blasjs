/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//intrinsic routines of fortran

const { max, min } = Math;

const { isInteger } = Number;

const { isArray } = Array;

export function isZero(v: Complex) {
    return v.im === 0 && v.re === 0;
}

export function isZeroE(re: number, im: number) {
    return re === 0 && im === 0;
}

export function isOne(v: Complex) {
    return v.re === 1 && v.im === 0;
}
/*
export function sign(a: number, b?: number): number {
    if (b === undefined) {
        return a;
    }
    const rc = Math.abs(a);
    return b >= 0 ? rc : -rc;
}
*/
//make JS arrays work like fortran ones
//1. export type
//2. export array wrapper factory

//1

type ArrayElt = { key: string | number, val: any };
export type Complex = { re: number, im: number };
export type fpArray = Float32Array | Float64Array;
export type FortranSetterGetter = (index: number) => (re?: number, im?: number) => number | Complex;
export type FortranArr = {
    base: number,
    r: fpArray,
    i?: fpArray,
    s: FortranSetterGetter,
    //assertComplex: (msg: string) => void | never;
    toArr: () => Complex[] | number[];

};

export type FortranArrEComplex = {
    base: number,
    r: fpArray,
    i: fpArray,
    s: FortranSetterGetter,
    //assertComplex: (msg: string) => void | never;
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
            s: (index: number) => (re?: number, im?: number): number | Complex => {

                if (index - startIndex >= r.length) {
                    throw new Error('Index out of bounds');
                }
                const oldRe = r[index - startIndex];
                const oldIm = i ? i[index - startIndex] : 0;
                const onlyGet = re === undefined && im === undefined;
                const onlyReal = r !== undefined && i === undefined;

                if (onlyGet) {
                    if (onlyReal) {
                        return oldRe;
                    }
                    return complex(oldRe, oldIm);
                }

                if (onlyReal && im !== undefined) {
                    throw new Error('You specified a complex number for a real array');
                }

                if (i !== undefined) {
                    i[index - startIndex] = im || 0;
                }


                r[index - startIndex] = re || 0;

                return { re: oldRe, im: oldIm || 0 };
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

function demuxComplex(...rest: (Complex | number)[]): { reals: number[], imags?: number[] } {

    const c = flatten(rest);

    const collect: { reals: number[], imags: number[] } = { reals: [], imags: [] };

    c.reduce((prev, v) => {
        if (typeof v === 'number') {
            prev.reals.push(v);
            return prev;
        }
        let re = ('re' in v) ? v.re : undefined;
        let im = ('im' in v) ? v.im : undefined;
        if (re !== undefined) {
            prev.reals.push(re);
        }
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

export function fortranArrComplex32(...rest: (number | number[] | Complex | Complex[])[]): (offset?: number) => FortranArr {

    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFArray(new Float32Array(collect.reals), i);
}

export function fortranArrComplex64(...rest: (number | number[] | Complex | Complex[])[]): (offset?: number) => FortranArr {

    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFArray(new Float64Array(collect.reals), i);
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

    let msg = errors.get(errNo) || ERROR_UKNOWN;
    /*if (!msg) {
        msg = ERROR_UKNOWN;
    }*/
    /*if (fmt.length === 0) {
        return msg;
    }*/
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

export type MatrixType = 'n' | 'b' | 'bu' | 'bl' | 'pl' | 'pu'; //untyped, normal, bandified, bandified-upper, bandified-lower, packed-lower, packed-upper 

export interface Matrix {
    readonly rowBase: number;
    readonly colBase: number;
    readonly nrCols: number; // inclusive note!!
    readonly nrRows: number;
    readonly r: fpArray; //[(ncols+1)*(nrows+1)]
    readonly i?: fpArray; //imaginary part of matrix [(ncols+1)*(nrows+1)]
    //readonly colOf: (number) => number;
    readonly colOfEx: (number) => number;
    //s: FortranMatrixSetterGetter
    // zap a row with fa valie
    coord(col): (row) => number;
    setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
    slice_used_for_test(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
    setLower(value?: number): Matrix;
    setUpper(value?: number): Matrix;
    upperBand(value: number): Matrix;
    lowerBand(value: number): Matrix;
    packedUpper(value?: number): FortranArr;
    packedLower(value?: number): FortranArr;
    real(): Matrix;
    imaginary(): Matrix;
    toArr(): Complex[] | number[];
}

export interface MatrixEComplex extends Matrix {
    readonly i: fpArray;
}


function bandifyMatrix(uplo: 'u' | 'l', k: number, A: Matrix): Matrix {
    const rowSize = (k + 1);
    const colSize = A.nrCols;
    const createArr = A.r instanceof Float64Array ? Float64Array : Float32Array;

    let re = new createArr(rowSize * colSize);
    let im: fpArray | undefined;
    if (A.i) {
        im = new createArr(rowSize * colSize);
    }
    for (let j = 1; j <= colSize; j++) {
        const m = uplo === 'u' ? k + 1 - j : 1 - j;
        const coorJ = A.colOfEx(j);
        const base = (j - 1) * rowSize - 1;
        const start = uplo === 'u' ? max(1, j - k) : j;
        const stop = uplo === 'u' ? j : min(colSize, j + k);
        for (let i = start; i <= stop; i++) {
            re[base + m + i] = A.r[coorJ + i];
            if (im && A.i) {
                im[base + m + i] = A.i[coorJ + i];
            }
        }
    }
    return mimicFMatrix(re, im)(rowSize, colSize);
}

function packedBandi_fied_Matrix(uplo: 'u' | 'l', k: number, A: Matrix): FortranArr {

    const packedSize = -(k + 1) * k / 2 + A.nrCols * (k + 1);
    const colSize = A.nrCols;
    const createArr = A.r instanceof Float64Array ? Float64Array : Float32Array;

    let re = new createArr(packedSize);
    let im: fpArray | undefined;
    if (A.i) {
        im = new createArr(packedSize);
    }
    let cursor = 0;
    for (let j = 1; j <= colSize; j++) {
        //const m = uplo === 'u' ? k + 1 - j : 1 - j;
        const coorJ = A.colOfEx(j);
        //const base = (j - 1) * rowSize - 1;
        const start = uplo === 'u' ? max(1, j - k) : j;
        const stop = uplo === 'u' ? j : min(colSize, j + k);
        for (let i = start; i <= stop; i++) {
            re[cursor] = A.r[coorJ + i];
            //re[base + m + i] = A.r[coorJ + i];
            if (im && A.i) {
                im[cursor] = A.i[coorJ + i];
                //im[base + m + i] = A.i[coorJ + i];
            }
            cursor++;
        }
    }
    return mimicFArray(re, im)();
}


export function mimicFMatrix(r: fpArray, i?: fpArray) {

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


        return Object.freeze<Matrix>({
            rowBase,
            colBase,
            nrCols,
            nrRows: lda,
            r,
            i,
            coord(col) {
                const tb = (col - colBase) * lda;
                return row => tb + (row - rowBase)
            },
            //  colOf: (col) => (col - colBase) * nrRows,
            colOfEx: (col) => (col - colBase) * lda - rowBase,
            setCol(col: number, rowStart: number, rowEnd: number, value: number = 0): void {
                // dont use colOf (avoid extra function)
                const coords = (col - this.colBase) * lda - rowBase;
                this.r.fill(value, coords + rowStart, coords + rowEnd + 1);
                if (this.i) {
                    this.i.fill(value, coords + rowStart, coords + rowEnd + 1);
                }
            },
            upperBand(k: number = lda - 1) {
                return bandifyMatrix('u', k, this);
            },
            lowerBand(k: number = lda - 1) {
                return bandifyMatrix('l', k, this);
            },
            packedUpper(k: number = lda - 1) {
                return packedBandi_fied_Matrix('u', k, this);
            },
            packedLower(k: number = lda - 1) {
                return packedBandi_fied_Matrix('l', k, this);
            },
            toArr(): Complex[] | number[] {
                const rc: any[] = new Array(lda * nrCols).fill(0);
                for (let j = 1; j <= nrCols; j++) {
                    const coor = (j - 1) * this.nrRows - 1;
                    for (let i = 1; i <= lda; i++) {
                        rc[coor + i] = (this.i === undefined) ? this.r[coor + i] : { re: this.r[coor + i], im: this.i[coor + i] };
                    }
                }
                return rc;
            },
            slice_used_for_test(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix {
                const _nrRows = (rowEnd - rowStart) + 1;
                const _nrCols = (colEnd - colStart) + 1;

                let re = new Float64Array(_nrRows * _nrCols);
                let im: Float64Array | undefined;
                //
                if (this.i) {
                    im = new Float64Array(_nrRows * _nrCols);
                }
                //
                for (let j = colStart; j <= colEnd; j++) {
                    const base = (j - colStart) * _nrRows - 1;
                    const coorJ = this.colOfEx(j);
                    for (let i = rowStart; i <= rowEnd; i++) {
                        // console.log(j, coorJ + i, base + i, r[coorJ + i]);
                        re[base + i] = r[coorJ + i];
                        if (this.i && im) {
                            im[base + (i - rowStart)] = this.i[coorJ + i];
                        }
                    }
                }
                return mimicFMatrix(re, im)(_nrRows, _nrCols);
            },
            setLower(value = 0) {
                const rc = this.r.slice();
                let ic = this.i ? this.i.slice() : undefined;
                for (let y = 1; y < nrCols; y++) {
                    const coor = this.colOfEx(y);
                    rc.fill(value, coor + y + 1, coor + lda + 1);
                    if (ic) {
                        ic.fill(value, coor + y + 1, coor + lda + 1);
                    }
                }
                return mimicFMatrix(rc, ic)(this.nrRows, this.nrCols, this.rowBase, this.colBase);
            },
            setUpper(value = 0) {
                const rc = this.r.slice();
                let ic = this.i ? this.i.slice() : undefined;
                for (let y = 2; y <= nrCols; y++) {
                    const coor = this.colOfEx(y);
                    rc.fill(value, coor + 1, coor + min(y, lda));
                    if (ic) {
                        ic.fill(value, coor + 1, coor + min(y, lda));
                    }
                }
                return mimicFMatrix(rc, ic)(this.nrRows, this.nrCols, this.rowBase, this.colBase);
            },
            real() {
                let reC: fpArray;
                reC = this.r.slice();
                return mimicFMatrix(reC)(this.nrRows, this.nrCols);
            },
            imaginary() {
                let imC: fpArray;
                if (!this.i) {
                    imC = new Float64Array(this.nrRows * this.nrCols);
                    imC.fill(0);
                }
                else {
                    imC = this.i.slice(); //copy
                }
                return mimicFMatrix(imC)(this.nrRows, this.nrCols);
            }
        });
    }
}

export function real(data?: number | number[] | Complex | Complex[]): number[] | number {

    if (typeof data === 'object' && 're' in data) {
        return data.re;
    }

    if (typeof data === 'number') {
        return data;
    }

    if (data instanceof Array) {
        //console.log('its an array');
        const copy: number[] = new Array<number>(data.length);
        for (let i = 0; i < data.length; i++) {
            copy[i] = real(data[i]) as any;
        }
        return copy;
    }
    return data as any;
}

export function imaginary(data?: number | number[] | Complex | Complex[]): number[] | number {

    if (typeof data === 'object' && 'im' in data) {
        //console.log(data.im);
        return data.im;
    }

    if (typeof data === 'number') {
        return 0;
    }

    if (data instanceof Array) {
        //console.log('its an array');
        const copy: number[] = new Array<number>(data.length);
        for (let i = 0; i < data.length; i++) {
            copy[i] = imaginary(data[i]) as any;
        }
        return copy;
    }
    return data as any;
}

export function fortranMatrixComplex32(...rest: (number | number[] | Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {

    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFMatrix(new Float32Array(collect.reals), i);
}

export function fortranMatrixComplex64(...rest: (number | number[] | Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {

    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFMatrix(new Float64Array(collect.reals), i);
}

export function lowerChar<T extends string>(c?: T): T {
    if (typeof c !== 'string') {
        return '' as any;
    }
    const cc = c.charCodeAt(0);
    if ((cc >= 65 && cc <= 90) || (cc >= 97 && cc <= 122)) {
        return String.fromCharCode(cc | 0x20) as any;
    }
    return c;
}

export const map = iter();
export const each = iter(false);

export function numberPrecision(prec: number = 6) {

    let runner: Function;
    function convert(x?: number | any): number {
        // try to loop over the object
        if (typeof x === 'object' && x !== null) {
            for (const key in x) {
                //recursion
                x[key] = runner(x[key]);
            }
            return x;
        }
        // this is a number!!
        if (typeof x === 'number') {
            return Number.parseFloat(x.toPrecision(prec));
        }
        //dont change the object, whatever it is
        return x;
    }
    runner = arrayrify(convert);
    return runner;
}

function iter<T>(wantMap = true) {
    return function n(xx: T): { (fn: (x: any, idx?: number | string) => any): any | any[] } {
        const fx: ArrayElt[] = coerceToArray(xx) as any;
        return function k(fn: (x: any, idx?: number | string) => any): any | any[] {
            return wantMap ? fx.map(o => fn(o.val, o.key)) : fx.forEach(o => fn(o.val, o.key));
        };
    }
}

export function arrayrify<T, R>(fn: (x: T, ...rest: any[]) => R) {
    return function n(x: T | T[], ...rest: any[]): R | R[] {
        const fp = Array.isArray(x) ? x : [x];
        const result = fp.map(p => fn(p, ...rest));
        return result.length === 1 ? result[0] : result.length === 0 ? undefined as any : result;
    };
}

export function coerceToArray(o: any): { key: string | number, val: any }[] {
    if (o === null || o === undefined) {
        return []; //throw new TypeError('Illegal argument excepton: input needs to NOT be "null" or "undefined".');
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
            return []; //throw new Error('Input argument is an Object with no properties');
        }
        return names.map(name => ({ key: name, val: o[name] })) as any;
    }
    throw new Error('unreachable code');
}


export function possibleScalar<T>(x: T[]): T | T[] {
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
        if (arg instanceof Function) {
            throw new Error('Sorry function arguments are not yet supported');
        }
        //if (arg instanceof Object) {
        analyzed.push([arg]);
        continue;
        //throw new Error('Sorry, looping over properties not yet supported');
        //}

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
        return rc; // possibleScalar(rc);
    };
}
/*
export function render(
    a: string,
    aix: string,
    b: string,
    bix: string,
    conjA = false,
    conjB = false): string {
 
    const re = `= ${a}.r[${aix} - ${a}.base]*${b}.r[${bix}-${b}.base] - ${a}.i[${aix} - ${a}.base]*${b}.i[${bix}-${b}.base];`;
    const im = `= ${a}.r[${aix} - ${a}.base]*${b}.i[${bix}-${b}.base] + ${a}.i[${aix} - ${a}.base]*${b}.r[${bix}-${b}.base];`;
    return `${re}\n${im}\n`;
}
*/

export function mul_rxc(re1: number, im1: number, c: Complex): Complex {
    const re = re1 * c.re - im1 * c.im;
    const im = re1 * c.im + im1 * c.re;
    return { re, im };
}

export function mul_rxr(re1: number, im1: number, re2: number, im2: number): Complex {
    const re = re1 * re2 - im1 * im2;
    const im = re1 * im2 + im1 * re2;
    return { re, im }
}

export function mul_cxr(c: Complex, re1: number, im1: number): Complex {
    const re = c.re * re1 - c.im * im1;
    const im = c.re * im1 + c.im * re1;
    return { re, im };
}

export function mul_cxc(c1: Complex, c2: Complex): Complex {
    const re = c1.re * c2.re - c1.im * c2.im;
    const im = c1.re * c2.im + c1.im * c2.re;
    return { re, im };
}

export function div_rxr(re1: number, im1: number, re2: number, im2: number): Complex {
    // (a+ib)/c+id) = (ac+bd)/(c*c +d*d) + i(-ad+bc)/(c*c+d*d)
    const norm = re2 * re2 + im2 * im2;
    let re = (re1 * re2 + im1 * im2) / norm;
    let im = (-re1 * im2 + im1 * re2) / norm;
    return { re, im };
}
