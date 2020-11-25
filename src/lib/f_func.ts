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

export function isZero(v: Complex): boolean {
    return v.im === 0 && v.re === 0;
}

export function isZeroE(re: number, im: number): boolean {
    return re === 0 && im === 0;
}

export function isOne(v: Complex): boolean {
    return v.re === 1 && v.im === 0;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any

export type Complex = { re: number; im: number };
export type fpArray = Float32Array | Float64Array;
export type FortranSetterGetter = (index: number) => (re?: number, im?: number) => number | Complex;
export type FortranArr = {
    base: number;
    r: fpArray;
    i?: fpArray;
    s: FortranSetterGetter;
    //assertComplex: (msg: string) => void | never;
    toArr: () => Complex[] | number[];
};

export type FortranArrEComplex = {
    base: number;
    r: fpArray;
    i: fpArray;
    s: FortranSetterGetter;
    //assertComplex: (msg: string) => void | never;
    toArr: () => Complex[] | number[];
};

export enum ERROR {
    ERR_MISSING_IMAGINARY = 1,
    ERR_MISSING_REAL = 2,
    ERR_WRONG_ARGUMENT = 3,
}

// private
const errors = new Map<ERROR, string>([
    [ERROR.ERR_MISSING_IMAGINARY, 'Missing imaginary part for %s'],
    [ERROR.ERR_MISSING_REAL, 'Missing real part for %s'],
    [ERROR.ERR_WRONG_ARGUMENT, 'function:[%s], argument [%s] is has invalid value.'],
]);

const ERROR_UKNOWN = 'Unkown Error code used! [%s]';

export function errorMsg(errNo: ERROR, ...fmt: string[]): string {
    const msg = errors.get(errNo) || ERROR_UKNOWN;
    return fmt.reduce((p, v) => p.replace('%s', v), msg);
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
    colOfEx(n: number): number;
    //s: FortranMatrixSetterGetter
    // zap a row with fa valie
    coord(col: number): (a: number) => number;
    setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
    slice_used_for_test(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
    slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
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

export function lowerChar<T extends string>(c?: T): T {
    if (typeof c !== 'string') {
        return '' as T;
    }
    const cc = c.charCodeAt(0);
    if ((cc >= 65 && cc <= 90) || (cc >= 97 && cc <= 122)) {
        return String.fromCharCode(cc | 0x20) as any;
    }
    return c;
}

export function possibleScalar<T>(x: T[]): T | T[] {
    return x.length === 1 ? x[0] : x;
}

export function mul_rxc(re1: number, im1: number, c: Complex): Complex {
    const re = re1 * c.re - im1 * c.im;
    const im = re1 * c.im + im1 * c.re;
    return { re, im };
}

export function mul_rxr(re1: number, im1: number, re2: number, im2: number): Complex {
    const re = re1 * re2 - im1 * im2;
    const im = re1 * im2 + im1 * re2;
    return { re, im };
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
    const re = (re1 * re2 + im1 * im2) / norm;
    const im = (-re1 * im2 + im1 * re2) / norm;
    return { re, im };
}
