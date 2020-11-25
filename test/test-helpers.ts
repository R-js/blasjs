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

import { Complex, FortranArr, fpArray, Matrix, errWrongArg } from '../src/lib/f_func';

const { max, min } = Math;

const { isInteger } = Number;

const { isArray } = Array;

type ArrayElt = { key: string | number; val: any };

// type guard
export function isComplex(a: any): a is Complex {
    return (
        a !== null &&
        typeof a === 'object' &&
        're' in a &&
        'im' in a &&
        typeof a['re'] === 'number' &&
        typeof a['im'] === 'number'
    );
}

export function approximatelyWithPrec(prec: number): (act: number | Complex, exp: number) => void {
    return (act: number | Complex, exp: number) => approximately(act, exp, prec);
}

export function approximately(act: number | Complex, exp: number | Complex, preca = 1e-9) {
    const prec = -Math.log10(preca);
    if (isComplex(act)) {
        if (isComplex(exp)) {
            expect(act.re).toBeCloseTo(exp.re, prec); //, 'real part numbers are NOT close');
            expect(act.im).toBeCloseTo(exp.im, prec);
            //assert.approximately(act.im, exp.im, prec, 'numbers are NOT close');
            return;
        }
        throw new Error('cannot compare complex type with non complex type');
    }

    switch (true) {
        case isNaN(act):
            expect(exp).toBeNaN();
            break;
        case isFinite(act):
            //console.log(act, exp, prec);
            expect(act).toBeCloseTo(<number>exp, prec);
            //assert.approximately(act, <number>exp, prec, 'numbers are NOT close');
            break;
        case !isFinite(act):
            expect(act).toBe(exp);
            //assert.equal(act, exp);
            break;
        default:
            throw new Error(`Icompatible values ${act}, ${exp}`);
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

export function fortranArrComplex32(
    ...rest: (number | number[] | Complex | Complex[])[]
): (offset?: number) => FortranArr {
    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFArray(new Float32Array(collect.reals), i);
}

export function fortranArrComplex64(
    ...rest: (number | number[] | Complex | Complex[])[]
): (offset?: number) => FortranArr {
    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFArray(new Float64Array(collect.reals), i);
}

function demuxComplex(...rest: (Complex | number)[]): { reals: number[]; imags?: number[] } {
    const c = flatten(rest);

    const collect: { reals: number[]; imags?: number[] } = { reals: [], imags: [] };

    c.reduce((prev, v) => {
        if (typeof v === 'number') {
            prev.reals.push(v);
            return prev;
        }
        const re = 're' in v ? v.re : undefined;
        const im = 'im' in v ? v.im : undefined;
        if (re !== undefined) {
            prev.reals.push(re);
        }
        if (im !== undefined && prev.imags) {
            prev.imags.push(im);
        }
        return prev;
    }, collect);
    //only reals?
    if (collect.reals.length > 0) {
        if (collect.imags && collect.imags.length === 0) {
            delete collect.imags;
        }
    }
    return collect;
}

export function flatten<T>(...rest: (T | T[])[]): T[] {
    const rc: T[] = [];
    for (const itm of rest) {
        if (isArray(itm)) {
            const rc2: T[] = flatten(...itm) as any;
            rc.push(...rc2);
            continue;
        }
        rc.push(itm as any);
    }
    return rc as any;
}

//2
export function mimicFArray(r: fpArray, i?: fpArray) {
    //Lets make some curry
    const func = function n(startIndex = 1): FortranArr {
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
                return multiplexer(Array.from(r), i && Array.from(i))((re, im) => (i ? { re, im } : re));
            },
        });
    };
    //   TODO: maybe keep reference later, who knows
    //    func['buffer'] = arr;
    return func;
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
    } //for
    // find the longest array
    const _max = max(...analyzed.map((a) => a.length));

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

export function complex(re = 0, im = 0): Complex {
    return { re, im };
}

export function muxCmplx(re: number[], im?: number[]): Complex[] {
    return multiplexer(re, im)((re, im) => ({ re, im }));
}

export function fortranMatrixComplex32(
    ...rest: (number | number[] | Complex | Complex[])[]
): (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {
    const collect = demuxComplex(rest as any);
    let i: Float32Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float32Array(collect.imags);
    }
    return mimicFMatrix(new Float32Array(collect.reals), i);
}

export function fortranMatrixComplex64(
    ...rest: (number | number[] | Complex | Complex[])[]
): (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix {
    const collect = demuxComplex(rest as any);
    let i: Float64Array | undefined;
    if (collect.imags !== undefined) {
        i = new Float64Array(collect.imags);
    }
    return mimicFMatrix(new Float64Array(collect.reals), i);
}

export function mimicFMatrix(r: fpArray, i?: fpArray) {
    return function c1(lda: number, nrCols: number, rowBase = 1, colBase = 1): Readonly<Matrix> {
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
            coord(col: number) {
                const tb = (col - colBase) * lda;
                return (row: number) => tb + (row - rowBase);
            },
            //  colOf: (col) => (col - colBase) * nrRows,
            colOfEx: (col) => (col - colBase) * lda - rowBase,
            setCol(col: number, rowStart: number, rowEnd: number, value = 0): void {
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
                        rc[coor + i] =
                            this.i === undefined ? this.r[coor + i] : { re: this.r[coor + i], im: this.i[coor + i] };
                    }
                }
                return rc;
            },
            slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix {
                const _nrRows = rowEnd - rowStart + 1;
                const _nrCols = colEnd - colStart + 1;

                const re = new Float64Array(_nrRows * _nrCols);
                let im: Float64Array | undefined;
                //
                if (this.i) {
                    im = new Float64Array(_nrRows * _nrCols);
                }
                //
                for (let j = colStart; j <= colEnd; j++) {
                    //  x x | x x
                    const targetBase = (j - colStart) * _nrRows;
                    const coorJ = this.colOfEx(j);
                    for (let i = rowStart; i <= rowEnd; i++) {
                        // console.log(j, coorJ + i, base + i, r[coorJ + i]);
                        re[targetBase + (i - rowStart)] = r[coorJ + i];
                        if (this.i && im) {
                            im[targetBase + (i - rowStart)] = this.i[coorJ + i];
                        }
                    }
                }
                return mimicFMatrix(re, im)(_nrRows, _nrCols);
            },
            slice_used_for_test(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix {
                const _nrRows = rowEnd - rowStart + 1;
                const _nrCols = colEnd - colStart + 1;

                const re = new Float64Array(_nrRows * _nrCols);
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
                const ic = this.i ? this.i.slice() : undefined;
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
                const ic = this.i ? this.i.slice() : undefined;
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
                const reC = this.r.slice();
                return mimicFMatrix(reC)(this.nrRows, this.nrCols);
            },
            imaginary() {
                let imC: fpArray;
                if (!this.i) {
                    imC = new Float64Array(this.nrRows * this.nrCols);
                    imC.fill(0);
                } else {
                    imC = this.i.slice(); //copy
                }
                return mimicFMatrix(imC)(this.nrRows, this.nrCols);
            },
        });
    };
}

function bandifyMatrix(uplo: 'u' | 'l', k: number, A: Matrix): Matrix {
    const rowSize = k + 1;
    const colSize = A.nrCols;
    const createArr = A.r instanceof Float64Array ? Float64Array : Float32Array;

    const re = new createArr(rowSize * colSize);
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
    const packedSize = (-(k + 1) * k) / 2 + A.nrCols * (k + 1);
    const colSize = A.nrCols;
    const createArr = A.r instanceof Float64Array ? Float64Array : Float32Array;

    const re = new createArr(packedSize);
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

export function coerceToArray(o: any): { key: string | number; val: any }[] {
    if (o === null || o === undefined) {
        return []; //throw new TypeError('Illegal argument excepton: input needs to NOT be "null" or "undefined".');
    }
    if (typeof o === 'number') {
        return [{ key: 0, val: o }] as any;
    }
    if (isArray(o)) {
        return o.map((x, idx) => ({ key: idx, val: x } as any));
    }
    if (typeof o === 'string') {
        return o.split('').map((x, idx) => ({ key: idx, val: x } as any));
    }
    if (typeof o === 'object') {
        const names = Object.getOwnPropertyNames(o);
        if (names.length === 0) {
            return []; //throw new Error('Input argument is an Object with no properties');
        }
        return names.map((name) => ({ key: name, val: o[name] })) as any;
    }
    throw new Error('unreachable code');
}

function iter<T>(wantMap = true) {
    return function n(xx: T): { (fn: (x: any, idx?: number | string) => any): any | any[] } {
        const fx: ArrayElt[] = coerceToArray(xx) as any;
        return function k(fn: (x: any, idx?: number | string) => any): any | any[] {
            return wantMap ? fx.map((o) => fn(o.val, o.key)) : fx.forEach((o) => fn(o.val, o.key));
        };
    };
}

export const map = iter();
export const each = iter(false);

export function numberPrecision(prec = 6) {
    const runner = arrayrify(convert);
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
    return runner;
}

export function arrayrify<T, R>(fn: (x: T, ...rest: any[]) => R) {
    return function n(x: T | T[], ...rest: any[]): R | R[] {
        const fp = Array.isArray(x) ? x : [x];
        const result = fp.map((p) => fn(p, ...rest));
        return result.length === 1 ? result[0] : result.length === 0 ? (undefined as any) : result;
    };
}
