//intrinsic routines of fortran

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
export type FortranSetterGetter = (index: number) => (value?: number) => number;
export type FortranArr = { base: number, arr: Float32Array | Float64Array, s: FortranSetterGetter };
//2
export function mimicFArray(arr: Float32Array | Float64Array) {
    //Lets make some curry
    let func = function n(startIndex: number = 0): FortranArr {

        return Object.freeze({
            base: startIndex,
            arr,
            s: (index: number) => (value?: number) => {
                const prevV = arr[index - startIndex];
                if (value !== undefined) {
                    arr[index - startIndex] = value;
                }
                return prevV;
            }
        });
    }
    //   TODO: maybe keep reference later, who knows
    //    func['buffer'] = arr;
    return func;
}

