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
export type FortranArr = (index: number) => (value?: number) => number;

//2
export function mimicFArray(arr: Float32Array | Float64Array) {
    //Lets make some curry
    let func = function n(offset: number = 0): FortranArr {

        return (index: number) => (value?: number) => {
            const prevV = arr[offset + index];
            if (value !== undefined) {
                arr[offset + index] = value;
            }
            return prevV;
        }
    }
    //   TODO: maybe keep reference later, who knows
    //    func['buffer'] = arr;
    return func;
}

