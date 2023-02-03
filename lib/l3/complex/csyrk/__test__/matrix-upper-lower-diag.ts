import { isFloat32Array } from "util/types";

// matrices are "column major"
export function upperPack<T extends Float32Array|Float64Array>(a: T, rank: number, complexElts: boolean): T { 
    const step = complexElts ? 2 : 1;
    const numElts = rank*(rank+1)/2*step;
    const b = isFloat32Array(a) ? new Float32Array(numElts) : new Float64Array(numElts);
    for (let i = 0, colBase = 0, colBaseSrc = 0; i < rank; i++, colBase += i*step, colBaseSrc += rank*step){
        b.set(a.slice(colBaseSrc, (colBaseSrc + (i+1)*step)), colBase);
    }
    return b as T;
}

// matrices are "column major", in place mutation
export function upper<T extends Float32Array|Float64Array>(a: T, rank: number, complexElts: boolean) { 
    const step = complexElts ? 2 : 1;
    for (let i = 0; i < rank; i++){
       a.fill(0, rank*i*step+(i+1)*step, rank*(i+1)*step);
    }
    return a;
}

/*
x x x x
x x x x
x x x x
x x x x
*/