import { isFloat64Array } from 'util/types';

export function copyMatricesToWasmMemory(forceHiFP: boolean, ...matrices: (ArrayLike<number> | Float64Array | Float32Array)[]) {
    let byteLength = 0;
    // count their length to get total mem usage
    const hiFP = forceHiFP ? true : matrices.some(isFloat64Array);
    for (let i = 0; i < matrices.length; i++) {
        const m = matrices[i];
        if (m.length === 0){
            continue;
        }
        const factor = hiFP ? Float64Array.BYTES_PER_ELEMENT : Float32Array.BYTES_PER_ELEMENT;
        byteLength += m.length * factor;
    }
    const initial = Math.ceil(byteLength / 65536);
    const mem = new WebAssembly.Memory({ initial });
    byteLength = 0;
    for (let i = 0; i < matrices.length; i++) {
        const m = matrices[i];
        if (m.length === 0){
            continue;
        }
        const factor = hiFP ? Float64Array.BYTES_PER_ELEMENT : Float32Array.BYTES_PER_ELEMENT;
        const storageType = hiFP ? Float64Array : Float32Array;
        const arr = new storageType(mem.buffer, byteLength, m.length);
        arr.set(m);
        byteLength += m.length * factor;
    }
    return { mem, byteLength, hiFP };
}

export function mapWasmMemoryToMatrices(forceHiFP: boolean, storage: WebAssembly.Memory, ...dimensions: number[]) {
    // count their length to get total mem usage
    let byteLength = dimensions.reduce((len, n) => {
        return len + n * (forceHiFP ? Float64Array.BYTES_PER_ELEMENT : Float32Array.BYTES_PER_ELEMENT);
    }, 0);

    if (byteLength > storage.buffer.byteLength){
        throw new Error(`You are asking for more bytes then there is available in storage ${byteLength} > ${storage.buffer.byteLength} `);
    }
    byteLength = 0;
    const matrices: (Float64Array|Float32Array)[] = Array.from({length: dimensions.length }); 
    for (let i = 0; i < dimensions.length; i++){
        const dim = dimensions[i];
        const StorageType = forceHiFP ? Float64Array : Float32Array;
        matrices[i] = new StorageType(storage.buffer, byteLength, dim);
        byteLength += dim*(forceHiFP ? Float64Array.BYTES_PER_ELEMENT : Float32Array.BYTES_PER_ELEMENT);
    }
    return matrices;
}


