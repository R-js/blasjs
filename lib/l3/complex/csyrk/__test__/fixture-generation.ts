import { rnormOne, rnorm, setSeed } from 'lib-r-math.js';

export default function generateMatrix(fp64: boolean, seed: number, n: number, k: number, complex: boolean, compPat = true) {
    setSeed(seed);
    const step = complex ? 2 : 1
    if (compPat) {
        const result = fp64 ? new Float64Array(n * k * step) : new Float32Array(n * k * step);
        for (let i = 0; i < n * k * step; i += step) {
            const real = rnormOne();
            result[i] = real;
        }
        if (step > 1) {
            for (let i = 1; i < n * k * 2; i += 2) {
                const imag = rnormOne();
                result[i] = imag;
            }
        }
        return result;
    }
    if (fp64) {
        return rnorm(n * k * step);
    }
    return new Float32Array(rnorm(n * k * step));
}
