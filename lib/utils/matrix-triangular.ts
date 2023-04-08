export function upperPack<T extends Float32Array | Float64Array>(
  a: T,
  rank: number,
  isComplex: boolean
): T {
  const step = isComplex ? 2 : 1;
  const numElts = ((rank * (rank + 1)) / 2) * step;
  const b = new (a.constructor as new (n: number) => T)(numElts);
  for (
    let i = 0, colBase = 0, colBaseSrc = 0;
    i < rank;
    i++, colBase += i * step, colBaseSrc += rank * step
  ) {
    b.set(a.slice(colBaseSrc, colBaseSrc + (i + 1) * step), colBase);
  }
  return b;
}

export function upper<T extends Float32Array | Float64Array>(
  a: T,
  rank: number,
  isComplex: boolean
): T {
  const step = isComplex ? 2 : 1;
  for (let i = 0; i < rank; i++) {
    a.fill(0, step * (rank * i + (i + 1)), rank * (i + 1) * step);
  }
  return a;
}

export function lower<T extends Float32Array | Float64Array>(
  a: T,
  rank: number,
  isComplex: boolean
): T {
  const step = isComplex ? 2 : 1;
  for (let i = 0; i < rank; i++) {
    a.fill(0, rank * i * step, step * (rank * i + i));
  }
  return a;
}

export function lowerPack<T extends Float32Array | Float64Array>(
  a: T,
  rank: number,
  isComplex: boolean
): T {
  const step = isComplex ? 2 : 1;
  const numElts = ((rank * (rank + 1)) / 2) * step;
  const b = new (a.constructor as new (n: number) => T)(numElts);
  for (
    let i = 0, colBase = 0, colBaseSrc = 0;
    i < rank;
    i++, colBase += (rank - (i - 1)) * step, colBaseSrc += (rank + 1) * step
  ) {
    b.set(a.slice(colBaseSrc, colBaseSrc + (rank - i) * step), colBase);
  }
  return b;
}

export function transposeStorage<T extends Float32Array | Float64Array>(
  a: T,
  n: number,
  k: number,
  isComplex: boolean
): T {
  const step = isComplex ? 2 : 1;
  const result = new (a.constructor as new (n: number) => T)(n * k * step);
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < n; j++) {
      const real = a[(i * n + j) * step];
      result[(j * k + i) * step] = real;
      if (isComplex) {
        const imag = a[(i * n + j) * step + 1];
        result[(j * k + i) * step + 1] = imag;
      }
    }
  }
  return result;
}
