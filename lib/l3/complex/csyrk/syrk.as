export function zsyrk(
  upper: bool,
  transPose: bool,
  n: u32,
  k: u32,
  alphaRe: f64,
  alphaIm: f64,
  betaRe: f64,
  betaIm: f64,
  isPacked: bool
): void {
  const alphaIsZero: bool = alphaRe == 0 && alphaIm == 0
  const betaIsOne: bool = betaRe == 1 && betaIm == 0
  const betaIsZero: bool = betaRe == 0 && betaIm == 0;
  const NN: u32 = n * 2;
  const KK: u32 = k * 2;

  if (n == 0 || ((alphaIsZero || k == 0) && betaIsOne)) return;

  // f64 buffer
  const BYTES_PER_ELEMENT: u32 = 8;

  const ABase: u32 = BYTES_PER_ELEMENT * (isPacked ? n * (n + 1) : n * n * 2);
  // let logIdx: u32 = ABase + BYTES_PER_ELEMENT * (n * k * 2);

  let packCursor: u32 = 0;

  if (alphaIsZero) {
    for (let j: u32 = 0, colBase: u32 = 0; j < n; j++, colBase += NN) { // j = column index
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1));
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN);
      const len: u32 = stop - start;
      if (betaIsZero) {
        if (isPacked) {
          memory.fill(packCursor, 0, len);
          packCursor += len;
        }
        else {
          memory.fill(start, 0, len);
        }
      } else {
        for (let i: u32 = start; i < stop; i += 2 * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {
          const idx: u32 = isPacked ? packCursor : i;
          const cRe: f64 = load<f64>(idx);
          const cIm: f64 = load<f64>(idx + BYTES_PER_ELEMENT);
          const re: f64 = betaRe * cRe - betaIm * cIm;
          const im: f64 = betaRe * cIm + betaIm * cRe;
          store<f64>(idx, re);
          store<f64>(idx + BYTES_PER_ELEMENT, im);
        }
      }
    }
    return;
  }

  if (!transPose) {
    // Form  C := alpha*A*A**T + beta*C.
    // j is the jth column of matrix C
    for (let j: u32 = 0, colBase: u32 = 0; j < n; j++, colBase += NN) {
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1));
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN);
      const len: u32 = stop - start;

      // at this point you want to do  C += alpha*A*A**T
      // A has "k" columns and "n" rows
      for (let i: u32 = start; i < stop; i += 2 * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {
        // column-major
        // because of transpose symmetry we only loop over [1..k]

        let matrixAStart: u32 = i % (NN * BYTES_PER_ELEMENT);
        let matrixATStart: u32 = (j << 1) * BYTES_PER_ELEMENT;

        // the column in C is the column in A transpose with is the row in A

        let tempRe: f64 = 0;
        let tempIm: f64 = 0;

        for (let l: u32 = 0; l < k; l++, matrixAStart += NN * BYTES_PER_ELEMENT, matrixATStart += NN * BYTES_PER_ELEMENT) {
          const aIdx: u32 = ABase + matrixAStart;
          const aTidx: u32 = ABase + matrixATStart;

          const aRe: f64 = f64.load(aIdx);
          const aIm: f64 = f64.load(aIdx + BYTES_PER_ELEMENT);

          const aTRe: f64 = f64.load(aTidx);
          const aTIm: f64 = f64.load(aTidx + BYTES_PER_ELEMENT);

          tempRe += aRe * aTRe - aIm * aTIm;
          tempIm += aRe * aTIm + aIm * aTRe;
        }
        const idx: u32 = isPacked ? packCursor : i;
        let re: f64 = alphaRe * tempRe - alphaIm * tempIm;
        let im: f64 = alphaRe * tempIm + alphaIm * tempRe;
        if (!betaIsZero) {
          let cRe = f64.load(idx);
          let cIm = f64.load(idx + BYTES_PER_ELEMENT)
          if (!betaIsOne) {
            const bre = betaRe * cRe - betaIm * cIm;
            const bim = betaRe * cIm + betaIm * cRe;
            cRe = bre;
            cIm = bim;
          }
          re += cRe;
          im += cIm;
        }
        f64.store(idx, re);
        f64.store(idx + BYTES_PER_ELEMENT, im);
      }
    }
    return;
  }
  else {
    //  Form  C := alpha*A**T*A + beta*C.
    // use Matrix C (nxn as a guide), let it move in sync with (n) columns of A**T

    for (let j: u32 = 0, colBase: u32 = 0, colBaseA_T: u32 = 0; j < n; j++, colBase += NN, colBaseA_T += KK * BYTES_PER_ELEMENT) {
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1)); // complex numbers take 2 positions so "n" is half a column height
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN); // exclusive end position , note again, complex numbers take 2 positions
      for (
        let i: u32 = start, rowBaseA_T: u32 = start % (NN * BYTES_PER_ELEMENT) * k;
        i < stop;
        i += 2 * BYTES_PER_ELEMENT, rowBaseA_T += KK * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {

        let tempRe: f64 = 0;
        let tempIm: f64 = 0;

        // load C 
        const idx: u32 = isPacked ? packCursor : i;
        const cRe: f64 = f64.load(idx);
        const cIm: f64 = f64.load(idx + BYTES_PER_ELEMENT);

        for (let l: u32 = 0; l < KK; l += 2) {

          const aTRowIdx: u32 = ABase + rowBaseA_T + BYTES_PER_ELEMENT * l;
          const aTColidx: u32 = ABase + colBaseA_T + BYTES_PER_ELEMENT * l;

          const aColRe: f64 = f64.load(aTColidx);
          const aColIm: f64 = f64.load(aTColidx + BYTES_PER_ELEMENT);

          const aRowRe: f64 = f64.load(aTRowIdx);
          const aRowIm: f64 = f64.load(aTRowIdx + BYTES_PER_ELEMENT);

          tempRe += aColRe * aRowRe - aColIm * aRowIm;
          tempIm += aColRe * aRowIm + aColIm * aRowRe;
        }

        //multiply with alpha
        let re: f64 = alphaRe * tempRe - alphaIm * tempIm;
        let im: f64 = alphaRe * tempIm + alphaIm * tempRe;

        if (!betaIsZero) {
          re += betaRe * cRe - betaIm * cIm;
          im += betaRe * cIm + betaIm * cRe;
        }
        f64.store(idx, re);
        f64.store(idx + BYTES_PER_ELEMENT, im);
      }
    }
  }
  return;
}

export function csyrk(
  upper: bool,
  transPose: bool,
  n: u32,
  k: u32,
  alphaRe: f32,
  alphaIm: f32,
  betaRe: f32,
  betaIm: f32,
  isPacked: bool
): void {
  const alphaIsZero: bool = alphaRe == 0 && alphaIm == 0
  const betaIsOne: bool = betaRe == 1 && betaIm == 0
  const betaIsZero: bool = betaRe == 0 && betaIm == 0;
  const NN: u32 = n * 2;
  const KK: u32 = k * 2;

  if (n == 0 || ((alphaIsZero || k == 0) && betaIsOne)) return;

  // f32 buffer
  const BYTES_PER_ELEMENT: u32 = 4;

  const ABase: u32 = BYTES_PER_ELEMENT * (isPacked ? n * (n + 1) : n * n * 2);
  // let logIdx: u32 = ABase + BYTES_PER_ELEMENT * (n * k * 2);

  let packCursor: u32 = 0;

  if (alphaIsZero) {
    for (let j: u32 = 0, colBase: u32 = 0; j < n; j++, colBase += NN) { // j = column index
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1));
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN);
      const len: u32 = stop - start;
      if (betaIsZero) {
        if (isPacked) {
          memory.fill(packCursor, 0, len);
          packCursor += len;
        }
        else {
          memory.fill(start, 0, len);
        }
      } else {
        for (let i: u32 = start; i < stop; i += 2 * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {
          const idx: u32 = isPacked ? packCursor : i;
          const cRe: f32 = load<f32>(idx);
          const cIm: f32 = load<f32>(idx + BYTES_PER_ELEMENT);
          const re: f32 = betaRe * cRe - betaIm * cIm;
          const im: f32 = betaRe * cIm + betaIm * cRe;
          store<f32>(idx, re);
          store<f32>(idx + BYTES_PER_ELEMENT, im);
        }
      }
    }
    return;
  }

  if (!transPose) {
    // Form  C := alpha*A*A**T + beta*C.
    // j is the jth column of matrix C
    for (let j: u32 = 0, colBase: u32 = 0; j < n; j++, colBase += NN) {
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1));
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN);
      const len: u32 = stop - start;

      // at this point you want to do  C += alpha*A*A**T
      // A has "k" columns and "n" rows
      for (let i: u32 = start; i < stop; i += 2 * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {
        // column-major
        // because of transpose symmetry we only loop over [1..k]

        let matrixAStart: u32 = i % (NN * BYTES_PER_ELEMENT);
        let matrixATStart: u32 = (j << 1) * BYTES_PER_ELEMENT;

        // the column in C is the column in A transpose with is the row in A

        let tempRe: f32 = 0;
        let tempIm: f32 = 0;

        for (let l: u32 = 0; l < k; l++, matrixAStart += NN * BYTES_PER_ELEMENT, matrixATStart += NN * BYTES_PER_ELEMENT) {
          const aIdx: u32 = ABase + matrixAStart;
          const aTidx: u32 = ABase + matrixATStart;

          const aRe: f32 = f32.load(aIdx);
          const aIm: f32 = f32.load(aIdx + BYTES_PER_ELEMENT);

          const aTRe: f32 = f32.load(aTidx);
          const aTIm: f32 = f32.load(aTidx + BYTES_PER_ELEMENT);

          tempRe += aRe * aTRe - aIm * aTIm;
          tempIm += aRe * aTIm + aIm * aTRe;
        }
        const idx: u32 = isPacked ? packCursor : i;
        let re: f32 = alphaRe * tempRe - alphaIm * tempIm;
        let im: f32 = alphaRe * tempIm + alphaIm * tempRe;
        if (!betaIsZero) {
          let cRe = f32.load(idx);
          let cIm = f32.load(idx + BYTES_PER_ELEMENT)
          if (!betaIsOne) {
            const bre = betaRe * cRe - betaIm * cIm;
            const bim = betaRe * cIm + betaIm * cRe;
            cRe = bre;
            cIm = bim;
          }
          re += cRe;
          im += cIm;
        }
        f32.store(idx, re);
        f32.store(idx + BYTES_PER_ELEMENT, im);
      }
    }
    return;
  }
  else {
    //  Form  C := alpha*A**T*A + beta*C.
    // use Matrix C (nxn as a guide), let it move in sync with (n) columns of A**T

    for (let j: u32 = 0, colBase: u32 = 0, colBaseA_T: u32 = 0; j < n; j++, colBase += NN, colBaseA_T += KK * BYTES_PER_ELEMENT) {
      const start: u32 = BYTES_PER_ELEMENT * (upper ? colBase : colBase + (j << 1)); // complex numbers take 2 positions so "n" is half a column height
      const stop: u32 = BYTES_PER_ELEMENT * (upper ? colBase + ((j + 1) << 1) : colBase + NN); // exclusive end position , note again, complex numbers take 2 positions
      for (
        let i: u32 = start, rowBaseA_T: u32 = start % (NN * BYTES_PER_ELEMENT) * k;
        i < stop;
        i += 2 * BYTES_PER_ELEMENT, rowBaseA_T += KK * BYTES_PER_ELEMENT, packCursor += 2 * BYTES_PER_ELEMENT) {

        let tempRe: f32 = 0;
        let tempIm: f32 = 0;

        // load C 
        const idx: u32 = isPacked ? packCursor : i;
        const cRe: f32 = f32.load(idx);
        const cIm: f32 = f32.load(idx + BYTES_PER_ELEMENT);

        for (let l: u32 = 0; l < KK; l += 2) {

          const aTRowIdx: u32 = ABase + rowBaseA_T + BYTES_PER_ELEMENT * l;
          const aTColidx: u32 = ABase + colBaseA_T + BYTES_PER_ELEMENT * l;

          const aColRe: f32 = f32.load(aTColidx);
          const aColIm: f32 = f32.load(aTColidx + BYTES_PER_ELEMENT);

          const aRowRe: f32 = f32.load(aTRowIdx);
          const aRowIm: f32 = f32.load(aTRowIdx + BYTES_PER_ELEMENT);

          tempRe += aColRe * aRowRe - aColIm * aRowIm;
          tempIm += aColRe * aRowIm + aColIm * aRowRe;
        }

        //multiply with alpha
        let re: f32 = alphaRe * tempRe - alphaIm * tempIm;
        let im: f32 = alphaRe * tempIm + alphaIm * tempRe;

        if (!betaIsZero) {
          re += betaRe * cRe - betaIm * cIm;
          im += betaRe * cIm + betaIm * cRe;
        }
        f32.store(idx, re);
        f32.store(idx + BYTES_PER_ELEMENT, im);
      }
    }
  }
  return;
}
