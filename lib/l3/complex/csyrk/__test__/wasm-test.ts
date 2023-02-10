import { resolve } from 'node:path';
import { copyMatricesToWasmMemory, mapWasmMemoryToMatrices } from '@utils/web-assembly-memory';
import {
    upperPack,
    lowerPack
} from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import { initWasmCSYRK64 } from '..';
import type { CSYRKfn } from '..';


const fp64 = true;
const packed = true;

describe('level 3 (64fp) zsyrk C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C', function () {
    let storage: WebAssembly.Memory;
    let zsyrk: CSYRKfn;
    beforeAll(() => {
        const { storage: _1, zsyrk: _2 } = initWasmCSYRK64();
        storage = _1;
        zsyrk = _2;
    });
    beforeEach(() => {
        const arr = new Uint8Array(storage.buffer);
        arr.fill(0);
    });
    describe("quick exit", () => {
        it('n = 0 | alpha = 0+i0 && beta = 1+0i| k = 0 && beta = 1', async () => {
            // n = 0 ,alpha != 0, beta != 0, k != 0
            const ci = new Float64Array(0);
            const ai = new Float64Array(0);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            let [betaRe, betaIm] = [1, 1];
            let [alphaRe, alphaIm] = [1, 1];

            zsyrk(true, false, 0, 2, alphaRe, alphaIm, betaRe, betaIm, false);
            expect(ci).toEqualFloatingPointBinary(co);
            expect(ai).toEqualFloatingPointBinary(ao);
            // alpha = 0 && beta = 1
            const ai2 = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci2 = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);

            const result2 = copyMatricesToWasmMemory(fp64, storage, ci2, ai2);
            expect(result2.storage).toBe(storage);
            const [co2, ao2] = mapWasmMemoryToMatrices(fp64, result2.storage, ci2.length, ai2.length);

            alphaRe = 0;
            alphaIm = 0;
            betaRe = 1;
            betaIm = 0;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(ci2).toEqualFloatingPointBinary(co2);
            expect(ai2).toEqualFloatingPointBinary(ao2);
            // k = 0 && beta = 1, alpha != 0
            const ci3 = ci2;
            const ai3 = new Float64Array(0);

            const result3 = copyMatricesToWasmMemory(fp64, storage, ci3, ai3);
            const [co3, ao3] = mapWasmMemoryToMatrices(fp64, result3.storage, ci3.length, ai3.length);

            alphaRe = 0;
            alphaIm = 1;
            betaRe = 1;
            betaIm = 0;

            zsyrk(true, false, 4, 0, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(ci3).toEqualFloatingPointBinary(co3);
            expect(ai3).toEqualFloatingPointBinary(ao3);
        });
        it('alpha = 0, upper=true', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, fp64);
            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(cCheck).toEqualFloatingPointBinary(co, 43); // matrix C changed
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, fp64);
            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            zsyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 42);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0 and beta = 0', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 0.0;
            const betaIm = 0.0;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);
            expect(co).toEqualFloatingPointBinary(0);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
    describe('fidelity', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            //const logging = new Uint32Array(result.storage.buffer, result.byteLength);
            //logging.fill(1);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;
            //logging;
            zsyrk(
                true, // upper
                false, // transpose
                4, // n 
                2, // k
                alphaRe,
                alphaIm,
                betaRe,
                betaIm,
                false // isPacked
            );

            expect(co).toEqualFloatingPointBinary(cCheck, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            zsyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 29);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            //const logging = new Uint32Array(result.storage.buffer, result.byteLength);
            //logging.fill(1);

            zsyrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);
            expect(co).toEqualFloatingPointBinary(cCheck, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            //const logging = new Uint32Array(result.storage.buffer, result.byteLength);
            //logging.fill(1);

            zsyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);


            expect(co).toEqualFloatingPointBinary(cCheck, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, fp64);

            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;

            zsyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 41);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
    describe("quick exit packed", () => {
        it('alpha = 0, upper=true', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const ciPacked = upperPack(ci, 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, fp64);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, true);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 43);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, fp64);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            zsyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, true);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 42);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0 and beta = 0', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const ciPacked = lowerPack(ci, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 0;
            const betaIm = 0;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, true);

            expect(co).toEqualFloatingPointBinary(0);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
    describe('fidelity packed', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, fp64);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.0;
            const betaIm = 0.0;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, fp64);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);


            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            zsyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, fp64);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 29);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, fp64);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            zsyrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);
            
            expect(co).toEqualFloatingPointBinary(cCheckPacked, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, fp64);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            zsyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 37);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, fp64);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, fp64);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);
            
            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);

            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;
         
            zsyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, packed);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 41);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
});