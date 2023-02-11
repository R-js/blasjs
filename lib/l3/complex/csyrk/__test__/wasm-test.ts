import { resolve } from 'node:path';
import { readFile } from 'node:fs/promises';

import generateMatrix from './fixture-generation';
import { copyMatricesToWasmMemory, mapWasmMemoryToMatrices } from '@utils/web-assembly-memory';
import {
    upperPack,
    lowerPack,
    lower,
    upper,
    transposeStorage
} from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import { initWasmCSYRK64 } from '..';
import type { CSYRKfn } from '..';
import ms from 'ms';


const fp64 = true;
const packed = true;
const globalA = generateMatrix(true, 1234, 4, 2, true);
const globalAT = transposeStorage(globalA, 4, 2, true);
const globalC = generateMatrix(true, 7894, 4, 4, true);


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
            const ai2 = globalA.slice();
            const ci2 = upper(globalC.slice(), 4, true);

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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, fp64);
            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            zsyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, false);

            expect(cCheck).toEqualFloatingPointBinary(co, 41); // matrix C changed
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);

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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 41);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
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
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
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
    describe.skip('huge-matrices test', () => {
        it('packed: C=2048x2048 complex, A=2048x1024', async () => {
            const n = 2048;
            const k = 1024;
            const t0 = performance.now();
            //const ai = generateMatrix(fp64, 123456, n, k, true, false);
            //await writeFile( resolve(__dirname, 'matrix-a-2048x1024.bin'), new Uint8Array(ai.buffer));
            const aiRaw = await readFile( resolve(__dirname, 'matrix-a-2048x1024.bin'));
            const ai = new Float64Array(aiRaw);
            //const ci = generateMatrix(fp64, 74231, n, n, true, false);
            //await writeFile( resolve(__dirname, 'matrix-c-2048x2024.bin'), new Uint8Array(ci.buffer));
            const ciRaw = await readFile( resolve(__dirname, 'matrix-c-2048x2048.bin'));
            const ci = new Float64Array(ciRaw);
            const ciPacked = lowerPack(ci, n, true);

            const result = copyMatricesToWasmMemory(fp64, storage, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ciPacked.length, ai.length);
            
            const t1 = performance.now();
            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            const t2 = performance.now();
            zsyrk(false, false, n, k, alphaRe, alphaIm, betaRe, betaIm, true);

            const t3 = performance.now();
            console.log(`packed: loading+preparing = ${ms((t1 - t0) / 1000)}, calculating=${ms((t3 - t2) / 1000)}`);
            console.log(ai.length, ci.length, ao.length, co.length);
        });
        it('C=2048x2048 complex, A=2048x1024', async () => {
            const t0 = Date.now();
            const n = 2048;
            const k = 1024;
            const ai = await loadData(resolve(__dirname, 'matrix-a-2048x1024.csv'), ',', true, true, true, fp64);
            const ci = await loadData(resolve(__dirname, 'matrix-c-2048x2048-lower.csv'), ',', true, true, true, fp64);
           
            return;
            const result = copyMatricesToWasmMemory(fp64, storage, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(fp64, result.storage, ci.length, ai.length);
            const t1 = Date.now();
            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const t2 = Date.now();
            zsyrk(false, false, n, k, alphaRe, alphaIm, betaRe, betaIm, false);
            const t3 = Date.now();
            console.log(`unpacked: loading+preparing = ${ms(t1 - t0)}, calculating=${ms(t3 - t2)}`);
            console.log(ao.length, co.length);
        });
    });
});