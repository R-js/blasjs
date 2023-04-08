import { resolve } from 'node:path';
import { copyMatricesToWasmMemory, mapWasmMemoryToMatrices } from '@utils/web-assembly-memory';
import {
    upperPack,
    lowerPack,
    upper,
    lower,
    transposeStorage
} from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import { syrk } from '..';
import generateMatrix from './fixture-generation';

const globalA = generateMatrix(false, 1234, 4, 2, true);
const globalAT = transposeStorage(globalA, 4, 2, true);
const globalC = generateMatrix(false, 7894, 4, 4, true);


describe('level 3 syrk C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C', function () {
    describe("quick exit", () => {
        it('n = 0 | alpha = 0+i0 && beta = 1+0i| k = 0 && beta = 1', () => {
            // n = 0 ,alpha != 0, beta != 0, k != 0
            const ci = new Float32Array(0);
            const ai = new Float32Array(0);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            let [betaRe, betaIm] = [1, 1];
            let [alphaRe, alphaIm] = [1, 1];

            syrk(true, false, 0, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(ci).toEqualFloatingPointBinary(co);
            expect(ai).toEqualFloatingPointBinary(ao);
            // alpha = 0 && beta = 1
            const ai2 = globalA.slice();
            const ci2 = upper(globalC.slice(), 4, true);
            globalA;
            globalC;
            const result2 = copyMatricesToWasmMemory(false, ci2, ai2);
            const [co2, ao2] = mapWasmMemoryToMatrices(false, result2.storage, ci2.length, ai2.length);

            alphaRe = 0;
            alphaIm = 0;
            betaRe = 1;
            betaIm = 0;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(ci2).toEqualFloatingPointBinary(co2);
            expect(ai2).toEqualFloatingPointBinary(ao2);
            // k = 0 && beta = 1, alpha != 0
            const ci3 = ci2;
            const ai3 = new Float32Array(0);

            const result3 = copyMatricesToWasmMemory(false, ci3, ai3);
            const [co3, ao3] = mapWasmMemoryToMatrices(false, result3.storage, ci3.length, ai3.length);

            alphaRe = 0;
            alphaIm = 1;
            betaRe = 1;
            betaIm = 0;

            syrk(true, false, 4, 0, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(ci3).toEqualFloatingPointBinary(co3);
            expect(ai3).toEqualFloatingPointBinary(ao3);
        });
        it('alpha = 0, upper=true', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const checkC = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(checkC).toEqualFloatingPointBinary(co, 18); // matrix C changed
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            syrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 18);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0 and beta = 0', () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 0;
            const betaIm = 0;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(0);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });

    describe('fidelity', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            syrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice(); // await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            syrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            syrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });

        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;

            syrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);

            expect(co).toEqualFloatingPointBinary(cCheck, 19);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });

    });
    describe("quick exit packed", () => {
        it('alpha = 0, upper=true', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const ciPacked = upperPack(ci, 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;

            syrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0 and beta = 0', () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const ciPacked = lowerPack(ci, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 0;
            const betaIm = 0;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(0);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
    describe('fidelity packed', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const ai = globalA.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);


            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            syrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const ai = globalA.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            syrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });

        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice();
            const ci = upper(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
            const ciPacked = upperPack(ci, 4, true);
            const cCheckPacked = upperPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            syrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const betaRe = 0;
            const betaIm = 0;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            syrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const ai = globalAT.slice();
            const ci = lower(globalC.slice(), 4, true); //await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);
            const ciPacked = lowerPack(ci, 4, true);
            const cCheckPacked = lowerPack(cCheck, 4, true);

            const result = copyMatricesToWasmMemory(false, ciPacked, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

            const betaRe = 0.6;
            const betaIm = 0.4;
            const alphaRe = 1.2;
            const alphaIm = 0.8;

            syrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

            expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });
});