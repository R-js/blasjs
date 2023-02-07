import { resolve } from 'node:path';
import { copyMatricesToWasmMemory, mapWasmMemoryToMatrices } from '@utils/web-assembly-memory';
//import { upperPack, lowerPack } from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import csyrk from '..';


describe('level 3 csyrk C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C', function () {
    describe("quick exit", () => {
        it('n = 0 | alpha = 0+i0 && beta = 1+0i| k = 0 && beta = 1', async () => {

            let [betaRe, betaIm] = [1, 1];
            let [alphaRe, alphaIm] = [1, 1];
            // n = 0 ,alpha != 0, beta != 0, k != 0
            const ci = new Float32Array(0);
            const ai = new Float32Array(0);
            const result = copyMatricesToWasmMemory(false, ci, ai);
            // console.log(result);
            csyrk(true, false, 0, 2, alphaRe, alphaIm, betaRe, betaIm, result.mem, false, false);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.mem, ci.length, ai.length);
            expect(ci).toEqualFloatingPointBinary(co);
            expect(ai).toEqualFloatingPointBinary(ao);
            // alpha = 0 && beta = 1
            const ai2 = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci2 = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const result2 = copyMatricesToWasmMemory(false, ci2, ai2);
            // console.log(result2);
            alphaRe = 0;
            alphaIm = 0;
            betaRe = 1;
            betaIm = 0;
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.mem, false, false);
            const [co2, ao2] = mapWasmMemoryToMatrices(false, result2.mem, ci2.length, ai2.length);
            expect(ci2).toEqualFloatingPointBinary(co2);
            expect(ai2).toEqualFloatingPointBinary(ao2);
            // k = 0 && beta = 1, alpha != 0
            const ci3 = ci2;
            const ai3 = new Float32Array(0);
            alphaRe = 0;
            alphaIm = 1;
            betaRe = 1;
            betaIm = 0;
            const result3 = copyMatricesToWasmMemory(false, ci3, ai3);
            // console.log(result3);
            csyrk(true, false, 4, 0, alphaRe, alphaIm, betaRe, betaIm, result.mem, false, false);
            const [co3, ao3] = mapWasmMemoryToMatrices(false, result3.mem, ci3.length, ai3.length);
            expect(ci3).toEqualFloatingPointBinary(co3);
            expect(ai3).toEqualFloatingPointBinary(ao3);
        });
        /*it('alpha = 0, upper=true', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);
            const beta = new Float32Array([1,1]);
            const alpha = new Float32Array(2);
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 18);
        });*/
        /*
        it('alpha = 0, upper=false', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
            const beta = new Float32Array([1,1]);
            const alpha = new Float32Array(2);
            csyrk(false, true, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 18);
        });
        it('alpha = 0 and beta = 0', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const beta = new Float32Array(2);
            const alpha = new Float32Array(2);
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(0);
        });*/
    });

    describe('fidelity', () => {
        /*it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
            const beta = new Float32Array(2);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
            const beta = new Float32Array(2);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(false, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);
            const beta = new Float32Array([0.6, 0.4]);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });
        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
            const beta = new Float32Array(2);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(true, true, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
            const beta = new Float32Array(2);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(false, true, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);
            const beta = new Float32Array([0.6, 0.4]);
            const alpha = new Float32Array([1.2, 0.8]);
            csyrk(false, true, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 15);
        });*/
    });
    describe("quick exit packed", () => {
        /* it('alpha = 0, upper=true', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);
             const co = upperPack(c, 4, true);
             const co2 = upperPack(c2, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array(2);
             beta.fill(1);
             csyrk(true, false, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 18);
         });
         it('alpha = 0, upper=false', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
             const co = lowerPack(c, 4, true);
             const co2 = lowerPack(c2, 4, true);
             const beta = new Float32Array([1, 1]);
             const alpha = new Float32Array(2);
             csyrk(false, true, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 18);
         });
         it('alpha = 0 and beta = 0', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const co = lowerPack(c, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array(2);
             csyrk(true, false, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(0);
         });*/
    });
    describe('fidelity packed', () => {
        /* it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
             const co = upperPack(c, 4, true);
             const co2 = upperPack(c2, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(true, false, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });
         it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
             const co = lowerPack(c, 4, true);
             const co2 = lowerPack(c2, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(false, false, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });
         it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);
             const co = upperPack(c, 4, true);
             const co2 = upperPack(c2, 4, true);
             const beta = new Float32Array([0.6, 0.4]);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(true, false, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });
         it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
             const co = upperPack(c, 4, true);
             const co2 = upperPack(c2, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(true, true, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });
         it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
             const co = lowerPack(c, 4, true);
             const co2 = lowerPack(c2, 4, true);
             const beta = new Float32Array(2);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(false, true, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });
         it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
             const a = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const c = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const c2 = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);
             const co = lowerPack(c, 4, true);
             const co2 = lowerPack(c2, 4, true);
             const beta = new Float32Array([0.6, 0.4]);
             const alpha = new Float32Array([1.2, 0.8]);
             csyrk(false, true, 4, 2, alpha, a, beta, co, true);
             expect(co).toEqualFloatingPointBinary(co2, 15);
         });*/
    });
});