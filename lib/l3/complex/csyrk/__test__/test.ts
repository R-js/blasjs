import { resolve } from 'node:path';
import { upperPack, lowerPack } from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import csyrk from '..';

describe('level 3 csyrk C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C', function () {
    describe("quick exit", () => {
        it('n = 0 | alpha = 0 && beta = 1| k = 0 && beta = 1', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const cloneC = new Float32Array(c);
            const beta = new Float32Array(2);
            const alpha = new Float32Array(2);
            // n = 0 ,alpha != 0, beta != 0, k != 0
            beta.fill(1);
            beta.fill(1);
            csyrk(true, false, 0, 2, alpha, a, beta, c);
            expect(cloneC).toEqualFloatingPointBinary(c);
            // alpha = 0 && beta = 1
            alpha.fill(0);
            beta[0] = 1;
            beta[1] = 0;
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(cloneC).toEqualFloatingPointBinary(c);
            // k = 0 && beta = 1
            alpha.fill(1);
            beta[0] = 1;
            beta[1] = 0;
            csyrk(true, false, 4, 0, alpha, a, beta, c);
            expect(cloneC).toEqualFloatingPointBinary(c);
        });
        it('alpha = 0, upper=true', async () => {
            const a = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const c = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const c2 = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);
            const beta = new Float32Array([1,1]);
            const alpha = new Float32Array(2);
            csyrk(true, false, 4, 2, alpha, a, beta, c);
            expect(c).toEqualFloatingPointBinary(c2, 18);
        });
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
        });
    });
    describe('fidelity', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
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
        });
    });
    describe("quick exit packed", () => {
        it('alpha = 0, upper=true', async () => {
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
        });
    });
    describe('fidelity packed', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
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
        });
    });
});