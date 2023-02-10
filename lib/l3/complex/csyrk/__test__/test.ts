import { resolve } from 'node:path';
import { copyMatricesToWasmMemory, mapWasmMemoryToMatrices } from '@utils/web-assembly-memory';
import {
   upperPack, 
   lowerPack 
} from '@utils/matrix-triangular';
import { loadData } from '@test-helpers/load';
import { csyrk } from '..';


describe('level 3 csyrk C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C', function () {
    describe("quick exit", () => {
        it('n = 0 | alpha = 0+i0 && beta = 1+0i| k = 0 && beta = 1', async () => {
            // n = 0 ,alpha != 0, beta != 0, k != 0
            const ci = new Float32Array(0);
            const ai = new Float32Array(0);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            let [betaRe, betaIm] = [1, 1];
            let [alphaRe, alphaIm] = [1, 1];

            csyrk(true, false, 0, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
           
            expect(ci).toEqualFloatingPointBinary(co);
            expect(ai).toEqualFloatingPointBinary(ao);
            // alpha = 0 && beta = 1
            const ai2 = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci2 = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);

            const result2 = copyMatricesToWasmMemory(false, ci2, ai2);
            const [co2, ao2] = mapWasmMemoryToMatrices(false, result2.storage, ci2.length, ai2.length);
            
            alphaRe = 0;
            alphaIm = 0;
            betaRe = 1;
            betaIm = 0;
            
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
           
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
           
            csyrk(true, false, 4, 0, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
           
            expect(ci3).toEqualFloatingPointBinary(co3);
            expect(ai3).toEqualFloatingPointBinary(ao3);
        });
        it('alpha = 0, upper=true', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const checkC = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;
           
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(checkC).toEqualFloatingPointBinary(co, 18); // matrix C changed
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0, upper=false', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 1;
            const betaIm = 1;
          
            csyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
           
            expect(co).toEqualFloatingPointBinary(cCheck, 18);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('alpha = 0 and beta = 0', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 0;
            const alphaIm = 0;
            const betaRe = 0;
            const betaIm = 0;
            
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(0);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
    });

    describe('fidelity', () => {
        it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;
            
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);

            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);

            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;
            
            csyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);
            
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);
            
            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;
                        
            csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
            
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);
            
            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;
            
            csyrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
            
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);
            
            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0;
            const betaIm = 0;
            
            csyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 15);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        
        it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
            const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
            const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
            const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);
            
            const result = copyMatricesToWasmMemory(false, ci, ai);
            const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ci.length, ai.length);
            
            const alphaRe = 1.2;
            const alphaIm = 0.8;
            const betaRe = 0.6;
            const betaIm = 0.4;
            
            csyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, false, false);
            
            expect(co).toEqualFloatingPointBinary(cCheck, 19);
            expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
        });
        
    });
    describe("quick exit packed", () => {
         it('alpha = 0, upper=true', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const ciPacked = upperPack(ci, 4, true);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper.csv'), ',', true, true, true, false);
             const cCheckPacked =  upperPack(cCheck, 4, true);
             
             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);
             
             const alphaRe = 0;
             const alphaIm = 0;
             const betaRe = 1;
             const betaIm = 1;
             
             csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);
             
             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('alpha = 0, upper=false', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const ciPacked = lowerPack(ci, 4, true);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower.csv'), ',', true, true, true, false);
             const cCheckPacked =  lowerPack(cCheck, 4, true);
             
             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);
             
             const alphaRe = 0;
             const alphaIm = 0;
             const betaRe = 1;
             const betaIm = 1;
             
             csyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);
             
             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('alpha = 0 and beta = 0', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const ciPacked = lowerPack(ci, 4, true);
             
             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);
             
             const alphaRe = 0;
             const alphaIm = 0;
             const betaRe = 0;
             const betaIm = 0;
             
             csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);
             
             expect(co).toEqualFloatingPointBinary(0);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
    });
    describe('fidelity packed', () => {
         it('|alpha| > 1 and beta = 0 (for upper triangular c)', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
             const ciPacked = upperPack(ci, 4, true);
             const cCheckPacked = upperPack(cCheck, 4, true);
             
             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

             const alphaRe = 1.2;
             const alphaIm = 0.8;
             const betaRe = 0;
             const betaIm = 0;

             csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);
             
             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('|alpha| > 1 and beta = 0 (for lower triangular c)', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
             const ciPacked = lowerPack(ci, 4, true);
             const cCheckPacked = lowerPack(cCheck, 4, true);

             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);


             const betaRe = 0;
             const betaIm = 0;
             const alphaRe = 1.2;
             const alphaIm = 0.8;

             csyrk(false, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('|alpha| > 1 and beta = 0.6+0.4i (for upper triangular c)', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha-beta.csv'), ',', true, true, true, false);
             const ciPacked = upperPack(ci, 4, true);
             const cCheckPacked = upperPack(cCheck, 4, true);

             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

             const betaRe = 0.6;
             const betaIm = 0.4;
             const alphaRe = 1.2;
             const alphaIm = 0.8;

             csyrk(true, false, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         
         it('|alpha| > 1 and beta = 0 (for upper triangular c), a is stored row-major', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-upper.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-upper-alpha.csv'), ',', true, true, true, false);
             const ciPacked = upperPack(ci, 4, true);
             const cCheckPacked = upperPack(cCheck, 4, true);

             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

             const betaRe = 0;
             const betaIm = 0;
             const alphaRe = 1.2;
             const alphaIm = 0.8;

             csyrk(true, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);
             
             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('|alpha| > 1 and beta = 0 (for lower triangular c), a is stored row-major', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha.csv'), ',', true, true, true, false);
             const ciPacked = lowerPack(ci, 4, true);
             const cCheckPacked = lowerPack(cCheck, 4, true);

             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

             const betaRe = 0;
             const betaIm = 0;
             const alphaRe = 1.2;
             const alphaIm = 0.8;

             csyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
         it('|alpha| > 1 and beta = 0.6+0.4i (for lower triangular c), a is stored row-major', async () => {
             const ai = await loadData(resolve(__dirname, 'matrix-a-row-major.csv'), ',', true, true, true, false);
             const ci = await loadData(resolve(__dirname, 'matrix-c-lower.csv'), ',', true, true, true, false);
             const cCheck = await loadData(resolve(__dirname, 'matrix-c2-lower-alpha-beta.csv'), ',', true, true, true, false);
             const ciPacked = lowerPack(ci, 4, true);
             const cCheckPacked = lowerPack(cCheck, 4, true);
             
             const result = copyMatricesToWasmMemory(false, ciPacked, ai);
             const [co, ao] = mapWasmMemoryToMatrices(false, result.storage, ciPacked.length, ai.length);

             const betaRe = 0.6;
             const betaIm = 0.4;
             const alphaRe = 1.2;
             const alphaIm = 0.8;
          
             csyrk(false, true, 4, 2, alphaRe, alphaIm, betaRe, betaIm, result.storage, true, false);

             expect(co).toEqualFloatingPointBinary(cCheckPacked, 15);
             expect(ao).toEqualFloatingPointBinary(ai); // matrix A did not change
         });
    });
});