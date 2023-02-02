
import parse from './parse-matrix-csv';
import { testData } from './parse-matrix-csv';

describe('level 3 csyrk', function() {
    it('caxpy', () => {
       const fla = new Float32Array([...parse(testData, ',')]);
       console.log(fla);
    });
});