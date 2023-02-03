
import parse from './matrix-csv-parse';
import { matrixC, matrixA } from './matrix-csv-parse';
import { upper, upperPack } from './matrix-upper-lower-diag';

describe('level 3 csyrk', function() {
    it('caxpy', () => {
       const flc = new Float32Array([...parse(matrixC, ',')]);
       const fla = new Float32Array([...parse(matrixA, ',')]);
       console.log(flc);
       console.log(upper(flc, 4, true));
       console.log(upperPack(flc, 4,true));
       console.log(fla);
    });
});