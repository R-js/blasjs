
import parse from './parse-matrix-from-csv';
import { matrixC, matrixA } from './parse-matrix-from-csv';
import { upper, upperPack, lower, lowerPack } from '../../../../utils/matrix-upper-lower-diag';

describe('level 3 csyrk', function () {
    xit('upper diagonal csyrk', () => {
        const flc = new Float32Array([...parse(matrixC, ',')]);
        const fla = new Float32Array([...parse(matrixA, ',')]);
        console.log(flc);
        console.log(upper(flc, 4, true));
        console.log(upperPack(flc, 4, true));
        console.log(fla);
    });
    it('lower diagonal csyrk', () => {
        const flc = new Float32Array([...parse(matrixC, ',')]);
        console.log(flc);
        console.log(lower(flc, 4, true));
        console.log(lowerPack(flc, 4, true));
    });
});