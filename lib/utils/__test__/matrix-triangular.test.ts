import { transformUpperPackedIdxToColumnRow, transformLowerPackedIdxToColumnRow } from '../transforms';

describe('coordinate transform testing', function () {
    describe('upper packed', () => {

        it.only('n=91*92/2 => col=91 row = 0w = 0', () => {
            /**0 1 2 3 4  etc
             * 0 1 3 6 10
             * . 2 4 7 11
             * . . 5 8 12
             * . . . 9 13
             * . . . . 14
             */
            const col = 91;
            expect(transformUpperPackedIdxToColumnRow(col * (col + 1) / 2)).toEqual({ col, row: 0 });
            expect(transformUpperPackedIdxToColumnRow(col * (col + 1) / 2 - 1)).toEqual({ col: col - 1, row: col - 1 });
        });

        describe('upper packed', () => {
            it('N=5, for idx= 7, 10 14, 0, 3 and 9', () => {
                /**
                     * 0 1 .2 .3 .4
                     * ----------
                     * 0 . .. .. ..
                     * 1 5 .. .. .. 
                     * 2 6 .9 .. ..
                     * 3 7 10 12 ..
                     * 4 8 11 13 14
                 */
                const N = 5;
                expect(transformLowerPackedIdxToColumnRow(7, N)).toEqual({ col: 1, row: 3 });
                expect(transformLowerPackedIdxToColumnRow(10, N)).toEqual({ col: 2, row: 3 });
                expect(transformLowerPackedIdxToColumnRow(14, N)).toEqual({ col: 4, row: 4 });
                expect(transformLowerPackedIdxToColumnRow(0, N)).toEqual({ col: 0, row: 0 });
            })
        });
    });
});