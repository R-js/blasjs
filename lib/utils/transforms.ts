import { isqrt, isqrtUpper } from './isqrt';

/*
  0 4 7 9  
  1 5 8 .
  2 6 . .
  3 . . .
*/

export function transformUpperPackedIdxToColumnRow(idx: number): { col: number, row: number } {
    const col = (isqrt((idx << 3) + 1) - 1) >> 1;
    /**
        * 0 1 2 3 4
        * ----------
        * 0 1 3 6 10
        * . 2 4 7 11
        * . . 5 8 12
        * . . . 9 13
        * . . . . 14
    */
    const row = idx - (col * (col + 1) >> 1);
    return { col, row };
}

/*
  0 . . .  
  1 4 . .
  2 5 7 .
  3 6 8 9
*/

export function transformLowerPackedIdxToColumnRow(idx: number, N: number): { col: number, row: number } {
    const b = ((N << 1) + 1);
    const bb = b * b;
    const ac4 = idx << 3;
    const D = bb - ac4;
    const col = (b - isqrtUpper(D)) >> 1;
    /**
        * 0 1 .2 .3 .4
        * ----------
        * 0 . .. .. ..
        * 1 5 .. .. .. 
        * 2 6 .9 .. ..
        * 3 7 10 12 ..
        * 4 8 11 13 14
    */
    const base = -col * (col - b) >> 1;
    const row = idx - base + col;
    return { col, row };
}

export function getCol(N: number, r: number): number {
    const b = 2 * N + 1;
    const D = b * b - 8 * r;
    console.log({ b, D, srD: Math.sqrt(D) });
    return (b - Math.ceil(Math.sqrt(D))) >> 1;
}

export function back(N: number, r: number): number {
    const b = ((N << 1) + 1);
    const bb = b * b;
    const ac4 = r << 3;
    const D = bb - ac4;
    //console.log(`isqrt:${isqrt(D)}, sqrt:${Math.sqrt(D)}, b-isqrt:${b - isqrt(D)}, b-sqrt:${b - Math.sqrt(D)}`);
    const col1 = (b - Math.sqrt(D)) >> 1;
    return col1;
}
