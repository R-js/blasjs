/**
  Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
  Bogers, J.K.F. blasjs 03/2018 (jkfbogers@gmail.com)

  REFERENCES

  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  Krogh, Basic linear algebra subprograms for Fortran
  usage, Algorithm No. 539, Transactions on Mathematical
  Software 5, 3 (September 1979), pp. 308-323.

*/


import { FortranArr } from '../../f_func';

export function sdotdot(
    n: number,
    sb: number,
    sx: FortranArr,
    incx: number,
    sy: FortranArr,
    incy: number): number {

    let sdot = sb;

    const bx = sx.base;
    const by = sy.base;


    if (n <= 0) {
        return sdot;
    }
    if (incx === incy && incx > 0) {
        let ns = n * incx;
        for (let i = 1; i <= ns; i += incx) {
            sdot += sx.r[i - bx] * sy.r[i - by];
        }
    }
    else {
        let kx = 1;
        let ky = 1;
        if (incx < 0) kx = 1 + (1 - n) * incx;
        if (incy < 0) ky = 1 + (1 - n) * incy;
        for (let i = 1; i <= n; i++) {
            sdot += sx.r[kx - bx] * sy.r[ky - by];
            kx += incx;
            ky += incy;
        }
    }
    return sdot;
}
