
import { cadd, cmul, complex } from './_helpers';
import { Complex } from './_types';
import { scabs1 as SCABS1 } from './atoms/scabs1';

export function caxpy(
      N: number,
      CA: Complex,
      CX: Complex[],
      INCX: number,
      CY: Complex[],
      INCY: number) {

      let IX, IY;

      if (N <= 0) return;
      if (SCABS1(CA) === 0) return;
      if (INCX === 1 && INCY === 1) {
            // code for both increments equal to 1
            for (let I = 0; I < N; I++) {
                  CY[I] = cadd(CY[I], cmul(CA, CX[I]));
            }
      }
      else {
            //code for unequal increments or equal increments
            // not equal to 1
            IX = 0;
            IY = 0;
            if (INCX < 0) IX = (-N + 1) * INCX /*+ 1*/;
            if (INCY < 0) IY = (-N + 1) * INCY /*+ 1*/;
            for (let I = 0; I < N; I++) {
                  CY[IY] = cadd(CY[IY], cmul(CA, CX[IX]));
                  IX = IX + INCX;
                  IY = IY + INCY;
            }
      }
}
