
/*
   Univ. of Tennessee
   Univ. of California Berkeley
   Univ. of Colorado Denver
   NAG Ltd.


   SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
   INCX is INTEGER
   SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
   INCY is INTEGER

     SPARAM is REAL array, dimension (5)
     SPARAM(1)=SFLAG
     SPARAM(2)=SH11
     SPARAM(3)=SH21
     SPARAM(4)=SH12
     SPARAM(5)=SH22
*/


import { FortranArr } from '../../f_func';

export function srotm(
      N: number,
      SX: FortranArr,
      INCX: number,
      SY: FortranArr,
      INCY: number,
      SPARAM: FortranArr,
): void {

      let SH11, SH12, SH21, SH22, W, Z;
      let I, KX, KY, NSTEPS;

      const bx = SX.base;
      const by = SY.base;
      const bp = SPARAM.base;

      let SFLAG = SPARAM.r[1 - bp];

      if (N <= 0 || SFLAG + 2 === 0) return;
      if (INCX === INCY && INCX > 0) {
            NSTEPS = N * INCX;
            if (SFLAG < 0) {
                  SH11 = SPARAM.r[2 - bp];
                  SH12 = SPARAM.r[4 - bp];
                  SH21 = SPARAM.r[3 - bp];
                  SH22 = SPARAM.r[5 - bp];

                  for (I = 1; I <= NSTEPS; I += INCX) {
                        let W = SX.r[I - bx];
                        let Z = SY.r[I - by];

                        SX.r[I - bx] = W * SH11 + Z * SH12;
                        SY.r[I - by] = W * SH21 + Z * SH22;
                  }
            }
            else if (SFLAG === 0) {//checked
                  SH12 = SPARAM.r[4 - bp];
                  SH21 = SPARAM.r[3 - bp];
                  for (I = 1; I <= NSTEPS; I += INCX) {
                        W = SX.r[I - bx];
                        Z = SY.r[I - by];
                        SX.r[I - bx] = W + Z * SH12;
                        SY.r[I - by] = W * SH21 + Z;
                  }
            }
            else {//checked
                  SH11 = SPARAM.r[2 - bp];
                  SH22 = SPARAM.r[5 - bp];
                  for (I = 1; I <= NSTEPS; I += INCX) {
                        W = SX.r[I - bx];
                        Z = SY.r[I - by];
                        SX.r[I - bx] = W * SH11 + Z;
                        SY.r[I - by] = -W + SH22 * Z;
                  }
            }
      }
      //INCX !== INCY || INCX <= 0
      else {
            KX = 1;
            KY = 1;
            if (INCX < 0) KX = 1 + (1 - N) * INCX;
            if (INCY < 0) KY = 1 + (1 - N) * INCY;

            if (SFLAG < 0) {
                  SH11 = SPARAM.r[2 - bp];
                  SH12 = SPARAM.r[4 - bp];
                  SH21 = SPARAM.r[3 - bp];
                  SH22 = SPARAM.r[5 - bp];
                  for (I = 1; I <= N; I++) {
                        W = SX.r[KX - bx];
                        Z = SY.r[KY - by];
                        SX.r[KX - bx] = W * SH11 + Z * SH12;
                        SY.r[KY - by] = W * SH21 + Z * SH22;
                        KX = KX + INCX;
                        KY = KY + INCY;
                  }
            }
            else if (SFLAG === 0) { //checked
                  SH12 = SPARAM.r[4 - bp];
                  SH21 = SPARAM.r[3 - bp];
                  for (I = 1; I <= N; I++) {
                        W = SX.r[KX - bx];
                        Z = SY.r[KY - by];

                        SX.r[KX - bx] = W + Z * SH12;
                        SY.r[KY - by] = W * SH21 + Z;
                        KX += INCX;
                        KY += INCY;
                  }
            }
            else {//checked
                  SH11 = SPARAM.r[2 - bp];
                  SH22 = SPARAM.r[5 - bp];
                  for (I = 1; I <= N; I++) {
                        W = SX.r[KX - bx];
                        Z = SY.r[KY - by];
                        SX.r[KX - bx] = W * SH11 + Z;
                        SY.r[KY - by] = -W + SH22 * Z;
                        KX = KX + INCX;
                        KY = KY + INCY;
                  }
            }
      }
}

