/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CSYRK  performs one of the symmetric rank k operations
*>
*>    C := alpha*A*A**T + beta*C,
*>
*> or
*>
*>    C := alpha*A**T*A + beta*C,
*>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*> in the second case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*>           triangular  part  of the  array  C  is to be  referenced  as
*>           follows:
*>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*>                                  is to be referenced.
*>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*>                                  is to be referenced.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry,  TRANS  specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
*>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the order of the matrix C.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*>           of  columns   of  the   matrix   A,   and  on   entry   with
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the
*>           matrix A.  K must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension ( LDA, ka ), where ka is
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by n  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*>           be at least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX
*>           On entry, BETA specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension ( LDC, N )
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*>           upper triangular part of the array C must contain the upper
*>           triangular part  of the  symmetric matrix  and the strictly
*>           lower triangular part of C is not referenced.  On exit, the
*>           upper triangular part of the array  C is overwritten by the
*>           upper triangular part of the updated matrix.
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*>           lower triangular part of the array C must contain the lower
*>           triangular part  of the  symmetric matrix  and the strictly
*>           upper triangular part of C is not referenced.  On exit, the
*>           lower triangular part of the array  C is overwritten by the
*>           lower triangular part of the updated matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*/

// multithreaded does not exist in javascript and this avoids creating a new FloatArray 

export default function csyrk<T extends Float32Array | Float64Array>(
    upper: boolean,
    transPose: boolean,
    n: number,
    k: number,
    alpha: T,
    a: T,
    beta: T,
    c: T,
    isPacked = false
): void {
    const alphaIsZero = alpha[0] === 0 && alpha[1] === 0;
    const betaIsOne = beta[0] === 1 && beta[1] === 0;
    const betaIsZero = beta[0] === 0 && beta[1] === 0;
    const NN = n * 2;
    const KK = k * 2;

    //  Quick return if possible.
    if (n === 0 || ((alphaIsZero || k === 0) && betaIsOne)) return;

    //  And when  alpha.eq.zero.
    let packCursor = 0;
    if (alphaIsZero) {
        // column major matrix layout
        for (let j = 0, colBase = 0; j < n; j++, colBase += NN) { // j = column index
            const start = upper ? colBase : colBase + (j << 1); // complex numbers take 2 positions so "n" is half a column height
            const stop = upper ? colBase + ((j + 1) << 1) : colBase + NN; // exclusive end position , note again, complex numbers take 2 positions
            if (betaIsZero) {
                if (isPacked){
                    const len = stop - start;
                    c.fill(0, packCursor, packCursor + len);
                    packCursor += len;
                 }
                 else {
                    c.fill(0, start, stop);
                 }
            } else {
                for (let i = start; i < stop; i += 2, packCursor += 2) {
                    const idx = isPacked ? packCursor : i;
                    const re = beta[0] * c[idx] - beta[1] * c[idx + 1];
                    const im = beta[0] * c[idx] + beta[1] * c[idx + 1];
                    c[idx] = re;
                    c[idx + 1] = im;
                }
            }
        }
        return;
    }

    // Start the operations.
    if (transPose === false) {
        // Form  C := alpha*A*A**T + beta*C.
        // j is the jth column of matrix C
        for (let j = 0, colBaseC = 0; j < n; j++, colBaseC += NN) {
            const packStart = packCursor;
            const start = upper ? colBaseC : colBaseC + (j << 1); // complex numbers take 2 positions so "n" is half a column height
            const stop = upper ? colBaseC + ((j + 1) << 1) : colBaseC + NN; // exclusive end position , note again, complex numbers take 2 positions
            if (betaIsZero) {
                // clear out matrix C to do later C := alpha*A*A**T 
                if (isPacked){
                    const len = stop - start;
                    c.fill(0, packCursor, packCursor + len);
                 }
                 else {
                    c.fill(0, start, stop);
                 }
            }
            // aka C := alpha*A*A**T + beta*C, 
            else if (!betaIsOne) {
                for (let i = start; i < stop; i += 2, packCursor += 2) {
                    const idx = isPacked ? packCursor : i;
                    const re = beta[0] * c[idx] - beta[1] * c[idx + 1];
                    const im = beta[0] * c[idx + 1] + beta[1] * c[idx];
                    c[idx] = re;
                    c[idx + 1] = im;
                }
              
            }
            // at this point you want to do  C += alpha*A*A**T
            // A has "k" columns and "n" rows
            packCursor = packStart; // rewind
            for (let i = start; i < stop; i += 2, packCursor += 2) {
                // because of transpose symmetry we only loop over [1..k]
                let matrixAStart = i % NN; // we dont devide by 2 because A is also complex, both i and NN take are vars in the space of the complex matrix C
                // the column in C is the column in A transpose with is the row in A
                let matrixATStart = j << 1;
                let tempRe = 0;
                let tempIm = 0;
                for (let l = 0; l < k; l++, matrixAStart += NN, matrixATStart += NN) {
                    tempRe += a[matrixAStart] * a[matrixATStart] - a[matrixAStart + 1] * a[matrixATStart + 1];
                    tempIm += a[matrixAStart] * a[matrixATStart + 1] + a[matrixAStart + 1] * a[matrixATStart]
                }
                const idx = isPacked ? packCursor : i;
                c[idx] += alpha[0] * tempRe - alpha[1] * tempIm;
                c[idx + 1] += alpha[0] * tempIm + alpha[1] * tempRe;
            }
        }
    } else {
        //  Form  C := alpha*A**T*A + beta*C.
        // use Matrix C (nxn as a guide), let it move in sync with (n) columns of A**T
        for (let j = 0, colBaseC = 0, colBaseA_T = 0; j < n; j++, colBaseC += NN, colBaseA_T += KK) {
            const start = upper ? colBaseC : colBaseC + (j << 1); // complex numbers take 2 positions so "n" is half a column height
            const stop = upper ? colBaseC + ((j + 1) << 1) : colBaseC + NN; // exclusive end position , note again, complex numbers take 2 positions
            for (let i = start, rowBaseA_T = (start % NN) * k; i < stop; i += 2, rowBaseA_T += KK, packCursor+=2) {
                // row in C = i % NN;
                // this is the column in AT because row-major
                // so (i % NN) * KK, 
                //
                // col in C = j
                // this is also the column in A_T (stored row major)
                let tempRe = 0;
                let tempIm = 0;
                for (let l = 0; l < KK; l += 2) {
                    tempRe += a[colBaseA_T + l] * a[rowBaseA_T + l] - a[colBaseA_T + l + 1] * a[rowBaseA_T + l + 1];
                    tempIm += a[colBaseA_T + l] * a[rowBaseA_T + l + 1] + a[colBaseA_T + l + 1] * a[rowBaseA_T + l];
                }
                //multiply with alpha
                let re = alpha[0] * tempRe - alpha[1] * tempIm;
                let im = alpha[0] * tempIm + alpha[1] * tempRe;
                const idx = isPacked ? packCursor : i;
                if (!betaIsZero) {
                    re += beta[0] * c[idx] - beta[1] * c[idx + 1];
                    im += beta[0] * c[idx + 1] + beta[1] * c[idx];
                }
                c[idx] = re;
                c[idx + 1] = im
            }
        }
    }
}
