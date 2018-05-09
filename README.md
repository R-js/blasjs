# BLASjs  <span style="font-size:small" ><span style="color:red; font-weight: bold;">B</span>asic <span style="color:red; font-weight: bold;">L</span>inear <span style="color:red; font-weight: bold;">A</span>lgebra <span style="color:red; font-weight: bold;">S</span>ubprograms.</span>

This is a 100% Pure Javascript ( TypeScript ) re-write of the reference implementation `Basic Linear Algebra SubPrograms` (BLAS) numerical library found [here][blas-site].
This is a manual re-write, ["emscripten"](https://kripken.github.io/emscripten-site) was not used.

#### summary

BLASjs contains all the functions (Complex, Real) of the reference implementation capable for `32 bit` and `64 bit` floating point arithmatic:

* 100% code coverage
* 1003 tests
* Output off all tests equal to the BLAS Fortran reference implementation.
* Level 1: all vector-vector operations implemented. 
* Level 2: all vector-matrix operations implemented.
* Level 3: all matrix-matrix operations implemented.
* Helper functions to easy porting of fortran BLAS usage to Javascript.

#### Node and Web

The library is an UMD library, it can be used in a web client
as in server side node environment.

## Installation

#### node

```bash
npm i blasjs
```

#### web

The module directory contains a minimized bundle for use in html `<script>` tag. The library is attached to the `window.BLAS` object after loading.

```html
<!-- script src="your_server_url/blasjs.min.js"></script -->
<!-- this example uses unpkg as CDN -->
<script src="https://unpkg.com/blasjs@latest/dist/lib/blasjs.min.js">
<script>
  const blas = window.BLAS; //UMD exposes it as BLAS

  //fetch some level3 matrix-matrix operations
  const {
      level3: { zsyrk, ztrmm, ztrsm } //level 3 double precision complex operations
   } = blas;
</script>
```

# Table of Contents

* [Differences with R](#differences-with-r)
* [*Read this first*: Helper functions](#helper-functions-for-porting-blas-programs)
* [Level 1](#level-1)
    * caxpy
    * 


| base functions |
| -------------- |  |  |  |  |  |
| caxpy.f        | level 1 | caxpy.ts |  |  | y = a*x + y |
| ccopy.f        |
| cdotc.f        |
| cdotu.f        |
| cgbmv.f        |
| cgemm.f        |
| cgemv.f        |
| cgerc.f        |
| cgeru.f        |
| chbmv.f        |
| chemm.f        |
| chemv.f        |
| cher.f         |
| cher2.f        |
| cher2k.f       |
| cherk.f        |
| chpmv.f        |
| chpr.f         |
| chpr2.f        |
| crotg.f        |
| cscal.f        |
| csrot.f        |
| csscal.f       |
| cswap.f        |
| csymm.f        |
| csyr2k.f       |
| csyrk.f        |
| ctbmv.f        |
| ctbsv.f        |
| ctpmv.f        |
| ctpsv.f        |
| ctrmm.f        |
| ctrmv.f        |
| ctrsm.f        |
| ctrsv.f        |
| dasum.f        |
| daxpy.f        |
| dcabs1.f       |
| dcopy.f        |
| ddot.f         |
| dgbmv.f        |
| dgemm.f        |
| dgemv.f        |
| dger.f         |
| dnrm2.f        |
| drot.f         |
| drotg.f        |
| drotm.f        |
| drotmg.f       |
| dsbmv.f        |
| dscal.f        |
| dsdot.f        |
| dspmv.f        |
| dspr.f         |
| dspr2.f        |
| dswap.f        |
| dsymm.f        |
| dsymv.f        |
| dsyr.f         |
| dsyr2.f        |
| dsyr2k.f       |
| dsyrk.f        |
| dtbmv.f        |
| dtbsv.f        |
| dtpmv.f        |
| dtpsv.f        |
| dtrmm.f        |
| dtrmv.f        |
| dtrsm.f        |
| dtrsv.f        |
| dzasum.f       |
| dznrm2.f       |
| icamax.f       |
| idamax.f       |
| isamax.f       |
| izamax.f       |
| lsame.f        |
| sasum.f        |
| saxpy.f        |
| scabs1.f       |
| scasum.f       |
| scnrm2.f       |
| scopy.f        |
| sdot.f         |
| sdsdot.f       |
| sgbmv.f        |
| sgemm.f        |
| sgemv.f        |
| sger.f         |
| snrm2.f        |
| srot.f         |
| srotg.f        | level:1 | srotg.ts | 27 feb 2018 | N/A | setup [given's rotation][srotg] |
| srotm.f        |
| srotmg.f       | level: 1 | srotlg.ts | 27 feb 2018 | N/A | setup [modified Givens rotation][givenmodified] |
|                |
| ssbmv.f        |
| sscal.f        |
| sspmv.f        |
| sspr.f         |
| sspr2.f        |
| sswap.f        |
| ssymm.f        |
| ssymv.f        |
| ssyr.f         |
| ssyr2.f        |
| ssyr2k.f       |
| ssyrk.f        |
| stbmv.f        |
| stbsv.f        |
| stpmv.f        |
| stpsv.f        |
| strmm.f        |
| strmv.f        |
| strsm.f        |
| strsv.f        |
| xerbla.f       |
| xerbla_array.f |
| zaxpy.f        |
| zcopy.f        |
| zdotc.f        |
| zdotu.f        |
| zdrot.f        |
| zdscal.f       |
| zgbmv.f        |
| zgemm.f        |
| zgemv.f        |
| zgerc.f        |
| zgeru.f        |
| zhbmv.f        |
| zhemm.f        |
| zhemv.f        |
| zher.f         |
| zher2.f        |
| zher2k.f       |
| zherk.f        |
| zhpmv.f        |
| zhpr.f         |
| zhpr2.f        |
| zrotg.f        |
| zscal.f        |
| zswap.f        |
| zsymm.f        |
| zsyr2k.f       |
| zsyrk.f        |
| ztbmv.f        |
| ztbsv.f        |
| ztpmv.f        |
| ztpsv.f        |
| ztrmm.f        |
| ztrmv.f        |
| ztrsm.f        |
| ztrsv.f        |

smldkjqsefhqdzer

SROTG - setup Givens rotation

SROTMG - setup modified Givens rotation

SROT - apply Givens rotation

SROTM - apply modified Givens rotation

SSWAP - swap x and y

SSCAL - x = a*x

SCOPY - copy x into y

SAXPY - y = a*x + y

SDOT - dot product

SDSDOT - dot product with extended precision accumulation

SNRM2 - Euclidean norm

SCNRM2- Euclidean norm

SASUM - sum of absolute values

ISAMAX - index of max abs value

DOUBLE

DROTG - setup Givens rotation

DROTMG - setup modified Givens rotation

DROT - apply Givens rotation

DROTM - apply modified Givens rotation

DSWAP - swap x and y

DSCAL - x = a*x

DCOPY - copy x into y

DAXPY - y = a*x + y

DDOT - dot product

DSDOT - dot product with extended precision accumulation

DNRM2 - Euclidean norm

DZNRM2 - Euclidean norm

DASUM - sum of absolute values

IDAMAX - index of max abs value

COMPLEX

CROTG - setup Givens rotation

CSROT - apply Givens rotation

CSWAP - swap x and y

CSCAL - x = a*x

CSSCAL - x = a*x

CCOPY - copy x into y

CAXPY - y = a*x + y

CDOTU - dot product

CDOTC - dot product, conjugating the first vector

SCASUM - sum of absolute values

ICAMAX - index of max abs value

DOUBLE COMLPEX

ZROTG - setup Givens rotation

ZDROTF - apply Givens rotation

ZSWAP - swap x and y

ZSCAL - x = a*x

ZDSCAL - x = a*x

ZCOPY - copy x into y

ZAXPY - y = a*x + y

ZDOTU - dot product

ZDOTC - dot product, conjugating the first vector

DZASUM - sum of absolute values

IZAMAX - index of max abs value

LEVEL 2
Single

SGEMV - matrix vector multiply

SGBMV - banded matrix vector multiply

SSYMV - symmetric matrix vector multiply

SSBMV - symmetric banded matrix vector multiply

SSPMV - symmetric packed matrix vector multiply

STRMV - triangular matrix vector multiply

STBMV - triangular banded matrix vector multiply

STPMV - triangular packed matrix vector multiply

STRSV - solving triangular matrix problems

STBSV - solving triangular banded matrix problems

STPSV - solving triangular packed matrix problems

SGER - performs the rank 1 operation A := alpha*x*y' + A

SSYR - performs the symmetric rank 1 operation A := alpha*x*x' + A

SSPR - symmetric packed rank 1 operation A := alpha*x*x' + A

SSYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

SSPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

Double

DGEMV - matrix vector multiply

DGBMV - banded matrix vector multiply

DSYMV - symmetric matrix vector multiply

DSBMV - symmetric banded matrix vector multiply

DSPMV - symmetric packed matrix vector multiply

DTRMV - triangular matrix vector multiply

DTBMV - triangular banded matrix vector multiply

DTPMV - triangular packed matrix vector multiply

DTRSV - solving triangular matrix problems

DTBSV - solving triangular banded matrix problems

DTPSV - solving triangular packed matrix problems

DGER - performs the rank 1 operation A := alpha*x*y' + A

DSYR - performs the symmetric rank 1 operation A := alpha*x*x' + A

DSPR - symmetric packed rank 1 operation A := alpha*x*x' + A

DSYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

DSPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

Complex

CGEMV - matrix vector multiply

CGBMV - banded matrix vector multiply

CHEMV - hermitian matrix vector multiply

CHBMV - hermitian banded matrix vector multiply

CHPMV - hermitian packed matrix vector multiply

CTRMV - triangular matrix vector multiply

CTBMV - triangular banded matrix vector multiply

CTPMV - triangular packed matrix vector multiply

CTRSV - solving triangular matrix problems

CTBSV - solving triangular banded matrix problems

CTPSV - solving triangular packed matrix problems

CGERU - performs the rank 1 operation A := alpha*x*y' + A

CGERC - performs the rank 1 operation A := alpha*x*conjg( y' ) + A

CHER - hermitian rank 1 operation A := alpha*x*conjg(x') + A

CHPR - hermitian packed rank 1 operation A := alpha*x*conjg( x' ) + A

CHER2 - hermitian rank 2 operation

CHPR2 - hermitian packed rank 2 operation

Double Complex

ZGEMV - matrix vector multiply

ZGBMV - banded matrix vector multiply

ZHEMV - hermitian matrix vector multiply

ZHBMV - hermitian banded matrix vector multiply

ZHPMV - hermitian packed matrix vector multiply

ZTRMV - triangular matrix vector multiply

ZTBMV - triangular banded matrix vector multiply

ZTPMV - triangular packed matrix vector multiply

ZTRSV - solving triangular matrix problems

ZTBSV - solving triangular banded matrix problems

ZTPSV - solving triangular packed matrix problems

ZGERU - performs the rank 1 operation A := alpha*x*y' + A

ZGERC - performs the rank 1 operation A := alpha*x*conjg( y' ) + A

ZHER - hermitian rank 1 operation A := alpha*x*conjg(x') + A

ZHPR - hermitian packed rank 1 operation A := alpha*x*conjg( x' ) + A

ZHER2 - hermitian rank 2 operation

ZHPR2 - hermitian packed rank 2 operation

LEVEL 3
Single

SGEMM - matrix matrix multiply

SSYMM - symmetric matrix matrix multiply

SSYRK - symmetric rank-k update to a matrix

SSYR2K - symmetric rank-2k update to a matrix

STRMM - triangular matrix matrix multiply

STRSM - solving triangular matrix with multiple right hand sides

Double

DGEMM - matrix matrix multiply

DSYMM - symmetric matrix matrix multiply

DSYRK - symmetric rank-k update to a matrix

DSYR2K - symmetric rank-2k update to a matrix

DTRMM - triangular matrix matrix multiply

DTRSM - solving triangular matrix with multiple right hand sides

Complex

CGEMM - matrix matrix multiply

CSYMM - symmetric matrix matrix multiply

CHEMM - hermitian matrix matrix multiply

CSYRK - symmetric rank-k update to a matrix

CHERK - hermitian rank-k update to a matrix

CSYR2K - symmetric rank-2k update to a matrix

CHER2K - hermitian rank-2k update to a matrix

CTRMM - triangular matrix matrix multiply

CTRSM - solving triangular matrix with multiple right hand sides

Double Complex

ZGEMM - matrix matrix multiply

ZSYMM - symmetric matrix matrix multiply

ZHEMM - hermitian matrix matrix multiply

ZSYRK - symmetric rank-k update to a matrix

ZHERK - hermitian rank-k update to a matrix

ZSYR2K - symmetric rank-2k update to a matrix

ZHER2K - hermitian rank-2k update to a matrix

ZTRMM - triangular matrix matrix multiply

ZTRSM - solving triangular matrix with multiple right hand sides

[srotg]: https://en.wikipedia.org/wiki/Givens_rotation
[givenmodified]: https://www.ibm.com/support/knowledgecenter/en/SSFHY8_5.5.0/com.ibm.cluster.essl.v5r5.essl100.doc/am5gr_srotm.htm

[caxpy]:  http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html

[blas-site]: http://www.netlib.org/blas/
[blas-source]: https://github.com/Reference-LAPACK/lapack/tree/master/BLAS



Some notes on Matrix symbols

https://fsymbols.com/generators/overline/

A̅ᵗ   Conjugate transpose
B̅    Conjugate

_`A̅ᵗ ∙ B̅`_

Aᵗ 
```
