# BLASjs  (<span style="font-size:small" ><span style="color:red; font-weight: bold;">B</span>asic <span style="color:red; font-weight: bold;">L</span>inear <span style="color:red; font-weight: bold;">A</span>lgebra <span style="color:red; font-weight: bold;">S</span>ubprograms.</span>)

This is a 100% Pure Javascript ( TypeScript ) re-write of the reference implementation `Basic Linear Algebra SubPrograms` (BLAS) numerical library found [here][blas-site].
This is a manual re-write, ["emscripten"](https://kripken.github.io/emscripten-site) was not used.

#### summary

BLASjs contains all the functions (Complex, Real) of the reference implementation capable for `32 bit` and `64 bit` floating point arithmatic:

* 100% code coverage
* 1003 tests
* Output off all tests equal to the BLAS FORTRAN reference implementation.
* Level 1: all vector-vector operations implemented.
* Level 2: all vector-matrix operations implemented.
* Level 3: all matrix-matrix operations implemented.
* Helper functions to easy porting of FORTRAN BLAS usage to Javascript.

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

* [Language differences with FORTRAN/BLAS](#language-differences-between-FORTRAN/BLAS-and- Javascript/blasjs)
* [*Read this first*: Helper functions](#helper-functions-for-working-with-blasjs)
* [Level 1 Functions](#level-1)
    * [`isamax`/`idamax`/`izamax`/`icamax` find maximum element of a vector]()
    * [`sasum`/`dasum` sum of the absolute vector element values]()
    * [`saxpy/daxpy/caxpy/zaxpy` ]
    * scnrm2
    * scopy/dcopy/ccopy/zcopy
    * sdot/ddot/cdotc/zdotc
    * sdsdot/dsdot  Compute the inner product of two vectors with extended precision.
    * snrm2/drnm2
    * srot/drot/zdrot/csrot
    * srotg/drotg/crotg/zrotg
    * srotm/drotm
    * srotmg/drotmg
    * sscal/dscal/cscal/csscal
    * sswap/dswap/cswap/zswap
    * cdotu/zdotu
    * scasum/dzasum

* [Level 2](#level-2) 
    * sgbmv/dgbmv/cgbmv/zgbmv
    * sgemv/dgemv/cgemv/zgemv
    * sger/dger
    * ssbmv/dspmv
    * sspr/dspr
    * ssymv/dsymv
    * ssyr/dsyr
    * sspr2/dspr2
    * ssymv, dsymv
    * ssyr/dsyr
    * ssyr2, dsyr2
    * stbmv/dtbmv/ctbmv/ztbmv
    * stbsv/dtbsv/ctbsv/ztbsv
    * stpmv/dtpmv/ctpmv/ztpmv
    * stpsv/dtpsv/ctpsv/ztpsv
    * strmv/dtrmv/ctrmv/ztrmv
    * strsv/dtrsv/ctrsv/ztrsv
    * strmv/dtrmv/ctrmv/ztrmv
    * strsv/dtrsv/ctrsv/ztrsv
    * cgerc/zgerc
    * cgeru/zgeru
    * chbmv/zhbmv
    * chemv/zhemv
    * cher/zher
    * cher2/zher2
    * chpmv/zhpmv
    * chpr/zhpr
    * chpr2/zhpr2

# Language differences between FORTRAN/BLAS and Javascript/blasjs.

FORTRAN language can instrinsicly work with non-zero based multidimensional arrays and complex numbers. Below are some examples from FORTRAN that have no Javascript counterpart. The reference implementation of BLAS functions expect inputs of these types.

`FORTRAN complex scalar, complex array, complex Matrix`

```f77
c    double precision Complex number
     COMPLEX*16 alpha
c
c    double precision Complex array with offset 2
     COMPLEX*16 vector(2,10)
c
c    double precision complex MultiDimensional Array (matrix)
c    rows 1 to 5 , columns 1 to 10
     COMPLEX*16 A(1:5,1:10)
```

To work with the concept of non-zero based arrays and complex numbers in JS, 
these FORTRAN constructs have equivalents in the `blasjs` library.

`blasjs complex scalar, complex array, complex Matrix`

```javascript
  const blas = require('blasjs');

  const {
      util:{

        complex, // function to create complex Object from to real numbers, 

        fortranArrComplex32, // Single precision Real/complex arrays, 

        fortranArrComplex64, // Double precision Real/Complex arrays

        fortranMatrixComplex32, // Single precision 2 dimensional Real/Complex arrays

        fortranMatrixComplex64, // Double precision 2 dimensional Real/Complex arrays
      }
  } = blas;
```

These functions are extensively documented in the [helper functions](#helper-functions-for-working-with-blasjs).
It is recommended you read this introductionary part of the documentation first.
before anything else.

### helper functions

#### arrayrify
  
Force any JS scalar or object into an array.

```javascript

```

#### complex
    each,
    fortranArrComplex32,
    fortranArrComplex64,
    fortranMatrixComplex32,
    fortranMatrixComplex64,
    map,
    multiplexer,
    muxCmplx,
    numberPrecision



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
