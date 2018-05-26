# BLASjs  (<span style="font-size:small" ><span style="color:red; font-weight: bold;">B</span>asic <span style="color:red; font-weight: bold;">L</span>inear <span style="color:red; font-weight: bold;">A</span>lgebra <span style="color:red; font-weight: bold;">S</span>ubprograms</span>)

This is a 100% Pure Javascript ( TypeScript ) re-write of the reference implementation `Basic Linear Algebra SubPrograms` (BLAS) numerical library found [here][blas-site].
This is a full manual re-write, ["emscripten"](https://kripken.github.io/emscripten-site) was not used.

#### summary

BLASjs contains all the functions (Complex, Real) of the reference implementation capable for `32 bit` and `64 bit` floating point arithmatic:

* :ok_hand: 100% code coverage
* 1005 tests
* Output off all tests equal to the BLAS FORTRAN reference implementation.
* Level 1: all vector-vector operations implemented.
* Level 2: all vector-matrix operations implemented.
* Level 3: all matrix-matrix operations implemented.
* Helper functions to ease the porting of FORTRAN BLAS usage to Javascript.

#### Node and Web

The resulting bundled `blasjs` file is an agnostic UMD library, it can be used in a web client
as-well as in a server side node environment.

## Installation

#### node

```bash
npm i blasjs
```

```javascript
//node
   const blas = require('blasjs');
//or typescript
   import * as blas from 'blasjs';
```

#### web

The module directory contains a minimized bundle for use in html `<script>` tag. The library is attached to the `window.BLAS` object after loading.

```html
<!-- script src="your_server_url/blasjs.min.js"></script -->
<!-- this example uses unpkg as CDN -->
<script src="https://unpkg.com/blasjs@latest/dist/lib/blasjs.min.js">
<script>
  const blas = window.BLAS; //UMD exposes it as BLAS

  //fetch some level3 complex 64 bit precision matrix-matrix operations
  const {
      level3: { zsyrk, ztrmm, ztrsm }
   } = blas;
</script>
```

# Table of Contents

* [Language differences with FORTRAN/BLAS](#language-differences-with-fortranblas)
* [*Read this first*: Helper functions](#helper-functions-for-working-with-blasjs)
    * [Types](#Types)
        * [fpArray](#fpArray)
        * [FortranArr](#fortranarr)
        * [Type Complex](#type-complex)
        * [Matrix](#matrix)
            * [Float[32/64]Array Complex number storage for Matrix](#float3264array-complex-number-storage-for-matrix)
            * [Handling FORTRAN matrices](#handling-fortran-matrices-multidimensional-arrays)
            * [Performance](#performance)
            * [Creating new transformed Matrix instances from existing ones](#creating-new-transformed-matrix-instances-from-existing-ones)
                * [Matrix.prototype.setlower](#matrixprototypesetlower)
                * [Matrix.prototype.setUpper](#matrixprototypesetupper)
                * [Matrix.prototype.upperBand](#matrixprototypeupperband)
                * [Matrix.prototype.lowerBand](#matrixprototypelowerband)
                * [Matrix.prototype.real](#matrixprototypereal)
                * [Matrix.prototype.imaginary](#matrixprototypeimaginary)
            * [Packed Matrices](#packed-matrices)
                * [Matrix.prototype.packedUpper](#matrixprototypepackedupper)
                * [Matrix.prototype.packedLower](#matrixprototypepackedlower)
            * [Convert Matrix object to Js Array](#convert-matrix-object-to-a-js-array)
                * [Matrix.prototype.toArr](#matrixprototypetoarr)
            * [Matrix Examples](#matrix-examples)
    * [General Helpers](#general-helpers)
        * [arrayrify](#arrayrify)
        * [complex](#complex)
        * [each](#each)
        * [map](#map)
        * [muxCmplx](#muxCmplx)
        * [numberPrecision](#numberPrecision)
    * [Vector constructors](#vector-constructors)
        * [fortranArrComplex32](#fortranarrcomplex32)
        * [fortranArrComplex64](#fortranarrcomplex64)
        * [Vector Creation Examples](#vector-creation-examples)
    * [Matrix constructors](#matrix-constructors)
        * [fortranMatrixComplex32](#fortranmatrixcomplex32)
        * [fortranMatrixComplex64](#fortranmatrixcomplex64)
        * [Matrix Creation Examples](#matrix-creation-examples)
* [Level 1 Routines](#level-1-routines)
    * [Euclidean norm: √(xᴴ·x) or √(xᵀ·x)]()
        * [`scnrm2`, `dznrm2`, `snrm2`, `dnrm2`](#naming)
    * [Construct a Givens plane rotation]()
        * [`srotg`, `drotg`, `crotg`,`zrotg`](#srotg-drotg-crotgzrotg)
    * [Construct the **modified** Givens rotation matrix `H`]()
        * [`srotmg`,`drotmg`](#srotmg-drotmg)
    * [Apply the modified Givens Transformation]()
        * [`srotm`, `drotm`](#srotm-drotm)
    * [Applies a plane rotation]()
        * [`srot`,`drot`, `csrot`, `zdrot`](#srot-drot-csrot-zdrot)
    * [Scale a vector by a constant]()
        * [`sscal`, `dscal`, `cscal`, `zscal`, `csscal`, `zdscal`](#sscal-dscal-cscal-zscal-csscal-zdscal)
    * [Takes the sum of the absolute values of the components of vector]()
        * [`sasum`, `dasum`, `scasum`, `dzasum`](#sasum-dasum-scasum-dzasum)
    * [Interchanges 2 vectors]()
        * [`cswap`, `dswap`, `sswap`, `zswap`](#cswap-dswap-sswap-zswap)
    * [Dot product of two complex vectors]()
        * [`cdotu`, `cdotc`, `zdotu`, `zdotc`](#cdotu-cdotc-zdotu-zdotc)
    * [Dot product of two vectors]()
        * [`sdot`, `ddot`, `sdsdot`, `dsdot`](#sdot-ddot-sdsdot-dsdot)
    * [Finds the index of the first element having maximum absolut value]()
        * [`isamax`, `idamax`, `icamax`, `izamax`](#isamax-idamax-icamax-izamax)
    * [Copy a vector x to a vector y]()
        * [`scopy`, `dcopy`, `ccopy`, `zcopy`](#scopy-dcopy-ccopy-zcopy)
    * [Constant times a vector plus a vector]()
        * [`saxpy`, `daxpy`, `caxpy`, `zaxpy`](#saxpy-daxpy-caxpy-zaxpy)
* [Level 2 Routines](#level-2-routines)    
    * [The hermitian rank 2 operation A ⟵ α·x·yᴴ + conjg( α )·y·xᴴ + A](#hermitian-rank-2-a--αxyᴴ--conjg-α-yxᴴ--a)
        * [`cher2`, `zher2`, `chpr2`, `zhpr2`](#cher2-zher2-chpr2-zhpr2)
    * [The symmetric rank 2 operation A ⟵ α·x·yᵀ + α·y·xᵀ + A](#symmetric-rank-2-a--αxyᵀ--αyxᵀ--a)
        * [`sspr2`, `dspr2`, `ssyr2`, `dsyr2`](#sspr2-dspr2-ssyr2-dsyr2)
    * [The rank 1 operation A ⟵ α·x·yᴴ + A or A ⟵ α·x·yᵀ + A]()
        * [`sger`, `dger`, `cgerc`, `zgerc`, `cgeru`, `zgeru`](#sger-dger-cgerc-zgerc-cgeru-zgeru)
    * [The hermitian rank 1 operation A ⟵ α·x·xᴴ + A]()
        * [`cher`, `zher`, `chpr`, `zhpr`]()
    * [The symmetric rank 1 operation A ⟵ α·x·xᵀ + A]()
        * [`sspr`, `dspr`, `ssyr`, `dsyr`]()
    * [The matrix-vector operation, y ⟵ α·A·x + β·y, or y ⟵ α·Aᵀ·x + β·y or y ⟵  α·Aᴴ·x + β·y]()
        * [`cgbmv`, `chbmv`, `dgbmv`, `dsbmv`, `sgbmv`, `ssbmv`, `zgbmv`, `zhbmv`]()
        * [`ctbmv`, `dtbmv`, `stbmv`, `ztbmv`]()
        * [`cgemv`, `chemv`, `dgemv`, `sgemv`, `zgemv`, `zhemv`]()
        * [`chpmv`, `dspmv`, `sspmv`, `zhpmv`]()
        * [`dsymv`, `ssymv`]()
    * [The matrix-vector operation, x ⟵ A·x, or x ⟵ Aᵀ·x, or x ⟵ Aᴴ·x]()
        * [`ctpmv`, `ctrmv`, `dtpmv`, `dtrmv`, `stpmv`, `strmv`, `ztpmv`, `ztrmv`]()
    * [Solves a systems of equations A·x = b, or Aᵀ·x = b, or Aᴴ·x = b]()
        * [`ctbsv`, `ctpsv`, `ctrsv`, `dtbsv`, `dtpsv`, `dtrsv`, `stbsv`, `stpsv`, `strsv`, `ztbsv`, `ztpsv`, `ztrsv`]()
* [Level 3 Routines](#level-3-routines)
    * [Hermitian rank 2k operations C ⟵ α·A·Bᴴ + conjg( α )·B·Aᴴ + β·C or C ⟵ α·Aᴴ·B + conjg( α )·Bᴴ·A + β·C]()
        * [`cher2k`, `zher2k`]()
    * [Symmetric rank 2k operations C ⟵ α·A·Bᵀ + α·B·Aᵀ + β·C, or C ⟵ α·Aᵀ·B +  α·Bᵀ·A + β·C]()
        * [`csyr2k`, `dsyr2k`, `ssyr2k`, `zsyr2k`]()
    * [Hermatian rank k operations C ⟵ α·A·Aᴴ + β·C, or C ⟵ α·Aᴴ·A + β·C]()
        * [`cherk`, `zherk`]()
    * [Symmetric rank k operations C ⟵ α·A·Aᵀ + β·C, or C ⟵ α·Aᵀ·A + β·C]()
        * [`csyrk`, `dsyrk`, `ssyrk`, `zsyrk`]()
    * [Matrix-matrix operations C ⟵ α·_f̅(A)_·_h̅(B)_ + β·C]()
        * [`cgemm`, `dgemm`, `sgemm`, `zgemm`]()
        * [`chemm`, `zhemm`]()
        * [`ctrmm`, `dtrmm`, `strmm`, `ztrmm`]()
        * [`csymm`, `dsymm`, `ssymm`, `zsymm`]()
    * [Solves the matrix equations _f̅( A )_·X = α·B, or X·_h̅( A )_ =  α·B]()
        * [`ctrsm`, `dtrsm`, `strsm`, `ztrsm`]()

# Level 1 Routines

Routines categorized as _"Level 1"_ perform scalar-vector and vector-vector operations.

## Euclidean norm of a vector


### `srotg`, `drotg`, `crotg`,`zrotg`
## `srotmg`, `drotmg`
## `srotm`, `drotm`
## `srot`, `drot`, `csrot`, `zdrot`
## `sscal`, `dscal`, `cscal`, `zscal`, `csscal`, `zdscal`
## `sasum`, `dasum`, `scasum`, `dzasum`
## `cswap`, `dswap`, `sswap`, `zswap`
## `cdotu`, `cdotc`, `zdotu`, `zdotc`
## `sdot`, `ddot`, `sdsdot`, `dsdot`
## `isamax`, `idamax`, `icamax`, `izamax`
## `scopy`, `dcopy`, `ccopy`, `zcopy`
## `saxpy`, `daxpy`, `caxpy`, `zaxpy`

# level 2 routines

## Hermitian rank 2: A ⟵ α·x·yᴴ + conjg( α )·y·xᴴ + A

### `cher2`, `zher2`, `chpr2`, `zhpr2`

## Symmetric rank 2: A ⟵ α·x·yᵀ + α·y·xᵀ + A

### `sspr2`, `dspr2`, `ssyr2`, `dsyr2`

## The rank 1: A ⟵ α·x·yᴴ + A or A ⟵ α·x·yᵀ + A

### `sger`, `dger`, `cgerc`, `zgerc`, `cgeru`, `zgeru`

## Hermitian rank 1: A ⟵ α·x·xᴴ + A

### `cher`, `zher`, `chpr`, `zpr`

## Symmetric rank 1: A ⟵ α·x·xᵀ + A

### `sspr`, `dspr`, `ssyr`, `dsyr`

## y ⟵ α·A·x + β·y ∨ y ⟵ α·Aᵀ·x + β·y ∨ y ⟵ α·Aᴴ·x + β·y

### `cgbmv`, `chbmv`, `dgbmv`, `dsbmv`, `sgbmv`, `ssbmv`, `zgbmv`, `zhbmv`

### `ctbmv`, `dtbmv`, `stbmv`, `ztbmv`
### `cgemv`, `chemv`, `dgemv`, `sgemv`, `zgemv`, `zhemv`
### `chpmv`, `dspmv`, `sspmv`, `zhpmv`
### `dsymv`, `ssymv`

## x ⟵ A·x ∨ x ⟵ Aᵀ·x ∨ x ⟵ Aᴴ·x

### `ctpmv`, `ctrmv`, `dtpmv`, `dtrmv`, `stpmv`, `strmv`, `ztpmv`, `ztrmv`

## solves A·x = b ∨ Aᵀ·x = b ∨ Aᴴ·x = b

### `ctbsv`, `ctpsv`, `ctrsv`, `dtbsv`, `dtpsv`, `dtrsv`, `stbsv`, `stpsv`, `strsv`, `ztbsv`, `ztpsv`, `ztrsv`

# level 3 routines

## Hermitian rank 2k: C ⟵ α·A·Bᴴ + conjg( α )·B·Aᴴ + β·C ∨ C ⟵ α·Aᴴ·B + conjg( α )·Bᴴ·A + β·C

### `cher2k`, `zher2k`

## Symmetric rank 2k: C ⟵ α·A·Bᵀ + α·B·Aᵀ + β·C ∨ C ⟵ α·Aᵀ·B +  α·Bᵀ·A + β·C

### `csyr2k`, `dsyr2k`, `ssyr2k`, `zsyr2k`

## Hermatian rank k: C ⟵ α·A·Aᴴ + β·C ∨ C ⟵ α·Aᴴ·A + β·C

### `csyrk`, `dsyrk`, `ssyrk`, `zsyrk`

## Matrix-matrix: C ⟵ α·_f̅(A)_·_h̅(B)_ + β·C

### `cgemm`, `dgemm`, `sgemm`, `zgemm`
### `chemm`, `zhemm`
### `ctrmm`, `dtrmm`, `strmm`, `ztrmm`
### `csymm`, `dsymm`, `ssymm`, `zsymm`

## Solves _f̅( A )_·X = α·B ∨ X·_h̅( A )_ =  α·B

### `ctrsm`, `dtrsm`, `strsm`, `ztrsm`



# Language differences with FORTRAN/BLAS.

FORTRAN language can instrinsicly work with non-zero based multidimensional arrays and complex numbers. Below are some examples from FORTRAN that have no Javascript counterpart. The reference implementation of BLAS functions expect inputs of these types.

_The FORTRAN complex scalar, complex array and complex "Matrix"_

```fortran
!    double precision Complex number
     COMPLEX*16 alpha
!
!    double precision Complex array with offset 2
     COMPLEX*16 vector(2,10)
!
!    double precision complex MultiDimensional Array (matrix)
!    rows 1 to 5 , columns 1 to 10
     COMPLEX*16 A(1:5,1:10)
```

To work with the concept of non-zero based arrays and complex numbers in JS,
these FORTRAN constructs have equivalents in the `blasjs` library.

_The `blasjs` helpers to create complex scalar, complex array and complex "Matrix" objects_

```javascript
  const blas = require('blasjs');

  const {
      helper:{
        /* create complex Object from 2 real numbers */
        complex,

        /* create single precision Real/complex arrays, */
        fortranArrComplex32,

        /* create double precision Real/Complex arrays */
        fortranArrComplex64,

        /* create single precision 2 dimensional Real/Complex arrays */
        fortranMatrixComplex32,

        /* Double precision 2 dimensional Real/Complex arrays */
        fortranMatrixComplex64,
      }
  } = blas;
```

These functions are extensively documented in the [helper functions](#helper-functions-for-working-with-blasjs).
It is recommended you read this introductory part of the documentation first.
before anything else.

# Helper functions

`blasjs`
uses "FORTRAN like" complex number 32/64 bit precision multidimensional complex/real data.
These helper functions have been designed to significantly ease the use of working with these
data types in JavaScript.

## Types

Typescript types/interfaces to mimic FORTRAN native (complex) multidimensional arrays.

### `fpArray`

Wraps JS types [Float32Array][float32-array] and [Float64Array][float64-array] into a single type.

_decl_:

```typescript
export type fpArray = Float32Array | Float64Array;
```

### `FortranArr`

Abstraction of a 1 dimensional single/double precision complex/real FORTRAN array.
Used by [level 1](#level-1) and [level 2](#level-2) `blasjs` functions.
`FortranArr`objects should be created by the [`fortranArrComplex32`][float32-array] and [`fortranArrComplex64`][float64-array] helper functions.

_decl_:

```typescript
export declare type FortranArr = {
    base: number;
    r: fpArray;
    i?: fpArray;
    s: (index: number) => (re?: number, im?: number) => number | Complex;
    toArr: () => Complex[] | number[];
};
```

fields:

* `base`: fortran by default has a 1-value based array. Mimiced by this property.
* `r`: See decl [fpArray](#fpArray). The Real part of complex array.
* `i`: (optional). See decl [fpArray](#fpArray). The Imaginary part of the complex array.
* `s`: set, get values of the array. Uses FORTRAN style array indexes taking the value of `base` into account.
* `toArr` generates an JavaScript array from the `r` and `i` (optional) data.

Usage:

```javascript
const blas = require('blasjs');

const { fortranArrComplex64 } = helper;

// You can also use the helper "complex" or "muxComplex"
// to generate JS complex arrays
const complexDataArr = [
    { re: 1.8, im: -0.2 },
    { re: 2.3, im: 0.6 }
];

// Create an object that mimics FORTRAN COMPLEX*16 SP(2:3)
//    and fill it with above data
const sp = fortranArrComplex64(complexArr)(2);

// fast! normal JS TypedArray access
let re = sp.r[ 2 - sp.base ];
// 1.8

let  im = sp.i[ 2 - sp.base ];
// -0.2

// not so fast, but easier syntax
let v = sp.s(2)(); // Terse syntax,
// { re: 1.8, im: -0.2 }

// sets the value at index 3 to complex: 0.11 - i0.9
//      and returns the old value: 2.3 + i0.6
let old = sp.s(3)(0.11, -0.9);

sp.toArr();
// [ { re:1.8, im: -0.2 },
//   { re:0.11, im: -0.9 } ]
```

_Usage TypeScript:_

```typescript
import {
    // pure types
    Complex,
    fpArray,
    FortranArr,
    // helper
    helper
} from 'blasjs';

const { fortranArrComplex64 } = helper;

const complexArr: Complex[] [
    { re: 1.8, im: -0.2 },
    { re: 2.3, im: 0.6 }
];

// Create an object that mimics FORTRAN COMPLEX*16 SP(2:3)
//    and fill it with above data
const sp: FortranArr = fortranArrComplex64(complexArr)(2);

let re = sp.r[ 2 - sp.base ]; //fastest! direct TypedArray access
// 1.8

let  im = sp.i[ 2 - sp.base ]; //fastest! direct TypedArray access
// -0.2

// not so fast, but easier syntax
let v = sp.s(2)(); // Terse syntax,
// { re: 1.8, im: -0.2 }

// sets the value at index 3 to complex: 0.11 - i0.9
//      and returns the old value: 2.3 + i0.6
let old = sp.s(3)(0.11, -0.9);
// {re: 2.3, im: 0.6 }
```

### `Type Complex`

Typescript definition of a complex scalar.

_decl_:

```typescript
declare type Complex = {
   re: number;
   im?: number;
}
```

Usage:

```typescript
import { Complex /* pure type */ } from 'blasjs';

const complexArr: Complex[] [
    { re: 1.8, im: -0.2 },
    { re: 2.3, im: 0.6 }
];
```

### `Matrix`

The `Matrix` object is the input of many level-2 and level-3 `blasjs` functions.
`Matrix` is created by the helpers [fortranMatrixComplex32](#fortranMatrixComplex32) and
[fortranMatrixComplex64](#fortranMatrixComplex64).
`Matrix` encapsulates objects of [Float32Array][float32-array] or [Float64Array][float64-array], the blasjs.

In this section the internals of `Matrix` are explained in detail and how `blasjs` accesses the data in the JS TypesArrays.

#### Float[32/64]Array Complex number storage for Matrix.

The `Matrix` object has 2 properties `r` and `i` for respectively real and imaginary parts of matrix elements. These are the actual aforementioned JS TypedArrays. The imaginary property part is optional if it is not defined the Matrix represents solely an array of real elements.

```typescript
declare type Matrix = { //Incomplete declaration
    .
    r: Float64Array|Float32Array;
    i: Float64Array|Float32Array;
    .
}
```

#### Handling FORTRAN matrices (multidimensional Arrays).

Contrary to languages like JavaScript. FORTRAN defines arrays ( aka `DIMENSIONS` in FORTRAN lingo ) as 1 based arrays by default.. This can be changed by specifying a different base in the declaration.

Some examples:

```fortran
       DOUBLE PRECISION A1(4)  ! array indexes 1,2,3,4
       DOUBLE PRECISION A2(-1:3)  ! array indexes -1,0,2,3
       DOUBLE PRECISION A3(0:3) ! Javascript like Array with 4 elements
```

This expands to 2-dimensional arrays (matrices).

```fortran
! (default) first index loops from 1 to 4(inclusive), second index loops from 1 to 5(inclusive)
       DOUBLE PRECISION A1(4,5)
! first index loops from -2 to 4(inclusive), second index loops from -5 to -7(inclusive)
       DOUBLE PRECISION A2(-2:4,-5:-7)
```

The values of the FORTRAN array basis are preserved as `rowBase` (first index) and `colBase` (second index).

```typescript
declare type Matrix = { //SHOW PARTIAL TYPE
    .
    rowBase: number;
    colBase: number;
    .
}
```

JavaScript doesn't have the notion of `typed 2-dimensional arrays`. The `Matrix` objects handles this by mapping 2 dimensional arrays to single 1-dimensional array, by serializing data on a column-first basis.

For example the elements 2x2 Matrix will be mapped in a TypedArray as:

```bash
matrix A =
 *           *
 | a11  a12  |
 | a21  a22  |
 *           *

# Stored in TypedArray as
At = [a11,a21, a12, a22]
```

In case of complex values for A, the real part will be stored in `r` and the imaginary part in `i` each in the same column-first manner.

#### performance

Direct access to TypedArrays within the `Matrix` object is the preferable way to get/set matrix data.
Since BLAS (and therefore `blasjs`) functions access matrices mostly to iterate over matrix row's first . It was decided to story 2 dimensional an a column-first basis.

To help with the calculation of finding/setting an element A(i,j) in `Matrix` the following helper member functions have been added to `Matrix`.

```typescript
declare type Matrix = { //SHOW PARTIAL TYPE
    .
    rowBase: number;
    colBase: number;
    nrCols: number;
    nrRows: number;
    .
    colOfEx(number): number;
    coord(col): (row) => number;
    setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
    .
}
```

Explanation:

* `nrRows`: The number of rows in the matrix.
* `nrCols`: The number of columns in the matrix.
* `colofEx`: Calculates the physical location of a `column offset` within the `TypedArray`. Taking int account the column base `colBase` and row base `colBase`. The index of  A(i,j) `=  (j - colBase)*nrRows + i - rowBase`.
* `coord`: Curried, emulates non-zero based FORTRAN index values for 2 dimensional Arrays. The index that is iterated over the least (usually) is used as the first to create the curried function.
* `setCol`: Uses underlying `TypedArray`, `fill` method to set multiple column elements to a single value.

_[See Example](#matrix-examples)_

#### Creating new transformed Matrix instances from existing ones

One can create/transform new Matrix instances form existing onces. A copy of all relevant data is made into the new `Matrix` instance.

#### `Matrix.prototype.slice`

Slices a rectangular piece of data out of an matrix into a new `Matrix` instance. **All arguments are FORTRAN-style non-zero based indexes**.

```typescript
declare type Matrix = { // only "slice" is shown
   .
   slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
   .
}
```

* `rowStart`: The row in the matrix to begin slicing.
* `rowEnd`: The last row to include in the slice.
* `colStart`: The column in the matrix to begin slicing.
* `colEnd`: The last column to include in the slice.

_[See Example](#matrix-examples)_

#### `Matrix.prototype.setLower`

Returns a new Matrix where everything below the matrix diagonal is set to a `value`.
Sets the real (and imaginary part, if it exist) to said value.

```typescript
declare type Matrix = { // only "setLower" is shown.
  .
  setLower(value = 0): Matrix;
  .
}
```

_[See Example](#matrix-examples)_

#### `Matrix.prototype.setUpper`

Returns a new Matrix where everything _below_ the matrix diagonal is set to a `value`.
Sets the real (and imaginary part, if it exist) to said value.

```typescript
declare type Matrix = { //only "setUpper" is shown
  .
  setUpper(value = 0): Matrix;
  .
}
```

_[See Example](#matrix-examples)_

#### `Matrix.prototype.upperBand`

Returns a new `Matrix` object where the `k` super-diagonals are retained into the new copy.
The efficient storage format of `BLAS` band matrices is used.

```typescript
declare type Matrix = { //only "upperBand" is shown
  .
  upperBand(k = nrRows - 1): Matrix;
  .
}
```

The default value for `k` is the the maximum size possible for the number of super-diagonals: ( `nrRows-1` )

_[See Example](#matrix-examples)_

#### `Matrix.prototype.lowerBand`

Returns a new `Matrix` object where the `k` sub-diagonals are retained into the new copy.
The efficient storage format of `BLAS` band matrices is used.

```typescript
declare type Matrix = { // Only "lowerBand" is shown
  .
  lowerBand(k = nrRows-1): Matrix;
  .
}
```

The default value for `k` is the the maximum size possible for the number of sub-diagonals: ( `nrRows-1` )

_[See Example](#matrix-examples)_

#### `Matrix.prototype.real`

Returns a new `Matrix` object where with only real elements (omits the imaginary part during copy).

```typescript
declare type Matrix = { // Only "real" is shown
  .
  real(): Matrix;
  .
}
```

_[See Example](#matrix-examples)_

#### `Matrix.prototype.imaginary`

Returns a new `Matrix` object where with only imaginary part of the element (omits the real part during copy).
**If there were now imaginary elements **

```typescript
declare type Matrix = { // Only "imaginary" is shown.
  .
  imaginary(): Matrix;
  .
}
```

_[See Example](#matrix-examples)_

#### Packed Matrices

BLAS ( and therefore `blasjs` ) can work with upper/lower-matrices and band-matrices in the most compacted form, aka `packed matrices`.
With `packed matrices` there are no unused elements in the matrix (no zeros). Packed matrices are instances of [FortranArr](#fortranarr). BLAS reference implementation in FORTRAN uses 1 dimensional arrays as an analog.

#### `Matrix.prototype.packedUpper`

Creates a packed array from a normal/upper Matrix only referencing the diagonal and super-diagonals.

```typescript
declare type Matrix = { // Only "packedUpper" is shown.
  .
  packedUpper(k = nrRows-1): FortranArr;
  .
}
```

_[See Example](#matrix-examples)_

The default value for `k` is the the maximum size possible for the number of super-diagonals: ( `nrRows-1` )

#### `Matrix.prototype.packedLower`

Creates a packed array from a normal/upper Matrix only referencing the diagonal and sub-diagonals.

```typescript
declare type Matrix = { // Only "packedUpper" is shown.
  .
  packedLower(k = nrRows-1): FortranArr;
  .
}
```

_[See Example](#matrix-examples)_

The default value for `k` is the the maximum size possible for the number of sub-diagonals: ( `nrRows-1` )

```typescript
declare type Matrix = { // Only "packedUpper" is shown.
  .
  packedLower(k = nrRows - 1): FortranArr;
  .
}
```

The default value for `k` is the the maximum size possible for the number of sub-diagonals: ( `nrRows - 1` )

#### Convert Matrix object to a JS array

The `Matrix` object can convert the underlying TypedArray(s) to real JavaScript arrays.

#### `Matrix.prototype.toArr`

Creates a normal JS Array with element of type 'number' or of type [Complex](#type-complex)

```typescript
declare type Matrix = { // Only "toArr" is shown.
  .
  toArr(): number[]|Complex[];
  .
}
```

_[See Example](#matrix-examples)_

#### Summary: Full type declaration of Matrix

Putting it all together, here is the full type declaration of `Matrix`:

```typescript
declare type Matrix = {
     rowBase: number;
     colBase: number;
     nrCols: number;
     nrRows: number;
     r: fpArray;
     i?: fpArray; //optional
     //
     // methods
     //
     colOfEx(column: number): void;
     coord(col: number): (row: number): void;
     setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
     //
     slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
     setLower(value?: number): Matrix;
     setUpper(value?: number): Matrix;
     upperBand(k: number): Matrix;
     lowerBand(k: number): Matrix;
     real(): Matrix;
     imaginary(): Matrix;
     //
     packedUpper(value?: number): FortranArr;
     packedLower(value?: number): FortranArr;
     //
     toArr(): Complex[] | number[];
}
```

#### Matrix Examples

Common usage of the Matrix type.

```javascript
const blas = require('../blasjs');
const { fortranMatrixComplex64 } = blas.helper;

// some matrix data 3x3 array  aka a_row_column

const a11 = { re: .2, im: -.11 };
const a21 = { re: .1, im: -.2 };
const a31 = { re: .3, im: .9 };
const a12 = { re: .4, im: .5 };
const a22 = { re: .9, im: -.34 };
const a32 = { re: -.2, im: .45 };
const a13 = { re: -.1, im: .89 };
const a23 = { re: .43, im: .23 };
const a33 = { re: .23, im: .56 };

//create Matrix A
const A = fortranMatrixComplex64([
    a11, a21, a31, a12, a22, a32, a13, a23, a33
])(3, 3);

// get the second column
const columnj = A.colOfEx(3); // formula: (j - colBase )* nrRows

A.r[A.coord(1, 2)] === a12.re // true

A.slice(1, 2,  2, 3);// creates new matrix with elements from A
/*[
    a12 a13
    a22 a23
]*/

A.setLower(0); // creates new Matrix object from A
/*[
    a11 a12 a13
    0   a22 a23
    0   0   a33
]*/

A.setUpper(0); //creates new Matrix object from A
/*[
    a11 0   0
    a21 a22 0
    a31 a32 a33
]*/

A.upperBand(1); // banded array storage for BLAS(js)
/*[
    0   a12   a23
    a11 a22   a33
]*/

A.lowerBand(1); // banded array storage for BLAS(js)
/*[
   a11  a22  a33
   a21  a32  0    
]*/

const Areal = A.real();
// Areal.i is undefined
// Areal.r =
/*[
    0.2 0.4  -0.1
    0.1 0.9   0.43
    0.3 -0.2, 0.23
]*/

const Aimag = A.imaginary();
// imaginary parts are copied to real side in new Matrix
// Aimag.i is undefined
// Aimag.r =
/*[
    -0.11   0.5,  0.89
     -0.2  -0.34  0.23
      0.9   0.45  0.56
]*/

A.packedUpper(1)
/* [ a11 a12 a22 a23 a 33] */

A.packedLower(1)
/* [ a11 a21 a22 a32 a33] */

A.toArr(); // returns JavaScript Array
/*[
  { re: 0.2, im: -0.11 },
  { re: 0.1, im: -0.2 },
  { re: 0.3, im: 0.9 },
  { re: 0.4, im: 0.5 },
  { re: 0.9, im: -0.34 },
  { re: -0.2, im: 0.45 },
  { re: -0.1, im: 0.89 },
  { re: 0.43, im: 0.23 },
  { re: 0.23, im: 0.56 }
]
*/
```

## General Helpers

Collection of helper function to manipulate common JS array and object types in a functional way.

### `arrayrify`

Creates a new function from an existing one, to add the ability to accept vectorized input.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { arrayrify } } = blas;
const PI = Math.PI;
//
const sin = arrayrify(Math.sin)

sin([PI/3, PI/4, PI/6]); // returns array aswell
// [ 0.866025, 0.7071067811, 0.5 ]

sin(PI/3); // returns scalar
sin( [ PI/3 ] ); // returns scalar
// 0.866025

sin([]) // edge case
// undefined

sin() //
//NaN  same as Math.sin()
```

### `complex`

Mimics the GNU Fortran extension [complex](https://gcc.gnu.org/onlinedocs/gfortran/COMPLEX.html).
Creates a JS object that represents a complex scalar number.
Used by `blasjs` for scalar input arguments.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { complex } } = blas;

const c1 = complex(0.1,0.3);
//c1 = { re: 0.1, im: 0.3 }

const c2 = complex();
//c2 = { re: 0, im: 0 }

const c3 = complex(0.5);
//c3 = { re: 0.5, im:0 }
```

### `each`

Curried functional analog to `Array.prototype.forEach`, but takes arbitrary input.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { each } } = blas;

//Iterates over an object like a map
const curry1 = each( {  hello: 'world', ts: new Date() })
curry1( (val, key) => console.log(`${val} ':'  ${key}`)))
//world: hello
//2018-05-10T13:57:08.923Z : ts

//Handles array also
each( ['a','b','c','d'])( (v,idx) =>console.log(v,idx, typeof idx))
//a 0 number
//b 1 number
//c 2 number
//d 3 number

//Edge cases
each()(console.log)
//nothing happens

each(null)(console.log)
//nothing happens

each([])(console.log)
//nothing happens
```

### `map`

Curried functional analog to `Array.prototype.map`, but takes arbitrary input.

:warning: Forces the output to be a an array regardless of the input.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { map } } = blas;

//trivial
map([1,2,3])(v=>v*2);
//[ 2, 4, 6 ]

//key properties
map({ a:'A', b:'B' })( (val, key) => key+'='+val);
//[ 'a=A', 'b=B' ]

map(null)( v => '/'+v);
//[]

map()( v => '/'+v);
//[]

map()()
//[]
```

### `muxCmplx`

Creates an array of complex numbers from arrayed input.
The result is always an array type.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { muxCmplx } } = blas;

const reals = [ 0.1, -0.2, 0.3, 0.45 ];
const imaginary = [ 0.1, -0.2, 0.3, 0.45 ];

// normal usage
muxCmplx(reals, imaginary)
/*[ { re: 0.1, im: 0.1 },
    { re: -0.2, im: -0.2 },
    { re: 0.3, im: 0.3 },
    { re: 0.45, im: 0.45 } ]*/

//R recycling rule is used
muxCmplx([1,2], imaginary)
/*^[ { re: 1, im: 0.1 },
     { re: 2, im: -0.2 },
     { re: 1, im: 0.3 },
     { re: 2, im: 0.45 } ]*/

//dont care about imaginary
muxCmplx(reals)
/*[ { re: 0.1, im: undefined },
    { re: -0.2, im: undefined },
    { re: 0.3, im: undefined },
    { re: 0.45, im: undefined } ]*/

muxCmplx() //
// [ { re: undefined, im: undefined } ]

muxCmplx(1) //
// [ { re: 1, im: undefined } ]

//3 specify real and imaginary
muxCmplx(1,-2)//
//[ { re: 1, im: -2 } ]
```

### `numberPrecision`

Enforces significant figure of a number, or on the properties of a JS object (deep search) with numeric values.

_Example_:

```javascript
const blas = require('blasjs');

const { helper: { numberPrecision } } = blas;

const _4 = numberPrecision(4);

_4(0.123456789);
//0.1235

_4(123456789)
//123500000

//enforce significance over properties
_4( { car: 'Mazda' , aux: { priceUSD: 24.3253E+3, maxWarpSpeed:3.42111E-4 } } );
//{ car: 'Mazda', aux: { priceUSD: 24330, maxWarpSpeed: 0.0003421 } }

_4([0.123456, 0.78901234]);
//[ 0.1235, 0.789 ]
```

## Vector Constructors

These constructors create the `FortranArr` object for working with single/double precision complex/real Arrays.

### `fortranArrComplex32`

Constructs a [FortranArr](#fortranArr) object using [Float32Array][float32-array] as the underlying array(s) (plural in the case of complex) elements.

```typescript
declare function fortranArrComplex32(
    ...rest: (number | number[] | Complex | Complex[])[]
    ): (offset = 1) => FortranArr;
```

`Argument list`:

* `rest`: takes as input.
    * A single numeric value.
    * A single [`Complex`](#type-complex) object.
    * An array of [`Complex`](#type-complex) objects.
    * An array of number values.
*  `offset`: the Fortran dimension offset (defaults to 1)

See _[Examples](#vector-creation-examples)_

### `fortranArrComplex64`

Constructs a [FortranArr](#fortranArr) object using [Float64Array][float64-array] as the underlying array(s) (plural in the case of complex) elements.

```typescript
declare function fortranArrComplex64(
    ...rest: (number | number[] | Complex | Complex[])[]
    ): (offset = 1) => FortranArr;
```

`Argument list`:

* `rest`: takes as input.
    * A single numeric value.
    * A single [`Complex`](#type-complex) object.
    * An array of [`Complex`](#type-complex) objects.
    * An array of number values.
*  `offset`: the Fortran dimension offset (defaults to 1)

#### Vector creation examples

```javascript
const blas = require('blasjs');

const { fortranArrComplex64, fortranArrComplex32 } = blas.helper;

const complexDataArr = [
    { re: 1.8, im: -0.2 },
    { re: 2.3, im: 0.6 }
];

const realData = [ 0.1, 2, 0.34, .56 ];

const sp1 = fortranArrComplex32(complexDataArr)();
//sp1.r = [ 1.7999999523162842, 2.299999952316284 ],
//sp1.i = [ -0.20000000298023224, 0.6000000238418579 ],

const sp2 = fortranArrComplex32(realData)();
//sp2.r =  [ 0.10000000149011612, 2, 0.3400000035762787, 0.5600000023841858 ]
//sp2.i = undefined

const sp3 = fortranArrComplex32({re:0.2, im:-0.3})();
//[ 0.20000000298023224 ]
//[ -0.30000001192092896 ]

const sp4 = fortranArrComplex32(123)(4);
/*{
  base: 4,
  r: Float32Array [ 123 ],
  i: undefined,
}*/

const sdp1 = fortranArrComplex64(complexDataArr)();
//sp1.r = [ 1.8, 2.3 ],
//sp1.i = [ -0.2, 0.6 ],

const sdp2 = fortranArrComplex64(realData)();
//sp2.r =  [ 0.1, 2, 0.34, 0.56 ]
//sp2.i = undefined

const sp3 = fortranArrComplex64({re:0.2, im:-0.3})();
//[ 0.2 ]
//[ -0.3 ]

const sp4 = fortranArrComplex64(123)(4);
/*{
  base: 4,
  r: Float32Array [ 123 ],
  i: undefined,
}*/
```

## Matrix Constructors

These constructors create the [`Matrix`](#matrix) object for working with single/double precision complex/real Matrices.

### `fortranMatrixComplex32`

Constructs a [Matrix](#matrix) object using [Float32Array][float32-array] as the underlying array(s) (plural in the case of complex) elements.

```typescript
declare function fortranMatrixComplex32(...rest: (Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix
```

`Argument list`:

* `rest`: takes as input.
    * A single numeric value.
    * A single [`Complex`](#type-complex) object.
    * An array of [`Complex`](#type-complex) objects.
    * An array of number values.
* `nrRows`: where nrRows is equal to `n` in the matrix A(m,n).
* `nrCols`: where nrCols is equal to `m` in the matrix A(m,n).
* `rowBase`: FORTRAN offset for the first dimension (rows) as explained in [Language differences][language-differences].
* `rowBase`: FORTRAN offset for the second dimension (columns) as explained in [Language differences][language-differences].

See _[Examples](#matrix-creation-examples)_


### `fortranMatrixComplex64`

Constructs a [Matrix](#matrix) object using [Float64Array][float64-array] as the underlying array(s) (plural in the case of complex) elements.

```typescript
declare function fortranMatrixComplex64(...rest: (Complex | Complex[])[]):
    (nrRows: number, nrCols: number, rowBase?: number, colBase?: number) => Matrix
```

`Argument list`:

* `rest`: takes as input.
    * A single numeric value.
    * A single [`Complex`](#type-complex) object.
    * An array of [`Complex`](#type-complex) objects.
    * An array of number values.
*  `nrRows`: where rnRows is equal to `n` in the matrix A(m,n).
*  `nrCols`: where nrCols is equal to `m` in the matrix A(m,n).
*  `rowBase`: FORTRAN offset for the first dimension (rows) as explained in [Language differences][language-differences].
*  `rowBase`: FORTRAN offset for the second dimension (columns) as explained in [Language differences][language-differences].

### Matrix Creation Examples

```javascript
const blas = require('blasjs');
const {
    fortranMatrixComplex64,
    fortranMatrixComplex32
} = blas.helper;

// some matrix data 3x3 array  aka a_row_column

const a11 = { re: .2, im: -.11 };
const a21 = { re: .1, im: -.2 };
const a31 = { re: .3, im: .9 };
const a12 = { re: .4, im: .5 };
const a22 = { re: .9, im: -.34 };
const a32 = { re: -.2, im: .45 };
const a13 = { re: -.1, im: .89 };
const a23 = { re: .43, im: .23 };
const a33 = { re: .23, im: .56 };

const {
    fortranMatrixComplex64,
    fortranMatrixComplex32
} = blas.helper;

// Some matrix data 3x3 array  aka a_row_column

const a11 = { re: .2, im: -.11 };
const a21 = { re: .1, im: -.2 };
const a31 = { re: .3, im: .9 };
const a12 = { re: .4, im: .5 };
const a22 = { re: .9, im: -.34 };
const a32 = { re: -.2, im: .45 };
const a13 = { re: -.1, im: .89 };
const a23 = { re: .43, im: .23 };
const a33 = { re: .23, im: .56 };

//functional curry to prepare for different mappings of A()
const A32 = fortranMatrixComplex64([
    a11, a21, a31, a12, a22, a32, a13, a23, a33
]);

//matrix 1
const m1 = A32(3, 3); // 3x3 matrix with rowBase=1, colBase=1

// mimic FORTRAN  "COMPLEX*8  A(-2:1, -3:0)"
const m2 = A32(3, 3, -2, -3);

//same as FORTRAN default COMPLEX*8 A(3,3) !aka A(1:3,1:3)
const m3 = A32(3, 3, 1, 1)

/* double precision */
/* double precision */
/* double precision */

const A64 = fortranMatrixComplex64([
    a11, a21, a31, a12, a22, a32, a13, a23, a33
]);

// matrix 1 FORTRAN "COMPLEX*16  A(-2:1, -3:0).
const m1 = A64(3, 3); // 3x3 matrix with rowBase=1, colBase=1

// mimic FORTRAN  "COMPLEX*16  A(-2:1, -3:0)"
const m2 = A64(3, 3, -2, -3);

// same as FORTRAN default COMPLEX*16 A(3,3) !aka A(1:3,1:3)
const m3 = A64(3, 3, 1, 1);
```

# A note on numeric precision

In `blasjs`, contrary to the FORTRAN reference implementation, the numeric precision of a routine, is not determined by its name but by [*how*](#vector-constructors) its arguments like [`FortranArr`](#fortranarr) and [`Matrix`](#matrix) are constructed before used as arguments in `blasjs` routines. The original FORTRAN names are kept for backwards compatibility to ease the porting of FORTRAN code toward `blasjs`.

# Mimicking FORTRAN OUT Arguments


In FORTRAN a subroutine can have IN, OUT and IN/OUT scalar arguments. In JavaScript only arguments of type `object` are passed by reference. To mimic OUT and IN/OUT FORTRAN arguments, scalars are wrapped in a JS object. See [Construct a Givens plane rotation](#construct-a-givens-plane-rotation) for an example.


# Level 1 routines

Routines categorized as _Level 1_ perform scalar-vector and vector-vector operations.

## Euclidean norm: √(xᴴ·x)  or √(xᵀ·x)

Calculates the norm of a (complex) vector.

xᴴ is the _conjugate_ of x

xᵀ is the _transpose_ of x

#### Naming

* `scrnm2`: complex, [single or double precision][precision-note]. See [blas ref][ref-scnrm2].
* `dznrm2`: complex, (alias for `scrnm2`). See [blas ref][ref-dznrm2].
* `snrm2`: real, [single or double precision][precision-note]. See [blas ref][ref-snrm2].
* `dnrm2`: real, (alias for `dnrm2`). See [blas ref][ref-dnrm2].

_decl_

```typescript
function scnrm2(n: number, x: FortranArr, incx: number): number;
function dznrm2(n: number, x: FortranArr, incx: number): number;
function snrm2(n: number, x: FortranArr, incx: number): number;
function dnrm2(n: number, x: FortranArr, incx: number): number;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { scnrm2, dznrm2, snrm2, dnrm2 } = BLAS.level1;
```

## Construct a Givens plane rotation

See [wiki][givens-rotation].

```Math
 |c  -s| x |a| =   |r |
 |s   c|   |b|     |0 |

 r =  √( a² + b² )

```

#### Naming

* `srotg`: real, (alias for `drotg`). See [blas ref][ref-srotg].
* `drotg`: real, [single or double precision][precision-note]. See [blas ref][ref-drotg].
* `crotg`: complex, [single or double precision][precision-note]. See [blas ref][ref-crotg].
* `zrotg`: complex, (alias for `crotg`). See [blas ref][ref-zrotg].

_decl_

```typescript
function srotg(p: { sa: number, sb: number, c: number, s: number } ): void;
function drotg(p: { sa: number, sb: number, c: number, s: number } ): void;
function crotg(ca: Complex, cb: Complex, c: { val: number }, s: Complex ): void
function zrotg(ca: Complex, cb: Complex, c: { val: number }, s: Complex ): void
```

Usage:

```javascript
const BLAS = require('blasjs');
const { srotg, drotg, crotg, zrotg } = BLAS.level1;
```

## Construct the **modified** Givens rotation matrix `H`

Construct the modified Givens transformation matrix H which zeros
the second component of the 2 vector  ( sx1*√(sd1) , sy1* √(sd2) )
See [researchgate.net][construct-modified-givens-transformation].

#### Naming

* `srotmg`: real, (alias for `drotmg`). See [blas ref][ref-srotmg].
* `drotmg`: real, [single or double precision][precision-note]. See [blas ref][ref-drotmg].

_decl_

```typescript
function srotmg(p: { sd1: number, sd2: number, sx1: number, sy1: number, sparam: FortranArr }): void;
function drotmg(p: { sd1: number, sd2: number, sx1: number, sy1: number, sparam: FortranArr }): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { srotmg, drotmg } = BLAS.level1;
```

## Apply the modified Givens Transformation

See [wiki][apply-modified-givens-transformation].

#### Naming

* `srotm`: real, (alias for `drotm`). See [blas ref][ref-srotm].
* `drotm`: real, [single or double precision][precision-note]. See [blas ref][ref-drotm].

_decl_

```typescript
function srotm(n: number, sy: FortranArr, incx: number, sy: FortranArr, incy: number, sparam: FortranArr)): void;

function drotm(n: number, sy: FortranArr, incx: number, sy: FortranArr, incy: number, sparam: FortranArr)): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { srotm, drotm } = BLAS.level1;
```

## Applies a plane rotation

See [researchgate.net][construct-modified-givens-transformation].

#### Naming

* `srot`: real, (alias for `drot`). See [blas ref][ref-srot].
* `drot`: real, [single or double precision][precision-note]. See [blas ref][ref-drot].
* `csrot`: complex, (alias for `zdrot`). See [blas ref][ref-csrot].
* `zdrot`: complex, [single or double precision][precision-note]. See [blas ref][ref-zdrot].

_decl_

```typescript
function srot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number, c: number, s: number): void;
function drot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number, c: number, s: number): void;

function csrot: (n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number, c: number, s: number): void;
function zdrot: (n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number, c: number, s: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { srot, drot, csrot, zdrot } = BLAS.level1;
```

## Scale a vector by a constant

x ⟵ α·x

#### Naming

* `sscal`: Alias for `dscal`. See [blas ref][ref-sscal].
* `dscal`:  by a REAL constant. See [blas ref][ref-dscal].
* `cscal`: Alias for `zscal`. See [blas ref][ref-cscal].
* `zscal`: Scales a COMPLEX vector with a COMPLEX constant. See [blas ref][ref-zscal].
* `csscal`: Alias for `zdscal`. [blas ref][ref-csscal].
* `zdscal`: Scales a COMPLEX vector with a REAL constant. See [blas ref][ref-zdscal].

_decl_

```typescript
function sscal(n: number, sa: number, sx: FortranArr, incx: number): void;
function dscal(n: number, sa: number, sx: FortranArr, incx: number): void;
function cscal(n: number, ca: Complex,cx: FortranArr, incx: number): void;
function zscal(n: number, ca: Complex,cx: FortranArr, incx: number): void;
function csscal(n: number, sa: number, cx: FortranArr, incx: number): void;
function zdscal(n: number, sa: number, cx: FortranArr, incx: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sscal, dscal, cscal, zscal, csscal, zdscal } = BLAS.level1;
```

## Takes the sum of the absolute values of the components of vector

s ⟵ ∑ ∥ Re( x ) ∥ + ∥ Im( x ) ∥

#### Naming

* `sasum`: Alias for `dasum`. See [blas ref][ref-sasum]
* `dasum`: uses REAL vector, ( [single or double precision][precision-note] ). See [blas-ref][ref-dasum].
* `scasum`: Alias for `dzasum`. See [blas ref][ref-scasum].
* `dzasum`: uses Complex vector, ( [single or double precision][precision-note] ). See [blas-ref][ref-dzasum].

_decl_

```typescript
function sasum(n: number, sx: FortranArr, incx: number): number;
function dasum(n: number, sx: FortranArr, incx: number): number;
function scasum(n: number, cx: FortranArr, incx: number): number;
function dzasum(n: number, cx: FortranArr, incx: number): number;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sasum, dasum, scasum, dzasum } = BLAS.level1;
```

## Interchanges 2 vectors

Swap 2 vectors.

#### Naming

* `sswap`: Alias for `dswap`. See [blas ref][ref-sswap].
* `dswap`: REAL vector, ( [single or double precision][precision-note] ). See [blas ref][ref-dasum].
* `cswap`: Alias for `zswap`. See [blas ref][ref-xcswap].  
* `zswap`: REAL vector, ( [single or double precision][precision-note] ). See [blas ref][ref-zswap].

_decl_

```javascript
function sswap(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number ): void;
function dswap(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number ): void;
function cswap(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number ): void;
function zswap(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number ): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sswap, dswap, cswap, zswap } = BLAS.level1;
```

## Dot product of two complex vectors

  xᵀ·y or xᴴ·y

### Naming

* `cdotu`: Alias for `zdotu`. See [blas ref][ref-zdotc].
* `cdotc`: Alias for `zdotc`. See [blas ref][ref-cdotc].
* `zdotu`: `xᵀ·y`. Complex arguments, ( [single or double precision][precision-note] ). See [blas-ref][ref-zdotu].
* `zdotc`: `xᴴ·y`. The fist complex vector argument is made conjugate, ( [single or double precision][precision-note] ). See [blas-ref][ref-zdotc].

_decl_

```javascript
function cdotu(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number): Complex;

// first argument sx is made conjugate
function cdotc(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number ): Complex;

function zdotu(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number ): Complex;

// first argument sx is made conjugate
function zdotc(n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number ): Complex;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { cdotu, cdotc, zdotu, zdotc } = BLAS.level1;
```

## Dot product of two non complex vectors
  
xᵀ·y

### Naming

* `sdot`: Alias for `dsdot`. See [blas ref][ref-sdot].
* `ddot`: Alias for `dsdot`. See [blas ref][ref-ddot].
* `sdsdot`: Alias for `dsdot`. See [blas ref][ref-sdsdot].
* `dsdot`: `xᵀ·y` Inner product of 2 vectors ( [single or double precision][precision-note] ). See [blas ref][ref-dsdot].

_decl_

```javascript
function sdot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): number;
function ddot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): number;
function sdsdot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): number;
function dsdot(n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): number;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sdot, ddot, sdsdot, dsdot } = BLAS.level1;
```

## Finds the index of the first element having maximum absolut value.

Find k for wich: ∥ xₖ ∥ > ∥ xₜ ∥ for all t ∈ [1, n].

### Naming


* `isamax`: Alias for `idamax`. See [blas ref]:[ref-isamax]
* `idamax`: Find the index of the maximum element of a REAL vector ( [single or double precision][precision-note] ). See [blas ref][ref-idamax]. 
* `icamax`: Alias for `izamax`. See [blas ref]:[ref-icamax]
* `izamax`: Find the index of the maximum element of a COMPLEX vector ( [single or double precision][precision-note] ). See [blas ref][ref-izamax].

_decl_

```javascript
function isamax: (n: number, sx: FortranArr, incx: number): number;
function idamax: (n: number, sx: FortranArr, incx: number): number;
function icamax: (n: number, sx: FortranArr, incx: number): number;
function izamax: (n: number, sx: FortranArr, incx: number): number;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { isamax, idamax, icamax, izamax } = BLAS.level1;
```

## Copy a vector x to a vector y

### Naming

* `scopy`: Alias for `dcopy`. See [blas ref]:[ref-scopy]
* `dcopy`: Copies a REAL vector ( [single or double precision][precision-note] ). See [blas ref][ref-dcopy].
* `ccopy`: Alias for `zcopy`. See [blas ref]:[ref-ccopy]
* `zcopy`: Copies a COMPLEX vector ( [single or double precision][precision-note] ). See [blas ref][ref-zcopy].

_decl_

```javascript
function scopy (n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): void;
function dcopy (n: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): void;
function ccopy (n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number): void;
function zcopy (n: number, cx: FortranArr, incx: number, cy: FortranArr, incy: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { scopy, dcopy, ccopy, zcopy } = BLAS.level1;
```

## Constant times a vector plus a vector

y ⟵ y + a·x  where y, a and x can be complex or a real number.

### Naming

* `saxpy`: Alias for `daxpy`. See [blas ref]:[ref-saxpy].
* `daxpy`: REAL constant used in multiplication with a vector ( [single or double precision][precision-note] ). See [blas ref]:[ref-daxpy].   
* `caxpy`: Alias for `zaxpy`. See [blas ref]:[ref-saxpy].
* `zaxpy`: Complex constant used in multiplication with a vector ( [single or double precision][precision-note] ). See [blas ref]:[ref-zaxpy].

_decl_

```javascript
function saxpy(n: number, sa: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): void;
function daxpy(n: number, sa: number, sx: FortranArr, incx: number, sy: FortranArr, incy: number): void;
function caxpy(n: number, ca: Complex, cx: FortranArr, incx: number, cy: FortranArr, incy: number): void;
function zaxpy(n: number, ca: Complex, cx: FortranArr, incx: number, cy: FortranArr, incy: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { saxpy, daxpy, caxpy, zaxpy } = BLAS.level1;
```

# Level 2 Routines

Routines categorized as _Level 2_ perform Matrix-vector operations.

## The hermitian rank 2 operation A ⟵ α·x·yᴴ + conjg( α )·y·xᴴ + A  

( ᴴ means conjugate transpose )

For the routines `chpr2` and `zhpr2` the matrix A is in packed form ( a [fortranArr](#vector-constructors) ).

For the routines `cher2` and `zher2` the matrix symmetry is exploited (use only upper/lower triangular part of the matrix). 

### Naming

* `cher2`: alias for `zher2`. See [blas ref][ref-cher2].
* `zher2`: The Matrix `A` is in upper or lower triangular form ( [single or double precision][precision-note] ). See [blas ref][ref-zher2].  
* `chpr2`: alias for `zhpr2`. See [blas ref][ref-chpr2].
* `zhpr2`: The matrix `A` is in [packed](#packed-matrices) form ( [single or double precision][precision-note] ). See [blas ref][ref-zhpr2].

_decl_

```typescript
function cher2(
     uplo: "u" | "l", 
     n: number, 
     alpha: Complex, 
     x: FortranArr, 
     incx: number, 
     y: FortranArr, 
     incy: number, 
     a: Matrix, 
     lda: number): void;

function zher2(
     uplo: "u" | "l", 
     n: number, 
     alpha: Complex, 
     x: FortranArr, 
     incx: number, 
     y: FortranArr, 
     incy: number, 
     a: Matrix, 
     lda: number): void; 

function chpr2(
     uplo: "u" | "l", 
     n: number, 
     alpha: Complex, 
     x: FortranArr, 
     incx: number, 
     y: FortranArr, 
     incy: number, 
     ap: FortranArr): void;

function zhpr2(
     uplo: "u" | "l", 
     n: number, 
     alpha: Complex, 
     x: FortranArr, 
     incx: number, 
     y: FortranArr, 
     incy: number, 
     ap: FortranArr): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

See: _[how to create Matrix](#matrix-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { cher2, zher2, chpr, zhpr } = BLAS.level2;
```

## The symmetric rank 2 operation A ⟵ α·x·yᵀ + α·y·xᵀ + A

For the routines `sspr2` and `dspr2` the matrix A is in [packed](#packed-matrices) form ( a [fortranArr](#vector-constructors) ).

For the routines `ssyr2` and `dsyr2` the matrix symmetry is exploited (use only upper/lower triangular part of the matrix). 

### Naming

* `sspr2`: Alias for dspr2. See [blas ref][ref-sspr2].
* `dspr2`: The matrix `A` is in [packed](#packed-matrices) form ( [single or double precision][precision-note] ). See [blas ref][ref-dspr2].
* `ssyr2`: Alias for dsyr2. See [blas ref][ref-ssyr2].
* `dsyr2`: The Matrix `A` is in upper or lower triangular form ( [single or double precision][precision-note] ). See [blas ref][ref-dsyr2].

 _decl_

```typescript
function sspr2(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    ap: FortranArr):void;

function dspr2(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    ap: FortranArr):void;

function ssyr2(
    uplo: 'u' | 'l', 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    A: Matrix, 
    lda: number): void;

function dsyr2(
    uplo: 'u' | 'l', 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    A: Matrix, 
    lda: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

See: _[how to create Matrix](#matrix-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sspr2, dspr2, ssyr2, dsyr2 } = BLAS.level2;
```

## The rank 1 operation A ⟵ α·x·yᴴ + A or A ⟵ α·x·yᵀ + A

( ᴴ means conjugate transpose )

The subroutines `sger` and `dger` perform  A ⟵ α·x·yᵀ + A. Where α is a REAL scalar,
 A, x, y are [single or double precision][precision-note] REAL Matrix and vectors.

The subroutines `cgerc` and `zgerc` perform  A ⟵ α·x·yᴴ + A. Where α is a COMPLEX scalar,
 A, x, y are [single or double precision][precision-note] COMPLEX Matrix and vectors.

The subroutines `cgeru` and `zgeru` perform  A ⟵ α·x·yᵀ + A. Where α is a COMPLEX scalar,
 A, x, y are [single or double precision][precision-note] COMPLEX Matrix and vectors.  

### Naming

* `sger`: alias for `dger`. See [blas ref][ref-sger].
* `dger`: See [blas ref][ref-dger].
* `cgerc`: alias for `zgerc`. See [blas ref][ref-cgerc].
* `zgerc`: See [blas ref][ref-zgerc].
* `cgeru`: alias for `zgeru`. See [blas ref][ref-cgeru].
* `zgeru`: See [blas ref][ref-zgeru].

_decl_

```typescript
function sger(
    m: number, 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number):void;

function dger(
    m: number, 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number):void;
    
function cgerc(
    m: number, 
    n: number, 
    alpha: Complex, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number): void;

function zgerc(
    m: number, 
    n: number, 
    alpha: Complex, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number): void;

function cgeru(
    m: number, 
    n: number, 
    alpha: Complex, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number): void;
    
function zgeru(
    m: number, 
    n: number, 
    alpha: Complex, 
    x: FortranArr, 
    incx: number, 
    y: FortranArr, 
    incy: number, 
    a: Matrix, 
    lda: number): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

See: _[how to create Matrix](#matrix-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sger, dger, cgerc, zgerc, cgeru, zgeru } = BLAS.level2;
```

## The hermitian rank 1 operation A ⟵ α·x·xᴴ + A

( ᴴ means conjugate transpose )

For the routines `cher` and `zher` α is a REAL scalar, the matrix symmetry of A is exploited (use only upper/lower triangular part of the matrix).

For the routines `chpr` and `zhpr` α is a REAL scalar, the matrix A is in [packed](#packed-matrices) form ( a [fortranArr](#vector-constructors) ).

### Naming

* `cher`: alias for `zher`. See [blas ref][ref-cher].
* `zher`: For [single or double precision][precision-note] complex `x` and `A`. See [blas ref][ref-cher].
* `chpr`: alias for `zher`. See [blas ref][ref-cher].
* `zhpr`: For [single or double precision][precision-note] complex `x` and `A`. See [blas ref][ref-cher].


```typescript
function cher(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    a: Matrix, 
    lda: number): void;

function zher(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    a: Matrix, 
    lda: number): void;

function chpr(u
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    ap: FortranArr): void;

function zhpr(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    ap: FortranArr): void;
```

See: _[how to create fortranArr](#vector-constructors)_.

See: _[how to create Matrix](#matrix-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { cher, zher, chpr, zhpr } = BLAS.level2;
```

## The symmetric rank 1 operation A ⟵ α·x·xᵀ + A

For the routines `ssyr` and `dsyr` α is a REAL scalar, the symmetry of the REAL matrix A is exploited (use only upper/lower triangular part of the matrix).

For the routines `sspr` and `dspr` α is a REAL scalar, the REAL matrix A is in [packed](#packed-matrices) form ( a [fortranArr](#vector-constructors) ).

### Naming

* `sspr`: alias for `dspr`. See [blas ref][ref-sspr].
* `dspr`: For [single or double precision][precision-note] REAL `α` , `x` and `A`. See [blas ref][ref-dspr].
* `ssyr`: alias for `ssyr`. See [blas ref][ref-ssyr].
* `dsyr`: For [single or double precision][precision-note] REAL `α` , `x` and `A`. See [blas ref][ref-dsyr].

_decl_

```typescript
function sspr(
   uplo: "u" | "l", 
   n: number, 
   alpha: number, 
   x: FortranArr, 
   incx: number, 
   ap: FortranArr): void;

function dspr(
   uplo: "u" | "l", 
   n: number, 
   alpha: number, 
   x: FortranArr, 
   incx: number, 
   ap: FortranArr): void;

function ssyr(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    a: Matrix, 
    lda: number): void;

function dsyr(
    uplo: "u" | "l", 
    n: number, 
    alpha: number, 
    x: FortranArr, 
    incx: number, 
    a: Matrix, 
    lda: number): void;    
```

See: _[how to create fortranArr](#vector-constructors)_.

See: _[how to create Matrix](#matrix-constructors)_.

Usage:

```javascript
const BLAS = require('blasjs');
const { sspr, dspr, ssyr, dsyr } = BLAS.level2;
```

## The matrix-vector operation, y ⟵ α·A·x + β·y, or y ⟵ α·Aᵀ·x + β·y or y ⟵  α·Aᴴ·x + β·y




#### Naming

* `cgbmv`:
* `chbmv`:
* `dgbmv`:
* `dsbmv`:
* `sgbmv`:
* `ssbmv`:
* `zgbmv`:
* `zhbmv`:
* `ctbmv`:
* `dtbmv`:
* `stbmv`:
* `ztbmv`:
* `cgemv`:
* `chemv`:
* `dgemv`:
* `sgemv`:
* `zgemv`:
* `zhemv`:
* `chpmv`:
* `dspmv`:
* `sspmv`:
* `zhpmv`:
* `dsymv`:
* `ssymv`:

[srotg]: https://en.wikipedia.org/wiki/Givens_rotation
[givenmodified]: https://www.ibm.com/support/knowledgecenter/en/SSFHY8_5.5.0/com.ibm.cluster.essl.v5r5.essl100.doc/am5gr_srotm.htm

[caxpy]:  http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html

[blas-site]: http://www.netlib.org/blas/
[blas-source]: https://github.com/Reference-LAPACK/lapack/tree/master/BLAS
[float32-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float32Array]

[float64-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float64Array]

[mimic-fortran-args]: [#mimicking-FORTRAN-OUT-Arguments]
[precision-note]: [#a-note-on-numeric-precision]
[language-differences]: [#language-differences-with-fortranblas]
[givens-rotation]: https://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation
[apply-givens-rotation]: https://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation
[ref-snrm2]: http://www.netlib.org/lapack/explore-html/d7/df1/snrm2_8f.html
[ref-dnrm2]: http://www.netlib.org/lapack/explore-html/da/d7f/dnrm2_8f.html
[ref-scnrm2]: http://www.netlib.org/lapack/explore-html/db/d66/scnrm2_8f.html
[ref-srotg]: http://www.netlib.org/lapack/explore-html/d7/d26/srotg_8f.html
[ref-zrotg]: http://www.netlib.org/lapack/explore-html/dc/dfe/zrotg_8f.html

[ref-srotmg]: http://www.netlib.org/lapack/explore-html/dd/d48/srotmg_8f.html
[ref-drotmg]: http://www.netlib.org/lapack/explore-html/df/deb/drotmg_8f.html

[construct-modified-givens-transformation]: https://www.researchgate.net/profile/JV_Mccanny/publication/224312422_Modified_Givens_rotations_and_their_application_to_matrix_inversion/links/55cdcefa08aee19936f80088/Modified-Givens-rotations-and-their-application-to-matrix-inversion.pdf

[ref-srotm]: http://www.netlib.org/lapack/explore-html/d6/d0f/srotm_8f.html
[ref-drotm]: http://www.netlib.org/lapack/explore-html/d8/d7b/drotm_8f.html

[ref-srot]: http://www.netlib.org/lapack/explore-html/db/d6c/srot_8f.html
[ref-drot]: http://www.netlib.org/lapack/explore-html/dc/d23/drot_8f.html
[ref-csrot]: http://www.netlib.org/lapack/explore-html/d1/dbb/csrot_8f.html
[ref-zdrot]: http://www.netlib.org/lapack/explore-html/d4/de9/zdrot_8f.html

[ref-sscal]: http://www.netlib.org/lapack/explore-html/d9/d04/sscal_8f.html
[ref-dscal]: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
[ref-cscal]: http://www.netlib.org/lapack/explore-html/dc/d81/cscal_8f.html
[ref-zscal]: http://www.netlib.org/lapack/explore-html/d2/d74/zscal_8f.html
[ref-csscal]: http://www.netlib.org/lapack/explore-html/de/d5e/csscal_8f.html
[ref-zdscal]: http://www.netlib.org/lapack/explore-html/dd/d76/zdscal_8f.html

[ref-sasum]: http://www.netlib.org/lapack/explore-html/df/d1f/sasum_8f.html
[ref-dasum]: http://www.netlib.org/lapack/explore-html/de/d05/dasum_8f.html
[ref-scasum]: http://www.netlib.org/lapack/explore-html/db/d53/scasum_8f.html
[ref-dzasum]: http://www.netlib.org/lapack/explore-html/df/d0f/dzasum_8f.html

[ref-cher2]: http://www.netlib.org/lapack/explore-html/db/d87/cher2_8f.html
[ref-zher2]: http://www.netlib.org/lapack/explore-html/da/d8a/zher2_8f.html
[ref-chpr2]: http://www.netlib.org/lapack/explore-html/d6/d44/chpr2_8f.html
[ref-zhpr2]: http://www.netlib.org/lapack/explore-html/d5/d52/zhpr2_8f.html

[ref-sspr2]: http://www.netlib.org/lapack/explore-html/db/d3e/sspr2_8f.html
[ref-dspr2]: http://www.netlib.org/lapack/explore-html/dd/d9e/dspr2_8f.html
[ref-ssyr2]: http://www.netlib.org/lapack/explore-html/db/d99/ssyr2_8f.html
[ref-dsyr2]: http://www.netlib.org/lapack/explore-html/de/d41/dsyr2_8f.html 

[ref-sger]: http://www.netlib.org/lapack/explore-html/db/d5c/sger_8f.html
[ref-dger]: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
[ref-cgerc]: http://www.netlib.org/lapack/explore-html/dd/d84/cgerc_8f.html
[ref-zgerc]: http://www.netlib.org/lapack/explore-html/d3/dad/zgerc_8f.html
[ref-cgeru]: http://www.netlib.org/lapack/explore-html/db/d5f/cgeru_8f.html
[ref-zgeru]: http://www.netlib.org/lapack/explore-html/d7/d12/zgeru_8f.html

[ref-cher]: http://www.netlib.org/lapack/explore-html/d3/d6d/cher_8f.html
[ref-zher]: http://www.netlib.org/lapack/explore-html/de/d0e/zher_8f.html
[ref-chpr]: http://www.netlib.org/lapack/explore-html/db/dcd/chpr_8f.html
[ref-zhpr]: http://www.netlib.org/lapack/explore-html/de/de1/zhpr_8f.html
 
[ref-sspr]: http://www.netlib.org/lapack/explore-html/d2/d9b/sspr_8f.html
[ref-dspr]: http://www.netlib.org/lapack/explore-html/dd/dba/dspr_8f.html
[ref-ssyr]: http://www.netlib.org/lapack/explore-html/d6/dac/ssyr_8f.html
[ref-dsyr]: http://www.netlib.org/lapack/explore-html/d3/d60/dsyr_8f.html 

 
[ref-cgbmv]: http://www.netlib.org/lapack/explore-html/d0/d75/cgbmv_8f.html
[ref-chbmv]: http://www.netlib.org/lapack/explore-html/db/dc2/chbmv_8f.html
[ref-dgbmv]: http://www.netlib.org/lapack/explore-html/d2/d3f/dgbmv_8f.html
[ref-dsbmv]: http://www.netlib.org/lapack/explore-html/d8/d1e/dsbmv_8f.html
[ref-sgbmv]: http://www.netlib.org/lapack/explore-html/d6/d46/sgbmv_8f.html
[ref-ssbmv]: http://www.netlib.org/lapack/explore-html/d3/da1/ssbmv_8f.html
[ref-zgbmv]: http://www.netlib.org/lapack/explore-html/d9/d46/zgbmv_8f.html
[ref-zhbmv]: http://www.netlib.org/lapack/explore-html/d3/d1a/zhbmv_8f.html
[ref-ctbmv]: http://www.netlib.org/lapack/explore-html/d3/dcd/ctbmv_8f.html
[ref-dtbmv]: http://www.netlib.org/lapack/explore-html/df/d29/dtbmv_8f.html
[ref-stbmv]: http://www.netlib.org/lapack/explore-html/d6/d7d/stbmv_8f.html
[ref-ztbmv]: http://www.netlib.org/lapack/explore-html/d3/d39/ztbmv_8f.html
[ref-cgemv]: http://www.netlib.org/lapack/explore-html/d4/d8a/cgemv_8f.html
[ref-chemv]: http://www.netlib.org/lapack/explore-html/d7/d51/chemv_8f.html
[ref-dgemv]: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
[ref-sgemv]: http://www.netlib.org/lapack/explore-html/db/d58/sgemv_8f.html 
[ref-zgemv]: http://www.netlib.org/lapack/explore-html/db/d40/zgemv_8f.html
[ref-zhemv]: http://www.netlib.org/lapack/explore-html/d0/ddd/zhemv_8f.html
[ref-chpmv]: http://www.netlib.org/lapack/explore-html/d3/d1a/zhbmv_8f.html
[ref-dspmv]: http://www.netlib.org/lapack/explore-html/d4/d85/dspmv_8f.html
[ref-sspmv]: http://www.netlib.org/lapack/explore-html/d8/d68/sspmv_8f.html
[ref-zhpmv]: http://www.netlib.org/lapack/explore-html/d0/d60/zhpmv_8f.html
[ref-dsymv]: http://www.netlib.org/lapack/explore-html/d8/dbe/dsymv_8f.html
[ref-ssymv]: http://www.netlib.org/lapack/explore-html/d2/d94/ssymv_8f.html
        


