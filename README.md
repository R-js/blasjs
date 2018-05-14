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
* Helper functions to easy porting of FORTRAN BLAS usage to Javascript.

#### Node and Web

The library js file is an agnostic UMD library, it can be used in a web client
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
        * [Matrix Creation Examples](###-matrix-creation-examples)
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

const { fortranArrComplex64 } = blas.helper;

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
*  `nrRows`: where rnRows is equal to `n` in the matrix A(m,n).
*  `nrCols`: where nrCols is equal to `m` in the matrix A(m,n).
*  `rowBase`: FORTRAN offset for the first dimension (rows) as explained in [Language differences][language-differences].
*  `rowBase`: FORTRAN offset for the second dimension (columns) as explained in [Language differences][language-differences].

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



[srotg]: https://en.wikipedia.org/wiki/Givens_rotation
[givenmodified]: https://www.ibm.com/support/knowledgecenter/en/SSFHY8_5.5.0/com.ibm.cluster.essl.v5r5.essl100.doc/am5gr_srotm.htm

[caxpy]:  http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html

[blas-site]: http://www.netlib.org/blas/
[blas-source]: https://github.com/Reference-LAPACK/lapack/tree/master/BLAS
[float32-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float32Array]

[float64-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float64Array]

[language-differences]: [#language-differences-with-fortranblas]

Some notes on Matrix symbols

https://fsymbols.com/generators/overline/

A̅ᵗ   Conjugate transpose
B̅    Conjugate

_`A̅ᵗ ∙ B̅`_

Aᵗ 
```
