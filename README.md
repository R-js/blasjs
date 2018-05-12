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
    * [General Helpers](#general-helpers)
        * [arrayrify](#arrayrify)
        * [complex](#complex)
        * [each](#each)
        * [map](#map)
        * [muxCmplx](#muxCmplx)
        * [numberPrecision](#numberPrecision)
    * [Blas array/matrix constructors](#blas-arraymatrix-constructors)
        * [fortranArrComplex32](#fortranArrComplex32)
        * [fortranArrComplex64](#fortranArrComplex64)
        * [](#)
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
declare type Matrix = { //Incomplete decleration
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
    colSize: number;
    .
    colOfEx(number): number;
    coord(col): (row) => number;
    setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
    .
}
```

Explanation:

* `colSize`: The number of rows in the matrix.
* `nrCols`: The number of columns in the matrix.
* `colofEx`: Calculates the physical location of a `column offset` within the `TypedArray`. Taking int account the column base `colBase` and row base `colBase`. The index of  A(i,j) `=  (j - colBase)*colSize + i - rowBase`.
* `coord`: Curried, emulates non-zero based FORTRAN index values for 2 dimensional Arrays. The index that is iterated over the least (usually) is used as the first to create the curried function.
* `setCol`: Uses underlying `TypedArray`, `fill` method to set multiple column elements to a single value.

Examples:

This code sample will create a Matrix object encapsulating a 3x3 complex array.

```javascript

const blas = require('blasjs');
const { fortranMatrixComplex64, muxComplx } = blas.util;

/**
 * Create a Martrix object and fill it with some data
*/
const a11 = { re: .2, im: -.11 };
const a21 = { re: .1, im: -.2 };
const a31 = { re: .3, im: .9 };
const a12 = { re: .4, im: .5 };
const a22 = { re: .9, im: -.34 };
const a32 = { re: -.2, im: .45 };
const a13 = { re: -.1, im: .89 };
const a23 = { re: .43, im: .23 };
const a33 = { re: .23, im: .56 };

const a = fortranMatrixComplex64(muxComplx(
   a11, a21, a31, a12, a22, a32, a13, a23 a33
]))(3,3);
/*
   The data in Matrix "a" is now stored in 2 TypedArray's.
*/

// Get the second column offset
// formula: "(j - colBase )* colSize - rowBase"
const column3 = a.colOfEx(3);
//5

a13.re === A.r[column3+1]; // Direct access to the TypedArray in A (real part)
//true

a23.im === A.i[column3+2]; // Direct access to the TypedArray in A (imaginary part)
//true

// use Fortran like indexes (taking into account "rowBase" and "colBase")
//
A.r[ A.coord(2)(1) ] === a12.re
//true
```

#### Creating new transformed Matrix instances from existing ones

One can create/transform new Matrix instances form existing onces. A copy of all relevant data is made into the new `Matrix` instance.

#### `Matrix.prototype.slice`

Slices a rectangular piece of data out of an matrix into a new `Matrix` instance. **All arguments are FORTRAN-style non-zero based indexes**.

```typescript
declare type Matrix = { ..
   slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
   .. }
```

* `rowStart`: The row in the matrix to begin slicing.
* `rowEnd`: The last row to include in the slice.
* `colStart`: The column in the matrix to begin slicing.
* `colEnd`: The last column to include in the slice.

_[See Example](#example-new-matrix)_

#### `Matrix.prototype.setLower`

Returns a new Matrix where everything below the matrix diagonal is set to a `value`.
Sets the real (and imaginary part, if it exist) to said value.

```typescript
declare type Matrix = { ..
  setLower(value?: number): Matrix;
   .. }
```

_[See Example](#example-new-matrix)_

#### `Matrix.prototype.setUpper`

Returns a new Matrix where everything below the matrix diagonal is set to a `value`.
Sets the real (and imaginary part, if it exist) to said value.

```typescript
declare type Matrix = { ..
  setLower(value?: number): Matrix;
   .. }
```

```typescript
 slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
    setLower(value?: number): Matrix;
    setUpper(value?: number): Matrix;
    upperBand(value: number): Matrix;
    lowerBand(value: number): Matrix;
    real(): Matrix;
    imaginary(): Matrix;


```typescript
declare Matrix {
    .
    setLower(value = m): Matrix;
    setUpper(value = m): Matrix;
    upperBand(value: number): Matrix;
    lowerBand(value: number): Matrix;
    .
}
```


The full type declaration of `Martrix`:

```typescript
interface Matrix {
    readonly rowBase: number;
    readonly colBase: number;
    readonly nrCols: number;
    readonly colSize: number;
    readonly r: fpArray;
    readonly i?: fpArray;
    readonly colOfEx: (number) => number;
    coord(col): (row) => number;
    setCol(col: number, rowStart: number, rowEnd: number, value: number): void;
    
    slice(rowStart: number, rowEnd: number, colStart: number, colEnd: number): Matrix;
    setLower(value?: number): Matrix;
    setUpper(value?: number): Matrix;
    upperBand(value: number): Matrix;
    lowerBand(value: number): Matrix;
    real(): Matrix;
    imaginary(): Matrix;
 
    packedUpper(value?: number): FortranArr;
    packedLower(value?: number): FortranArr;
    toArr(): Complex[] | number[];
}
```

* `rowBase`:




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
//world : hello
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

## Vector and Matrix constructors

JS Analog for working with single/double precision complex/real Arrays and Matrices. 

### `fortranArrComplex32`

Constructs a [FortranArr](#fortranArr) object using 2 [Float32Array][float32-array] as a single array of complex numbers. If only REAL numbers are
specified only a single `Float32Array` object will be used.

_decl_:

```typescript
function fortranArrComplex32(
    ...rest: (number | number[] | Complex | Complex[])[]
    ): (offset = 1) => FortranArr;
```

* `rest`: takes as input an array or a single value of type number or [Complex](#type-complex).
* `offset`: the Fortran offset (defaults to 1)
* returns an object of type [FortranArr](#fortranarr)

Usage:

```javascript
const blas = require('blasjs');

const { fortranArrComplex32 } = helper;

// You can also use the helper "complex" or "muxComplex"
// to generate JS complex arrays
const complexDataArr = [
    { re: 1.8, im: -0.2 },
    { re: 2.3, im: 0.6 }
];

// Create an object that mimics FORTRAN COMPLEX*8 SP(1:2)
//    and fill it with above data
const sp = fortranArrComplex32(complexArr)();
```

### `fortranArrComplex64`

Encapsulates two [Float64Array][float64-array] as a single array of complex numbers. If only REAL numbers are
used only a single `Float64Array` object will be used.

_decl_:

```typescript
function fortranArrComplex64(
    ...rest: (number | number[] | Complex | Complex[])[]
    ): (offset = 1) => FortranArr;
```

* `rest`: takes as input an array or a single value of type number or [Complex](#type-complex).
* `offset`: the Fortran offset (defaults to 1)
* returns an object of type (FortranArr)[#fortranarr]

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

// Create an object that mimics FORTRAN COMPLEX*16 SP(1:2)
//    and fill it with above data
const sp = fortranArrComplex64(complexArr)();
```

### `fortranMatrixComplex32`

Constructs an analog for a single precision 2-dimensional FORTRAN array (a matrix).

_decl_

```typescript

```

### `fortranMatrixComplex64`

Constructs an analog for a double precision 2-dimensional FORTRAN array (a matrix).


[srotg]: https://en.wikipedia.org/wiki/Givens_rotation
[givenmodified]: https://www.ibm.com/support/knowledgecenter/en/SSFHY8_5.5.0/com.ibm.cluster.essl.v5r5.essl100.doc/am5gr_srotm.htm

[caxpy]:  http://www.netlib.org/lapack/explore-html/da/df6/group__complex__blas__level1_ga9605cb98791e2038fd89aaef63a31be1.html

[blas-site]: http://www.netlib.org/blas/
[blas-source]: https://github.com/Reference-LAPACK/lapack/tree/master/BLAS
[float32-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float32Array]

[float64-array]: [https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Float64Array]

Some notes on Matrix symbols

https://fsymbols.com/generators/overline/

A̅ᵗ   Conjugate transpose
B̅    Conjugate

_`A̅ᵗ ∙ B̅`_

Aᵗ 
```
