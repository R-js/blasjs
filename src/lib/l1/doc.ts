const desc = {
    CAXPY: 'constant times a vector plus a vector.',
    ZAXPY: 'constant times a vector plus a vector.',
    DAXPY: 'constant times a vector plus a vector.',
    SAXPY: 'constant times a vector plus a vector.',
    //
    CCOPY: 'copies a vector x to a vector y.',
    ZCOPY: 'copies a vector, x, to a vector, y.',
    SCOPY: 'copies a vector, x, to a vector, y.',
    DCOPY: 'copies a vector, x, to a vector, y.',

    CDOTC: 'forms the dot product of two complex vectors',
    ZDOTC: 'forms the dot product of two complex vectors.',
    ZDOTU: 'forms the dot product of two complex vectors,  ZDOTU = X^T * Y',
    CDOTU: 'forms the dot product of two complex vectors',
    SDOT: 'forms the dot product of two vectors.',
    SDSDOT: 'Compute the inner product of two vectors with extended precision accumulation.',
    DSDOT: 'Compute the inner product of two vectors with extended precision accumulation and result.',
    DDOT: 'forms the dot product of two vectors.',

    CSCAL: 'scales a complex vector by a complex constant.',
    ZSCAL: 'scales a cplx vector by a cplx constant.',
    CSSCAL: 'scales a complex vector by a real constant.',
    ZDSCAL: 'scales a vector by a constant.',
    SSCAL: 'scales a vector by a constant.',
    DSCAL: 'scales a vector by a constant.',

    CSWAP: 'interchanges two vectors.',
    ZSWAP: 'interchanges two vectors.',
    SSWAP: 'interchanges two vectors.',
    DSWAP: 'interchanges two vectors.',

    ICAMAX: 'finds the index of the first',
    IZAMAX: 'finds the index of the first element having maximum |Re(.)| + |Im(.)|',
    ISAMAX: 'finds the index of the first element having maximum absolute value.',
    IDAMAX: 'finds the index of the first element having maximum absolute value.',

    SCASUM: 'takes the sum of the (|Re(.)| + |Im(.)|)\'s of a complex vector and returns a single precision result.',
    DZASUM: 'takes the sum of the (|Re(.)| + |Im(.)|)\'s of a complex vector and returns a single precision result.',
    SASUM: 'takes the sum of the absolute values.',
    DASUM: 'takes the sum of the absolute values.',

    SCNRM2: 'returns the euclidean norm of a vector via the function name, so that SCNRM2 := sqrt( x**H*x )',
    SNRM2: 'returns the euclidean norm of a vector via the function name, so that SNRM2 := sqrt( x\'*x ).',
    DZNRM2: 'returns the euclidean norm of a vector via the function name, so that',
    DNRM2: 'returns the euclidean norm of a vector via the function name, so that DNRM2 := sqrt( x\'*x )',

    SROT: 'applies a plane rotation.',
    DROT: 'applies a plane rotation.',
    CSROT: 'applies a plane rotation, where the cos and sin (c and s) are real',
    ZDROT: 'Applies a plane rotation, where the cos and sin (c and s) are real and the vectors cx and cy are complex.',

    SROTG: 'construct givens plane rotation.',
    CROTG: 'determines a complex Givens rotation.',
    ZROTG: 'determines a double complex Givens rotation.',
    DROTG: 'construct givens plane rotation.',

    DROTM: 'APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX.',
    SROTM: 'APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX',

    DROTMG: 'CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H',
    SROTMG: 'CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H'
};

(function i() {

    const props = Object.getOwnPropertyNames(desc);
    const reverted = props.map(name => ({ o: name, n: name.split('').reverse().join('') }));
    const reverseSorted = reverted.sort((o1, o2) => o1.n > o2.n ? 1 : o1.n < o2.n ? -1 : 0);
    reverseSorted.forEach(o => console.log(`${o.o}`, desc[o.o]));


})();

/**
x SCNRM2 returns the euclidean norm of a vector via the function name, so that SCNRM2 := sqrt( x**H*x )
x DNRM2 returns the euclidean norm of a vector via the function name, so that DNRM2 := sqrt( x'*x )
x SNRM2 returns the euclidean norm of a vector via the function name, so that SNRM2 := sqrt( x'*x ).
x DZNRM2 returns the euclidean norm of a vector via the function name, so that

x CROTG determines a complex Givens rotation.
x DROTG construct givens plane rotation.
x SROTG construct givens plane rotation.
x ZROTG determines a double complex Givens rotation.

x DROTMG CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H
x SROTMG CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H

x DROTM APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX.
x SROTM APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX



x DROT applies a plane rotation.
x ZDROT Applies a plane rotation, where the cos and sin (c and s) are real and the vectors cx and cy are complex.
x SROT applies a plane rotation.
x CSROT applies a plane rotation, where the cos and sin (c and s) are real


x CSCAL scales a complex vector by a complex constant.
x DSCAL scales a vector by a constant.
x ZDSCAL scales a vector by a constant.
x SSCAL scales a vector by a constant.
x CSSCAL scales a complex vector by a real constant.
x ZSCAL scales a cplx vector by a cplx constant.


x SCASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and returns a single precision result.
x DASUM takes the sum of the absolute values.
x SASUM takes the sum of the absolute values.
x DZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and returns a single precision result.

x CSWAP interchanges two vectors.
x DSWAP interchanges two vectors.
x SSWAP interchanges two vectors.
x ZSWAP interchanges two vectors.

x CDOTU forms the dot product of two complex vectors
x ZDOTU forms the dot product of two complex vectors,  ZDOTU = X^T * Y
x CDOTC forms the dot product of two complex vectors
x ZDOTC forms the dot product of two complex vectors.

x DDOT forms the dot product of two vectors.
x SDOT forms the dot product of two vectors.
x DSDOT Compute the inner product of two vectors with extended precision accumulation and result.
x SDSDOT Compute the inner product of two vectors with extended precision accumulation.


ICAMAX finds the index of the first
IDAMAX finds the index of the first element having maximum absolute value.
ISAMAX finds the index of the first element having maximum absolute value.
IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|

CCOPY copies a vector x to a vector y.
DCOPY copies a vector, x, to a vector, y.
SCOPY copies a vector, x, to a vector, y.
ZCOPY copies a vector, x, to a vector, y.

CAXPY constant times a vector plus a vector.
DAXPY constant times a vector plus a vector.
SAXPY constant times a vector plus a vector.
ZAXPY constant times a vector plus a vector.
*/
