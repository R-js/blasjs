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
