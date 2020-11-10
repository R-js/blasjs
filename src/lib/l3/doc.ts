export const desc: { [index: string]: string } = {
    //
    CGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H',
    CTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B,   or   B := alpha*B*op( A )',
    CTRSM: 'solves one of the matrix equations op( A )*X = alpha*B,   or   X*op( A ) = alpha*B',
    CHEMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    CHER2K: 'performs one of the hermitian rank 2k operations C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C or C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C',
    CHERK: 'performs one of the hermitian rank k operations C := alpha*A*A**H + beta*C, or C := alpha*A**H*A + beta*C',
    CSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    CSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    CSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',

    ZGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X  or   op( X ) = X**T   or   op( X ) = X**H',
    ZHEMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or  C := alpha*B*A + beta*C',
    ZHER2K: 'performs one of the hermitian rank 2k operations C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C, or C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C',
    ZHERK: ' performs one of the hermitian rank k operations C := alpha*A*A**H + beta*C, or C := alpha*A**H*A + beta*C',
    ZSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    ZSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    ZSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    ZTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B,   or   B := alpha*B*op( A ) where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or non-unit,  upper or lower triangular matrix  and  op( A )  is one  of op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.',
    ZTRSM: 'solves one of the matrix equations op( A )*X = alpha*B,   or   X*op( A ) = alpha*B',

    DGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C',
    DSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    DSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    DSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    DTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B,   or   B := alpha*B*op( A )',
    DTRSM: 'solves one of the matrix equations op( A )*X = alpha*B,   or   X*op( A ) = alpha*B',

    SGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C',
    SSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    SSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    SSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    STRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B,   or   B := alpha*B*op( A )',
    STRSM: 'solves one of the matrix equations op( A )*X = alpha*B,   or   X*op( A ) = alpha*B'
};


export const x = {

    //CHER2K: 'performs one of the hermitian rank 2k operations C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C or C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C',
    //ZHER2K: 'performs one of the hermitian rank 2k operations C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C, or C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C',
    //
    //CSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    //DSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    //SSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    //ZSYR2K: 'performs one of the symmetric rank 2k operations C := alpha*A*B**T + alpha*B*A**T + beta*C, or C := alpha*A**T*B + alpha*B**T*A + beta*C',
    //
    //CHERK: 'performs one of the hermitian rank k operations C := alpha*A*A**H + beta*C, or C := alpha*A**H*A + beta*C',
    //ZHERK: ' performs one of the hermitian rank k operations C := alpha*A*A**H + beta*C, or C := alpha*A**H*A + beta*C',
    //
    //CSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    //DSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    //SSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    //ZSYRK: 'performs one of the symmetric rank k operations C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C',
    //
    //CGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C, where op( X ) is one of op( X ) = X or op( X ) = X**T or op( X ) = X**H',
    //DGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C',
    //SGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C',
    //ZGEMM: 'performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C, where op( X ) is one of op( X ) = X or op( X ) = X**T or op( X ) = X**H',
    //
    //CHEMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //ZHEMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //
    //CTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B, or B := alpha*B*op( A )',
    //DTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B, or B := alpha*B*op( A )',
    //STRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B, or B := alpha*B*op( A )',
    //ZTRMM: 'performs one of the matrix-matrix operations B := alpha*op( A )*B, or B := alpha*B*op( A ) where alpha is a scalar, B is an m by n matrix, A is a unit, or non-unit, upper or lower triangular matrix and op( A ) is one of op( A ) = A or op( A ) = A**T or op( A ) = A**H.',
    //
    //CSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //DSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //SSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //ZSYMM: 'performs one of the matrix-matrix operations C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C',
    //
    CTRSM: 'solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B',
    DTRSM: 'solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B',
    STRSM: 'solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B',
    ZTRSM: 'solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B',
    //
};

