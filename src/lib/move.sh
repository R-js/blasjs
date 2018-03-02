mkdir -p l1/single
mkdir -p l1/double
mkdir -p l1/complex
mkdir -p l1/dcomplex
mkdir -p l2/single
mkdir -p l2/double
mkdir -p l2/complex
mkdir -p l2/dcomplex
mkdir -p l3/single
mkdir -p l3/double
mkdir -p l3/complex
mkdir -p l3/dcomplex

# level 1

# single l1

mv srotg.f l1/single
mv srotmg.f l1/single
mv srot.f l1/single
mv srotm.f l1/single
mv sswap.f l1/single
mv sscal.f l1/single
mv scopy.f l1/single
mv saxpy.f l1/single
mv sdot.f l1/single
mv sdsdot.f l1/single
mv snrm2.f l1/single
mv scnrm2.f l1/single
mv sasum.f l1/single
mv isamax.f l1/single

# l1 double

mv drotg.f l1/double
mv drotmg.f l1/double
mv drot.f l1/double
mv drotm.f l1/double
mv dswap.f l1/double
mv dscal.f l1/double
mv dcopy.f l1/double
mv daxpy.f l1/double
mv ddot.f l1/double
mv dsdot.f l1/double
mv dnrm2.f l1/double
mv dznrm2.f l1/double
mv dasum.f l1/double
mv idamax.f l1/double

# l1 complex

mv crotg.f l1/complex
mv csrot.f l1/complex
mv cswap.f l1/complex
mv cscal.f l1/complex
mv csscal.f l1/complex
mv ccopy.f l1/complex
mv caxpy.f l1/complex
mv cdotu.f l1/complex
mv cdotc.f l1/complex
mv scasum.f l1/complex
mv icamax.f l1/complex

# l1 double comlpex

mv zrotg.f l1/dcomplex
mv zdrotf.f l1/dcomplex
mv zdrot.f l1/dcomplex
mv zswap.f l1/dcomplex
mv zscal.f l1/dcomplex
mv zdscal.f l1/dcomplex
mv zcopy.f l1/dcomplex
mv zaxpy.f l1/dcomplex
mv zdotu.f l1/dcomplex
mv zdotc.f l1/dcomplex
mv dzasum.f l1/dcomplex
mv izamax.f l1/dcomplex

# level 2

mv sgemv.f l2/single
mv sgbmv.f l2/single
mv ssymv.f l2/single
mv ssbmv.f l2/single
mv sspmv.f l2/single
mv strmv.f l2/single
mv stbmv.f l2/single
mv stpmv.f l2/single
mv strsv.f l2/single
mv stbsv.f l2/single
mv stpsv.f l2/single
mv sger.f l2/single
mv ssyr.f l2/single
mv sspr.f l2/single
mv ssyr2.f l2/single
mv sspr2.f l2/single

mv dgemv.f l2/double
mv dgbmv.f l2/double
mv dsymv.f l2/double
mv dsbmv.f l2/double
mv dspmv.f l2/double
mv dtrmv.f l2/double
mv dtbmv.f l2/double
mv dtpmv.f l2/double
mv dtrsv.f l2/double
mv dtbsv.f l2/double
mv dtpsv.f l2/double
mv dger.f l2/double
mv dsyr.f l2/double
mv dspr.f l2/double
mv dsyr2.f l2/double
mv dspr2.f l2/double

mv cgemv.f l2/complex
mv cgbmv.f  l2/complex
mv chemv.f  l2/complex
mv chbmv.f l2/complex
mv chpmv.f  l2/complex
mv ctrmv.f  l2/complex
mv ctbmv.f  l2/complex
mv ctpmv.f  l2/complex
mv ctrsv.f  l2/complex
mv ctbsv.f  l2/complex
mv ctpsv.f  l2/complex
mv cgeru.f  l2/complex
mv cgerc.f  l2/complex
mv cher.f l2/complex
mv chpr.f l2/complex
mv cher2.f l2/complex
mv chpr2.f l2/complex

# double complex

mv zgemv.f l2/dcomplex
mv zgbmv.f l2/dcomplex
mv zhemv.f l2/dcomplex
mv zhbmv.f l2/dcomplex
mv zhpmv.f l2/dcomplex
mv ztrmv.f l2/dcomplex
mv ztbmv.f l2/dcomplex
mv ztpmv.f l2/dcomplex
mv ztrsv.f l2/dcomplex
mv ztbsv.f l2/dcomplex
mv ztpsv.f l2/dcomplex
mv zgeru.f  l2/dcomplex
mv zgerc.f  l2/dcomplex
mv zher.f  l2/dcomplex
mv zhpr.f  l2/dcomplex
mv zher2.f  l2/dcomplex
mv zhpr2.f  l2/dcomplex

# level 3

mv sgemm.f l3/single
mv ssymm.f l3/single
mv ssyrk.f l3/single
mv ssyr2k.f l3/single
mv strmm.f l3/single
mv strsm.f l3/single

mv dgemm.f l3/double
mv dsymm.f l3/double
mv dsyrk.f l3/double
mv dsyr2k.f l3/double
mv dtrmm.f l3/double
mv dtrsm.f l3/double

mv cgemm.f l3/complex
mv csymm.f l3/complex
mv chemm.f l3/complex
mv csyrk.f l3/complex
mv cherk.f l3/complex
mv csyr2k.f l3/complex
mv cher2k.f l3/complex
mv ctrmm.f l3/complex
mv ctrsm.f l3/complex

mv zgemm.f l3/dcomplex
mv zsymm.f l3/dcomplex
mv zhemm.f l3/dcomplex
mv zsyrk.f l3/dcomplex
mv zherk.f l3/dcomplex
mv zsyr2k.f l3/dcomplex
mv zher2k.f l3/dcomplex
mv ztrmm.f l3/dcomplex
mv ztrsm.f l3/dcomplex

