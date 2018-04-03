      PROGRAM dsbmvcase0
      
      EXTERNAL DSPMV, COPY, ZEROA

      CHARACTER UPLO
      INTEGER N,K, INCX,INCY, APSIZE 
      DOUBLE PRECISION ALPHA, BETA, Y(8), X(8), YC(8), XC(8)
      PARAMETER (N=8)
      PARAMETER (K=7)
      PARAMETER (APSIZE=N*(K+1)-K*(K+1)/2)
      DOUBLE PRECISION AP(APSIZE)
          

      DATA   (AP(I), I=1,APSIZE)/2.024761390077735,
     +               -1.0457176517867426,
     +               -0.8962112639605789,
     +               -0.060634778120552485,
     +               -0.5013783177950893,
     +               0.9260627253330198,
     +               -1.4577072101223822,
     +               0.09505622698864298,
     +               0.8476649635960255,
     +               -1.6243645297605886,
     +               1.5761581812187344,
     +               -1.475547635260523,
     +               -0.14460820731054116,
     +               -1.0750101908181724,
     +               0.40654273194492646,
     +               -0.14727079003896576,
     +               1.5415930688269568,
     +               -0.9818556688038707,
     +               0.4965781726616992,
     +               1.6969478807230851,
     +               -0.26073630856812263,
     +               0.5013218277237272,
     +               -1.0135396704947914,
     +               1.6147522354681292,
     +               0.005641984852495504,
     +               -2.9048990603455724,
     +               -1.1071648189687495,
     +               1.5475669326182715,
     +               -0.10150344763172572,
     +               0.042650249796697896,
     +               -1.5967180142971973,
     +               0.490967372597059,
     +               0.421603365384753,
     +               1.8739038985953016,
     +               1.0345143239443348,
     +               0.08181031035401386/

      DATA     (YC(I),I=1,8)/
     +     0.7021167106675735,
     +     2.5071111484833684,
     +    -1.890027143624024,
     +    -0.5898127901911715,
     +    -1.7145022968458246,
     +    -0.4209978978166964,
     +     0.310141376504687,
     +     1.7025705860144955/
     
      DATA     (XC(I),I=1,8)/ 
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523,
     + -1.1343309685168443,
     +  1.4623515387464268/
     
      uplo = 'U'
      alpha = 1
      beta =0.25
      incx = 1
      incy = 1


      PRINT *, "apsize",APSIZE, "n=",n, "k=",k

      PRINT *,'case0'

      CALL COPY(XC,X,N)
      CALL COPY(YC,Y,N)

      
     
      PRINT *, "beta=", beta, "alpha=", alpha, "incx=", incx
      PRINT *,  "incy=", incy, "x=",x
      call DSPMV(UPLO,N,ALPHA,AP,x,INCX,BETA,Y,INCY)


      print *,I,"|", Y

      CALL COPY(XC,X,N)
      CALL COPY(YC,Y,N)
      uplo = 'U'
      alpha = 1
      beta =0.25
      incx = -1
      incy = -1

      PRINT *,'case1'
      PRINT *, "beta=", beta, "alpha=", alpha, "incx=", incx
      PRINT *,  "incy=", incy, "x=",x
      call DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
c     UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY
      print *,I,"|", Y
     
      CALL COPY(XC,X,N)
      CALL COPY(YC,Y,N)
      uplo = 'U'
      alpha = 0
      beta =1
      incx = -1
      incy = -1

      PRINT *,'case2'
      PRINT *, "beta=", beta, "alpha=", alpha, "incx=", incx
      PRINT *,  "incy=", incy, "x=",z
      call DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
c     UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY
      print *,I,"|", Y

      CALL COPY(XC,X,N)
      CALL COPY(YC,Y,N)
      uplo = 'U'
      alpha = 0
      beta =0
      incx = -1
      incy = -1

      PRINT *,'case3'
      PRINT *, "beta=", beta, "alpha=", alpha, "incx=", incx
      PRINT *,  "incy=", incy, "x=",z
      call DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
c     UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY
      print *,I,"|", Y


c     CASE4

      CALL COPY(XC,X,N)
      CALL COPY(YC,Y,N)
c      CALL ZEROA(X,N);
c      X(1)=1.0

      uplo = 'L'
      alpha = 1
      beta =1
      incx = 1
      incy = 1

      PRINT *, "case4"
      PRINT *, "beta=", beta, "alpha=", alpha, "incx=", incx
      PRINT *, "incy=", incy, "x=",z
      call DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
c     UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY
      print *,I,"|", Y

      
     
      end

      SUBROUTINE COPY(X,Y,N)

       INTEGER N
       DOUBLE PRECISION X(*), Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE

      END SUBROUTINE

      SUBROUTINE ZEROA(X,N)

       INTEGER N
       DOUBLE PRECISION X(*)

       DO 20 I=1,N
            X(I)=0.0
20     CONTINUE

      END SUBROUTINE
     