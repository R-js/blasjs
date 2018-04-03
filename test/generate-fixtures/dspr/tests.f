      PROGRAM dsprtest
c      
c     A := alpha*x*x**T + A
c
c     HELPERS
c
c
      INTEGER APSIZE
      INTEGER k

      EXTERNAL DSPR, COPY, ZEROA


      CHARACTER UPLO
      
      INTEGER N, INCX


      DOUBLE PRECISION ALPHA, X(8), XC(8)
      PARAMETER (N=8)
      PARAMETER (K=7)
      PARAMETER (APSIZE=N*(K+1)-K*(K+1)/2)
      DOUBLE PRECISION A(APSIZE), AC(APSIZE)
          


      DATA   (AC(I), I=1,APSIZE)/2.024761390077735,
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


     
      DATA     (XC(I),I=1,8)/ 
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523,
     + -1.1343309685168443,
     +  1.4623515387464268/
     
      PRINT * , "==CASE 0======="
      uplo = 'U'
      alpha = 1
      incx = 1


    
      CALL COPY(XC,X,N)
      CALL COPY(AC,A,N);
  
      CALL  DSPR(UPLO,N,ALPHA,X,INCX,A)

      PRINT *, A

      PRINT * , "==CASE 1======="


      uplo = 'L'
      alpha = 1
      incx = -1


    
      
      CALL COPY(XC,X,N);
      CALL COPY(AC,A,APSIZE);
  
      
      CALL  DSPR(UPLO,N,ALPHA,X,INCX,A)

      PRINT *, A

      PRINT * , "==CASE 2======="


      uplo = 'L'
      alpha = 1
      incx = -1
           
        X(1)=1
        X(2)=0
        X(3)=2
        X(4)=0
        X(5)=3
        X(6)=0
        X(7)=4
        X(8)=0
    
c       x: [1, 0, 2, 0, 3, 0, 4, 0],

      CALL COPY(AC,A,APSIZE);
  
      
      CALL  DSPR(UPLO,N,ALPHA,X,INCX,A)

      PRINT *, A

       PRINT * , "==CASE 3======="


      uplo = 'U'
      alpha = 1
      incx = -1
           
        X(1)=1
        X(2)=0
        X(3)=2
        X(4)=0
        X(5)=3
        X(6)=0
        X(7)=4
        X(8)=0
c       x: [1, 0, 2, 0, 3, 0, 4, 0],

      CALL COPY(AC,A,APSIZE);
  
      
      CALL  DSPR(UPLO,N,ALPHA,X,INCX,A)

      PRINT *, A

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

