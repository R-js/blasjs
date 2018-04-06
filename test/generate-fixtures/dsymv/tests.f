      PROGRAM dsymvtest
c      
c    y := alpha*A*x + beta*y,
c
c     HELPERS
c
c
c SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      EXTERNAL DSYMV, COPY, ZEROA, COPYMA


      CHARACTER UPLO
      
      INTEGER N, INCX, INCY


      DOUBLE PRECISION ALPHA, BETA
      DOUBLE PRECISION X(6), XC(6), Y(6), YC(6), A(6,6) , AC(6,6)
  
          


      DATA   ((AC(I,J), I=1,6),J=1,6)/1.053750863028617, 0,0,0, 0, 0,
     +  0.197684262345795,-1.068692711254786,0, 0,0, 0,
     +   0.2626454586627623,-1.2329011995712644,-0.00372353379218051,
     +   0,0,0,
     +   -0.9740025611125269,
     +   0.6893726977654734,
     +   -0.955839103276798,
     +   -1.2317070584140966,
     +   0,
     +   0,
     +   -0.9106806824932887,
     +   0.7412763052602079,
     +   0.06851153327714439,
     +   -0.3237507545879617,
     +   -1.0865030469936974,
     +   0,
     +   -0.767790184730859,
     +   -1.1197200611269833,
     +   -0.4481742366033955,
     +   0.47173637445323024,
     +   -1.180490682884277,
     +   1.4702569970829857/

      DATA     (YC(I),I=1,6)/
     +     0.7021167106675735,
     +     2.5071111484833684,
     +    -1.890027143624024,
     +    -0.5898127901911715,
     +    -1.7145022968458246,
     +    -0.4209978978166964/

     
      DATA     (XC(I),I=1,6)/ 
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523/
     
      PRINT * , "==CASE 0======="
      uplo = 'U'
      alpha = 1
      beta=1
      incx = 1
      incy = 1
      lda=6
      n=6
      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
      CALL  DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      PRINT *, Y

      PRINT * , "==CASE 1======="
      uplo = 'U'
      alpha = 0
      beta = 0.35
      incx = 1
      incy = 1
      lda = 6
      n = 6
      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
      CALL  DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      PRINT *, Y

       PRINT * , "==CASE 2======="
      uplo = 'U'
      alpha = 0.5
      beta = 0
      incx = -1
      incy = -1
      lda = 6
      n = 6
      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
      CALL  DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      PRINT *, Y

      PRINT * , "==CASE 3======="
      uplo = 'L'
      alpha = 0.5
      beta = 0
      incx = -1
      incy = -1
      lda = 6
      n = 6
      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
      CALL  DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      PRINT *, Y


      PRINT * , "==CASE 4======="
      uplo = 'L'
      alpha = 0
      beta = 1
      incx = -1
      incy = -1
      lda = 6
      n = 6
      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
      CALL  DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      PRINT *, Y


      end





      SUBROUTINE COPY(X,Y,N)

       INTEGER N
       DOUBLE PRECISION X(*), Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE
      END SUBROUTINE

      SUBROUTINE COPYMA(AC,A,N,LDA)

        INTEGER N, LDA
        DOUBLE PRECISION AC(LDA, *), A(LDA,*)

       DO 40 J = 1, N
         DO 50 I = 1, LDA
              A(I,J)=AC(I,J);
50       CONTINUE
40     CONTINUE
      END SUBROUTINE

      SUBROUTINE ZEROA(X,N)

       INTEGER N
       DOUBLE PRECISION X(*)

       DO 20 I=1,N
            X(I)=0.0
20     CONTINUE

      END SUBROUTINE

