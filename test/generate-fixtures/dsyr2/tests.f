      PROGRAM dsymvtest
c      
c   A := alpha*x*x**T + A
c
c     HELPERS
c
c

c SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)

      EXTERNAL DSYR2, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC



      CHARACTER UPLO
      
      INTEGER N, INCX


      DOUBLE PRECISION ALPHA, BETA, DUMMY

      DOUBLE PRECISION X(6), XC(6), Y(6), YC(6)  
      
      DOUBLE PRECISION A(6,6) , AC(6,6)
  
          


      DATA   ((AC(I,J), I=1,6),J=1,6)/1.053750863028617, 0,0,0, 0, 0,
     +  0.197684262345795,-1.068692711254786,0, 0,0, 0,
     +   0.2626454586627623,-1.2329011995712644,
     + -0.00372353379218051,
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

    

     
      DATA     (XC(I),I=1,6)/ 
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523/

       DATA     (YC(I),I=1,6)/ 
     + -7.7995958197861910,
     + -0.90987740363925695,
     + -1.3999755270779133,
     + 14.987896066159010,
     + 6.2202929537743330,
     + 11.557443037629128/

      PRINT * , "==CASE 0======="
      uplo = 'U'
      N=6
      alpha=1
      lda=6
      incx = 1
      incy = 1

      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
c     DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)    
      CALL DSYR2(UPLO,N,ALPHA,X,INCX,Y, INCY, A,LDA)
      PRINT *,"A="
      CALL PRNMATR(A,LDA,N)

      PRINT * , "==CASE 1======="
      uplo = 'L'
      N=6
      alpha=1
      lda=6
      incx = -1
      incy = -1

      
      CALL COPY(XC,X,6)
      CALL COPY(YC,Y,6)
      CALL COPYMA(AC,A,6,6)  
c     DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)    
      CALL DSYR2(UPLO,N,ALPHA,X,INCX,Y, INCY, A,LDA)
      PRINT *,"A="
      CALL PRNMATR(A,LDA,N)
      

      PRINT * , "==CASE 2======="
     
      uplo = 'L'
      N=6
      alpha=1
      lda=6
      incx = -1
      incy = -1

      
      CALL COPY(XC,X,6)
      x(2)=0
      x(5)=0
      CALL COPY(YC,Y,6)
      y(1)=0
      y(2)=0
      y(4)=0
      CALL COPYMA(AC,A,6,6)  
c     DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)    
      CALL DSYR2(UPLO,N,ALPHA,X,INCX,Y, INCY, A,LDA)
      PRINT *,"A="
      CALL PRNMATR(A,LDA,N)

        PRINT * , "==CASE 3======="
     
      uplo = 'u'
      N=6
      alpha=1
      lda=6
      incx = -1
      incy = -1


      
      CALL COPY(XC,X,6)
      x(2)=0
      x(5)=0
      CALL COPY(YC,Y,6)
      y(1)=0
      y(2)=0
      y(4)=0
      CALL COPYMA(AC,A,6,6)  
c     DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)    
      CALL DSYR2(UPLO,N,ALPHA,X,INCX,Y, INCY, A,LDA)
      PRINT *,"A="
      CALL PRNMATR(A,LDA,N)


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

      SUBROUTINE FILL(X,N, V)

       INTEGER N
       DOUBLE PRECISION X(*),V

       DO 20 I=1,N
            X(I)=V
20     CONTINUE

      END SUBROUTINE

      SUBROUTINE FILLM(MM,N,M, V)

       INTEGER N,M
       DOUBLE PRECISION MM(M,*),V

       DO 20 J=1,N
         DO I=1,M
            MM(I,J)=V
40       END DO            
20     CONTINUE

      END SUBROUTINE

      SUBROUTINE PRNMATR(A, N, M)

      INTEGER N,M
      DOUBLE PRECISION A(M,*)
      PRINT *, '['
      DO J=1,N
         DO I=1,M
            PRINT *, A(I,J), ","
         END DO
      END DO
     
      PRINT *, ']'

      END SUBROUTINE

      SUBROUTINE PRNVEC(X, N)

      INTEGER N
      DOUBLE PRECISION X(*)

      PRINT *, '['
      DO J=1,N
            PRINT *, X(J), ", "
      END DO
     
      PRINT *, ']'

      END SUBROUTINE

