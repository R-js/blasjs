      PROGRAM dsymvtest
c      
c    A*x = b,   or   A**T*x = b,
c SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
c 

      EXTERNAL DTPSV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC

      CHARACTER UPLO, TRANS, DIAG
      
      INTEGER N, K, LDA, INCX, SIZE
      PARAMETER (K=5)
      PARAMETER (N=6)

      PARAMETER (SIZE=-(K + 1) * K / 2 + N * (K + 1))

      DOUBLE PRECISION X(6), XC(6) 
      
      DOUBLE PRECISION APUPC(SIZE) , APLOC(SIZE), AP(SIZE)
            
      DATA (APUPC(I),I=1,SIZE)/
     +   1.053750863028617,
     +   0.197684262345795,
     +   -1.068692711254786,
     +   0.2626454586627623,
     +   -1.2329011995712644,
     +   -0.00372353379218051,
     +   -0.9740025611125269,
     +   0.6893726977654734,
     +   -0.955839103276798,
     +   -1.2317070584140966,
     +   -0.9106806824932887,
     +   0.7412763052602079,
     +   0.06851153327714439,
     +   -0.3237507545879617,
     +   -1.0865030469936974,
     +   -0.767790184730859,
     +   -1.1197200611269833,
     +   -0.4481742366033955,
     +   0.47173637445323024,
     +   -1.180490682884277,
     +   1.4702569970829857/      

      DATA (APLOC(I),I=1,SIZE)/
     + 1.4702569970829857,
     +   -1.180490682884277,
     +   0.47173637445323024,
     +   -0.4481742366033955,
     +   -1.1197200611269833,
     +   -0.767790184730859,
     +   -1.0865030469936974,
     +   -0.3237507545879617,
     +   0.06851153327714439,
     +   0.7412763052602079,
     +   -0.9106806824932887,
     +   -1.2317070584140966,
     +   -0.955839103276798,
     +   0.6893726977654734,
     +   -0.9740025611125269,
     +   -0.00372353379218051,
     +   -1.2329011995712644,
     +   0.2626454586627623,
     +   -1.068692711254786,
     +   0.197684262345795,
     +   1.053750863028617/
     
      DATA     (XC(I),I=1,6)/ 
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523/

      PRINT * , "==CASE 0======="

      uplo = 'U'
      trans='N'
      diag='N'
      incx = 1
                 
      CALL COPY(XC,X,6)
      CALL COPY(APUPC,AP,SIZE)

      X(3)=0
      X(6)=0
c      
c     DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)      
c
      CALL DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
c           
      PRINT *,"B="
      CALL PRNVEC(X,N)

      PRINT * , "==CASE 1======="

      uplo = 'U'
      trans='N'
      diag='U'
      incx = -1

            
      CALL COPY(XC,X,6)
      CALL COPY(APUPC,AP,SIZE)  
      X(3)=0
c     DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)      
      CALL DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
c           
      PRINT *,"B="
      CALL PRNVEC(X,N)
            

      PRINT * , "==CASE 2======="
     
      uplo = 'U'
      trans='T'
      diag='N'
      incx = 1
     
            
      CALL COPY(XC,X,6)
      CALL COPY(APUPC,AP,SIZE)  
      X(3)=0
c     DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)      
      CALL DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)
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

