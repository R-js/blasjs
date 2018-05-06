      PROGRAM test
c      ZTRSV  solves one of the systems of equations
c      SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
c      A*x = b,   or   A**T*x = b,   or   A**H*x = b,

      EXTERNAL DTRSV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC
      EXTERNAL COPYCC

      INTEGER INCX,N, LDA

                 
      CHARACTER UPLO, DIAG, TRANS

      DOUBLE PRECISION ACLO(6, 6), ACUP(6, 6), A(6,6)
                
      DOUBLE PRECISION V6(6), X(6) 

        DATA ((ACUP(I,J),I=1,6),J=1,6)/
     +    1.053750863028617, 
     +    0,
     +    0,
     +    0,
     +    0,
     +    0,
     +    0.197684262345795,
     +    -1.068692711254786,
     +    0,
     +    0,
     +    0,
     +    0,
     +    0.2626454586627623,
     +    -1.2329011995712644,
     +    -0.00372353379218051,
     +    0,
     +    0,
     +    0,
     +    -0.9740025611125269,
     +    0.6893726977654734,
     +   -0.955839103276798,
     +   -1.2317070584140966,
     +    0,
     +    0,
     +    -0.9106806824932887,
     +    0.7412763052602079,
     +    0.06851153327714439,
     +    -0.3237507545879617,
     +    -1.0865030469936974,
     +    0,
     +   -0.767790184730859,
     +   -1.1197200611269833,
     +   -0.4481742366033955,
     +    0.47173637445323024,
     +   -1.180490682884277,
     +    1.4702569970829857/
            
            DATA ((ACLO(I,J),I=1,6),J=1,6)/
     +   1.2629542848807933,
     +   -0.3262333607056494,
     +   1.3297992629225006,
     +   1.2724293214294047,
     +   0.4146414344564082,
     +   -1.5399500419037095,
     +   0,
     +   -0.2947204467905602,
     +   -0.005767172747536955,
     +   2.404653388857951,
     +   0.7635934611404596,
     +   -0.7990092489893682,
     +   0,
     +   0,
     +   -0.29921511789731614,
     +   -0.411510832795067,
     +   0.2522234481561323,
     +   -0.8919211272845686,
     +   0,
     +   0,
     +   0,
     +   0.37739564598170106,
     +   0.1333363608148414,
     +   0.8041895097449078,
     +   0,
     +   0,
     +   0,
     +   0,
     +   -1.2845993538721883,
     +   0.04672617218835198,
     +   0,
     +   0,
     +   0,
     +   0,
     +   0,
     +   1.1519117540872/

     
       DATA (V6(I),I=1,6)/
     + -0.08252376201716412,
     +  0.6060734308621007,
     + -0.8874201453170976,
     +  0.10542139019376515,
     +  0.3528744733184766,
     +  0.5503933584550523/
     
c        b := A*x,   or   b := A**T*x,   or   b := A**H*x,
      
      PRINT * , "======CASE 4======="
       
        LDA = 6
        N = 6
        INCX = 1
        
        UPLO = "L"
        TRANS = "N"
        DIAG = "N"

        CALL COPY(V6,X,6)
        X(1)=0
        CALL COPYMA(ACLO, A, 6,6)

        CALL DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      
        CALL PRNVEC(X,N)
    
      end


      SUBROUTINE COPY(X,Y,N)

       INTEGER N
       DOUBLE PRECISION X(*), Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE
      END SUBROUTINE

       SUBROUTINE CCPLX(X,Y,N)

       INTEGER N
       COMPLEX*16 X(*), Y(*)

       DO 22 I=1,N
            Y(I)=X(I)
22     CONTINUE
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


      SUBROUTINE COPYCC(AC,A,N,LDA)

        INTEGER N, LDA
        COMPLEX*16 AC(LDA, *), A(LDA,*)

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

      SUBROUTINE PRNVECC(X, N)

      INTEGER N
      COMPLEX*16 X(*)

      PRINT *, '['
      DO J=1,N
            PRINT *, "complex",X(J), ", "
      END DO
     
      PRINT *, ']'

      END SUBROUTINE

      SUBROUTINE CRCMPLX(AC,A,LDA,N)

        INTEGER N, LDA
        COMPLEX*16 AC(LDA,N), A(LDA,N)

       DO 40 J = 1, N
         DO 50 I = 1, LDA
              A(I,J)=AC(I,J);
50       CONTINUE
40     CONTINUE
      END SUBROUTINE

      SUBROUTINE PRNMATRC(A, N, M)

      INTEGER N,M
      COMPLEX *16 A(M,*)
      PRINT *, '['
      DO J=1,N
         DO I=1,M
            PRINT *, "complex", A(I,J), ","
         END DO
      END DO
     
      PRINT *, ']'

      END SUBROUTINE
