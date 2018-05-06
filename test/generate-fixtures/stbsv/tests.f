      PROGRAM dsymvtest
c      
c     x := A*x,   or   x := A**T*x,
c SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 

      EXTERNAL DTBSV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC



      CHARACTER UPLO, TRANS, DIAG
      
      INTEGER N, K, LDA, INCX


      DOUBLE PRECISION X(6), XC(6) 
      
      DOUBLE PRECISION A(6,6) , AC(6,6), UPA(6,6), LPA(6,6)
            
      DATA ((LPA(I,J),I=1,6),J=1,6)/
     + 1.2629542848807933,
     +  -0.3262333607056494,
     +  1.3297992629225006,
     +  1.2724293214294047,
     +  0.4146414344564082,
     +  -1.5399500419037095,
     +  -0.2947204467905602,
     +  -0.005767172747536955,
     +  2.404653388857951,
     +  0.7635934611404596,
     +  -0.7990092489893682,
     +  0,
     +  -0.29921511789731614,
     +  -0.411510832795067,
     +  0.2522234481561323,
     +  -0.8919211272845686,
     +  0,
     +  0,
     +  0.37739564598170106,
     +  0.1333363608148414,
     +  0.8041895097449078,
     +  0,
     +  0,
     +  0,
     +  -1.2845993538721883,
     +  0.04672617218835198,
     +  0,
     +  0,
     +  0,
     +  0,
     +  1.1519117540872,
     +  0,
     +  0,
     +  0,
     +  0,
     +  0/      

      DATA ((UPA(I,J),I=1,6),J=1,6)/
     +  0,
     +   0,
     +   0,
     +   0,
     + 0,
     +  1.053750863028617,
     + 0,
     + 0,
     + 0,
     + 0,
     + 0.197684262345795,
     + -1.068692711254786,
     + 0,
     + 0,
     + 0,
     + 0.2626454586627623,
     + -1.2329011995712644,
     + -0.00372353379218051,
     + 0,
     + 0,
     + -0.9740025611125269,
     + 0.6893726977654734,
     + -0.955839103276798,
     + -1.2317070584140966,
     + 0,
     + -0.9106806824932887,
     + 0.7412763052602079,
     + 0.06851153327714439,
     + -0.3237507545879617,
     + -1.0865030469936974,
     + -0.767790184730859,
     + -1.1197200611269833,
     + -0.4481742366033955,
     + 0.47173637445323024,
     + -1.180490682884277,
     + 1.4702569970829857/

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

      PRINT * , "==CASE 0======="
     
      uplo = 'U'
      trans='N'
      diag='N'
      N=6
      K=5
      lda=6
      incx = -1
     
            
      CALL COPY(XC,X,6)
      CALL COPYMA(UPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 1======="
     
      uplo = 'U'
      trans='N'
      diag='U'
      N=6
      K=5
      lda=6
      incx = 1
     
            
      CALL COPY(XC,X,6)
      CALL COPYMA(UPA,A,6,6)  
      X(3)=0
      X(6)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 2======="
     
      uplo = 'U'
      trans='T'
      diag='N'
      N=6
      K=5
      lda=6
      incx = 1
     
            
      CALL COPY(XC,X,6)
      CALL COPYMA(UPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 3======="
     
      uplo = 'U'
      trans='T'
      diag='U'
      N=6
      K=5
      lda=6
      incx = 1
     
            
      CALL COPY(XC,X,6)
      CALL COPYMA(UPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 4======="
     
      uplo = 'L'
      trans='N'
      diag='N'
      N=6
      K=5
      lda=6
      incx = 1
     
      CALL COPY(XC,X,6)
      CALL COPYMA(LPA,A,6,6)  
      X(3)=0
      X(1)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 5======="
     
      uplo = 'L'
      trans='N'
      diag='U'
      N=6
      K=5
      lda=6
      incx = 1
     
      CALL COPY(XC,X,6)
      CALL COPYMA(LPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

       PRINT * , "==CASE 6======="
     
      uplo = 'L'
      trans='T'
      diag='N'
      N=6
      K=5
      lda=6
      incx = 1
     
      CALL COPY(XC,X,6)
      CALL COPYMA(LPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c 
      PRINT *,"B="
      CALL PRNVEC(X,N)

        PRINT * , "==CASE 7======="
     
      uplo = 'L'
      trans='T'
      diag='U'
      N=6
      K=5
      lda=6
      incx = 1
     
      CALL COPY(XC,X,6)
      CALL COPYMA(LPA,A,6,6)  
      X(3)=0
      CALL DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
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

