      PROGRAM test
    
* SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

        EXTERNAL dsyr2k, COPY, ZEROA, copyma, FILL
        EXTERNAL FILLM, PRNMATR, PRNVEC

         DOUBLE PRECISION ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,N
        CHARACTER TRANS,UPLO

        DOUBLE PRECISION A(6,6),B(6,6),C(6,6), M6x6(6,6)

        DATA ((M6X6(I,J),I=1,6), J=1,6)/
     + 1.2629542848807933,
     + -0.3262333607056494,
     + 1.3297992629225006,
     + 1.2724293214294047,
     + 0.4146414344564082,
     + -1.5399500419037095,
     + -0.9285670347135381,
     + -0.2947204467905602,
     + -0.005767172747536955,
     + 2.404653388857951,
     + 0.7635934611404596,
     + -0.7990092489893682,
     + -1.1476570092363514,
     + -0.28946157368822334,
     + -0.29921511789731614,
     + -0.411510832795067,
     + 0.2522234481561323,
     + -0.8919211272845686,
     + 0.43568329935571865,
     + -1.237538421929958,
     + -0.22426788527830935,
     + 0.37739564598170106,
     + 0.1333363608148414,
     + 0.8041895097449078,
     + -0.057106774383808755,
     + 0.5036079722337261,
     + 1.085769362145687,
     + -0.6909538396968303,
     + -1.2845993538721883,
     + 0.04672617218835198,
     + -0.23570655643950122,
     + -0.5428882550102544,
     + -0.4333103174567822,
     + -0.6494716467962331,
     + 0.726750747385451,
     + 1.1519117540872/
c
*     DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 
c  
c
      PRINT *, "=====CASE 2======="
      
       ALPHA=0
       BETA=0.25
       K=3
       LDA=6
       LDB=6
       LDC=6
       N=6
       TRANS="N"
       UPLO="U"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)

      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6) 

      PRINT *, "=====CASE 3======="
      
       ALPHA=0
       BETA=0
       K=3
       LDA=6
       LDB=6
       LDC=6
       N=6
       TRANS="N"
       UPLO="U"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)

      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6) 

        PRINT *, "=====CASE 4======="
      
       ALPHA=0
       BETA=0
       K=3
       LDA=6
       LDB=6
       LDC=6
       N=6
       TRANS="N"
       UPLO="L"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)

      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6) 

       PRINT *, "=====CASE 5======="
      
       ALPHA=0.3
       BETA=-0.5
       K=3
       LDA=6
       LDB=6
       LDC=6
       N=6
       TRANS="N"
       UPLO="U"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0

      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)
    
         PRINT *, "=====CASE 6======="
      
       ALPHA=0.3
       BETA=0
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       TRANS="N"
       UPLO="U"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

        PRINT *, "=====CASE 7======="
      
       ALPHA=0.3
       BETA=1
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       TRANS="N"
       UPLO="U"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

       PRINT *, "=====CASE 8======="
      
       ALPHA=0.3
       BETA=0
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       TRANS="N"
       UPLO="L"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

        PRINT *, "=====CASE 9======="
      
       ALPHA=0.3
       BETA=0.5
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       TRANS="N"
       UPLO="L"

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

        PRINT *, "=====CASE 10======="
      
       ALPHA=0.3
       BETA=1
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="L"
       TRANS="N"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

       PRINT *, "=====CASE 11======="
      
       ALPHA=0.3
       BETA=1
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="U"
       TRANS="T"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

      PRINT *, "=====CASE 12======="
      
       ALPHA=0.3
       BETA=0
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="U"
       TRANS="T"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

        PRINT *, "=====CASE 13======="
      
       ALPHA=0.3
       BETA=0.4
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="U"
       TRANS="T"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

       PRINT *, "=====CASE 14======="
      
       ALPHA=0.3
       BETA=0
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="L"
       TRANS="T"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)

      PRINT *, "=====CASE 15======="
      
       ALPHA=0.3
       BETA=0.2
       N=6
       K=3
       LDA=6
       LDB=6
       LDC=6
      
       UPLO="L"
       TRANS="T"
     

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      A(1,1)=0;
      B(1,1)=0
      CALL DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) 

      CALL PRNMATR(C, 6, 6)
    
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

