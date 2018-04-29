      PROGRAM test
    
c     x := A*x,   or   x := A**T*x,
c
c     HELPERS
c
c  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  

        EXTERNAL dgemm, COPY, ZEROA, copyma,FILL
        EXTERNAL FILLM, PRNMATR, PRNVEC

*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

      DOUBLE PRECISION ALPHA,BETA
      DOUBLE PRECISION A(6,6),B(6,6),C(6,6), M6X6(6,6)

      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB

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
c    DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
c     C := alpha*op( A )*op( B ) + beta*C,
c
      PRINT *, "=====CASE 0======="

      TRANSA = "N"
      TRANSB = "N"
      M = 4
      N = 6
      K = 3
      ALPHA = 0
      BETA=1

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      LDA=6
      LDB=6
      LDC=6

      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

      CALL PRNMATR(C, 6, 6) 

      PRINT *, "=====CASE 4======="

      TRANSA = "N"
      TRANSB = "C"
      M = 6
      N = 6
      K = 3
      ALPHA = 0
      BETA=0.75

      CALL COPYMA(M6X6, A, 6, 6)
      CALL COPYMA(M6X6, B, 6, 6)
      CALL COPYMA(M6X6, C, 6, 6)
      LDA=6
      LDB=6
      LDC=6

      CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      CALL PRNMATR(C, 6, 6)


                  PRINT *, "=====CASE 5======="

                  TRANSA = "N"
                  TRANSB = "N"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=-1.2

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

      PRINT *, "=====CASE 6======="

                  TRANSA = "N"
                  TRANSB = "N"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=0

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

            CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                  PRINT *, "=====CASE 7======="

                  TRANSA = "N"
                  TRANSB = "N"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=1

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)


                  PRINT *, "=====CASE 8======="

                  TRANSA = "C"
                  TRANSB = "N"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=1

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                      PRINT *, "=====CASE 9======="

                  TRANSA = "C"
                  TRANSB = "N"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=0

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                      PRINT *, "=====CASE 10======="

                  TRANSA = "N"
                  TRANSB = "C"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=0

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                     PRINT *, "=====CASE 11======="

                  TRANSA = "N"
                  TRANSB = "C"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA=1

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                   PRINT *, "=====CASE 12======="

                  TRANSA = "N"
                  TRANSB = "C"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA = 0.5

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                     PRINT *, "=====CASE 13======="

                  TRANSA = "C"
                  TRANSB = "C"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA = 0.5

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

                     PRINT *, "=====CASE 14======="

                  TRANSA = "C"
                  TRANSB = "C"
                  M = 4
                  N = 6
                  K = 3
                  ALPHA = 0.3
                  BETA = 0

                  CALL COPYMA(M6X6, A, 6, 6)
                  CALL COPYMA(M6X6, B, 6, 6)
                  CALL COPYMA(M6X6, C, 6, 6)
                  LDA=6
                  LDB=6
                  LDC=6

       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
                  CALL PRNMATR(C, 6, 6)

*> DGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,      
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.      

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

