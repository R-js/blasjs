      PROGRAM test

*> ZTRMM  performs one of the matrix-matrix operations
*>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*>    ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

      EXTERNAL ZTRMM, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC
      EXTERNAL COPYCC

      INTEGER  N, M, LDA

  
                 
      CHARACTER SIDE, UPLO, TRANSA, DIAG

      COMPLEX*16 M6X6(6, 6), A(6, 6), B(6,6)

      COMPLEX*16 ALPHA
                
c  
c     put Upper and lower matrix in the same 6x6 structure
c     at times we will snap 6x6 or 4x4 because the routine skips over elements it doenst use
c  
        DATA ((M6X6(I,J),I=1,6),J=1,6)/
     + (0,0),
     + (-0.3262333607056494,-0.42951310949188126),
     + (1.3297992629225006,1.2383041008533804),
     + (1.2724293214294047,-0.2793462818542693),
     + (0.4146414344564082,1.7579030898107073),
     + (-1.5399500419037095,0.5607460908880562),
     + (-0.9285670347135381,-0.4527839725531578),
     + (-0.2947204467905602,-0.8320432961178319),
     + (-0.005767172747536955,-1.166570547084707),
     + (2.404653388857951,-1.0655905803882961),
     + (0.7635934611404596,-1.563782051071005),
     + (-0.7990092489893682,1.1565369971501793),
     + (-1.1476570092363514,0.8320471285723897),
     + (-0.28946157368822334,-0.22732869142475534),
     + (-0.29921511789731614,0.2661373616721048),
     + (-0.411510832795067,-0.3767027185836281),
     + (0.2522234481561323,2.4413646288945894),
     + (-0.8919211272845686,-0.7953391172553718),
     + (0.43568329935571865,-0.054877473711578625),
     + (-1.237538421929958,0.2501413228541527),
     + (-0.22426788527830935,0.6182432935662469),
     + (0.37739564598170106,-0.17262350264585732),
     + (0.1333363608148414,-2.2239002740099374),
     + (0.8041895097449078,-1.263614384970583),
     + (-0.057106774383808755,0.3587288959713519),
     + (0.5036079722337261,-0.011045478465663564),
     + (1.085769362145687,-0.9406491626186084),
     + (-0.6909538396968303,-0.11582532215695436),
     + (-1.2845993538721883,-0.8149687088699175),
     + (0.04672617218835198,0.24226348085968588),
     + (-0.23570655643950122,-1.4250983947324998),
     + (-0.5428882550102544,0.36594112304921983),
     + (-0.4333103174567822,0.2484126488725964),
     + (-0.6494716467962331,0.06528818167162072),
     + (0.726750747385451,0.01915639166027384),
     + (1.1519117540872,0.2573383771555333)/
     
c     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
      PRINT * , "==CASE 0======="
       
c       FULL COPY
        CALL CRCMPLX(M6X6, A, 6,6)
        CALL CRCMPLX(M6X6, B, 6,6)
c        
c       B IS ALWAYS, MXN matrix, regardless if it is left or right multiplied
c            
c       IF SIDE = "L" A IS MXM matrix , 4x4
c       IF SIDE = "R" A iS NXN matrix , 6x6
        SIDE = 'L'
c       A IS UPPER OR LOWER TRIANGULAR MATRIX of (kxk)        
        UPLO = 'U'

c       [N]O (TRANSLATION or CONJUGATION of A, note!! A is k x k matrix so!!!)        
        TRANSA = 'N'
c       DIAGONALS ARE NOT UNIT  
        DIAG = 'N'
c       IF SIDE=L, A= (mxm) otherwise M,N,ALPHA,A,LDA,B,LDB
        M=4  
        N=6
        LDB=6
        LDA=6
        alpha = (0.2,0.8)
        
        CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, N, N)

         PRINT * , "==CASE 1======="
       
c       FULL COPY
        CALL CRCMPLX(M6X6, A, 6,6)
        CALL CRCMPLX(M6X6, B, 6,6)
c        
c       B IS ALWAYS, MXN matrix, regardless if it is left or right multiplied
c            
c       IF SIDE = "L" A IS MXM matrix , 4x4
c       IF SIDE = "R" A iS NXN matrix , 6x6
        SIDE = 'L'
c       A IS UPPER OR LOWER TRIANGULAR MATRIX of (kxk)        
        UPLO = 'U'

c       [N]O (TRANSLATION or CONJUGATION of A, note!! A is k x k matrix so!!!)        
        TRANSA = 'N'
c       DIAGONALS ARE NOT UNIT  
        DIAG = 'U'
c       IF SIDE=L, A= (mxm) otherwise M,N,ALPHA,A,LDA,B,LDB
        M=4  
        N=6
        LDB=6
        LDA=6
        alpha = (0.2,0.8)
        
        CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, N, N)

         PRINT * , "==CASE 2======="
       
c       FULL COPY
        CALL CRCMPLX(M6X6, A, 6,6)
        CALL CRCMPLX(M6X6, B, 6,6)
c        
c       B IS ALWAYS, MXN matrix, regardless if it is left or right multiplied
c            
c       IF SIDE = "L" A IS MXM matrix , 4x4
c       IF SIDE = "R" A iS NXN matrix , 6x6
        SIDE = 'L'
c       A IS UPPER OR LOWER TRIANGULAR MATRIX of (kxk)        
        UPLO = 'L'

c       [N]O (TRANSLATION or CONJUGATION of A, note!! A is k x k matrix so!!!)        
        TRANSA = 'N'
c       DIAGONALS ARE NOT UNIT  
        DIAG = 'N'
c       IF SIDE=L, A= (mxm) otherwise M,N,ALPHA,A,LDA,B,LDB
        M=4  
        N=6
        LDB=6
        LDA=6
        alpha = (0.2,0.8)

        CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, N, N)

         PRINT * , "==CASE 3======="
       
c       FULL COPY
        CALL CRCMPLX(M6X6, A, 6,6)
        CALL CRCMPLX(M6X6, B, 6,6)
c        
c       B IS ALWAYS, MXN matrix, regardless if it is left or right multiplied
c            
c       IF SIDE = "L" A IS MXM matrix , 4x4
c       IF SIDE = "R" A iS NXN matrix , 6x6
        SIDE = 'L'
c       A IS UPPER OR LOWER TRIANGULAR MATRIX of (kxk)        
        UPLO = 'L'

c       [N]O (TRANSLATION or CONJUGATION of A, note!! A is k x k matrix so!!!)        
        TRANSA = 'N'
c       DIAGONALS ARE NOT UNIT  
        DIAG = 'U'
c       IF SIDE=L, A= (mxm) otherwise M,N,ALPHA,A,LDA,B,LDB
        M=4  
        N=6
        LDB=6
        LDA=6
        alpha = (0.2,0.8)

        CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, N, N)
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
