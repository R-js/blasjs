      PROGRAM test

c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDC,N
*       CHARACTER TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),C(LDC,*)
*       ..

      EXTERNAL ZSYRK, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC, CCOPYMA

          
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDC,N
      CHARACTER TRANS,UPLO

      COMPLEX*16 M6x6(6,6)
      COMPLEX*16 A(6,6) ,B(6,6), C(6,6)
                 
      
        DATA ((M6X6(I,J),I=1,6), J=1,6)/
     + (1.2629542848807933, 0.9921603654457979),
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


        PRINT *, "==CASE 1=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        TRANS='N'
        UPLO='U'
      
        ALPHA = (0, 0)
        BETA =  (0.2, -0.3)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

       PRINT *, "==CASE 2=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        TRANS='N'
        UPLO='L'
      
        ALPHA = (0, 0)
        BETA =  (0.2, -0.3)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

          PRINT *, "==CASE 3=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        TRANS='N'
        UPLO='L'
      
        ALPHA = (0, 0)
        BETA =  (0, 0)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

          PRINT *, "==CASE 4=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        A(1,1)=0

        TRANS='N'
        UPLO='U'
      
        ALPHA = (0.2, 0.8)
        BETA =  (0.2, -0.3)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

        PRINT *, "==CASE 5=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        A(1,1)=0

        TRANS='N'
        UPLO='L'
      
        ALPHA = (0.2, 0.8)
        BETA =  (0, 0)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

          PRINT *, "==CASE 6=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        A(1,1)=0

        TRANS='N'
        UPLO='L'
      
        ALPHA = (0.2, 0.8)
        BETA =  (1, 0)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

        PRINT *, "==CASE 7=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        A(1,1)=0

        TRANS='T'
        UPLO='U'
      
        ALPHA = (0.2, 0.8)
        BETA =  (-.3, .7)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

        PRINT *, "==CASE 8=="
c  SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL CCOPYMA(M6X6, A, 6,6)
        
        CALL CCOPYMA(M6X6, C, 6,6)

        A(1,1)=0

        TRANS='T'
        UPLO='L'
      
        ALPHA = (0.2, 0.8)
        BETA =  (0,0)
        
        N=4
        K=6

        LDA=6
        LDB=6
        LDC=6
c  ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)       
        CALL ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
        CALL PRNMATRC(C, 6, 6)

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

       SUBROUTINE CCOPYMA(AC,A,N,LDA)

        INTEGER N, LDA
        COMPLEX *16 AC(LDA, *), A(LDA,*)

       DO 400 J = 1, N
         DO 500 I = 1, LDA
              A(I,J)=AC(I,J);
500       CONTINUE
400     CONTINUE
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

      SUBROUTINE CRCMPLX(RE,IM,A,N,LDA)

        INTEGER N, LDA, CUR
        COMPLEX*16 A(LDA,N)

        DOUBLE PRECISION RE(LDA*N), IM(LDA*N)

       DO 40 J = 1, N
         DO 50 I = 1, LDA
              CUR = (J-1)*LDA+I
              A(I,J)=COMPLEX(RE(CUR),IM(CUR));
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
