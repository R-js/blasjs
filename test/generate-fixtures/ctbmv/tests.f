      PROGRAM test
c       SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
c        x := A*x,   or   x := A**T*x,   or   x := A**H*x,
      EXTERNAL ZTBMV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC
      EXTERNAL COPYCC


      INTEGER INCX,N,K,LDA
     
           
      CHARACTER UPLO, DIAG, TRANS

      COMPLEX*16 AC(6, 6), A(6, 6)
     
           
      COMPLEX*16 V6(6), X(6) 

       DATA ((AC(I,J),I=1,6),J=1,6)/
     + (1.2629542848807933,0.9921603654457979),
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

     
       DATA (V6(I),I=1,6)/
     + (-0.6490100777088978, 0.7721421858045301),
     +  ( -0.11916876241803812, -0.21951562675343952),
     +  ( 0.6641356998941105,-0.4248102833772871),
     +  ( 1.100969102194087,-0.418980099421959),
     +  ( 0.14377148075806995,  0.9969868609091059),
     +  ( -0.11775359816595128 ,-0.27577802908802723)/
     
        N=6
        LDA=6
        K=3
      PRINT * , "==CASE 0======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='N'
        UPLO='U'
        DIAG='N'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 1======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='N'
        UPLO='U'
        DIAG='U'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

            PRINT * , "==CASE 2======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='N'
        UPLO='L'
        DIAG='N'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

        PRINT * , "==CASE 3======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='N'
        UPLO='L'
        DIAG='U'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 4======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='T'
        UPLO='U'
        DIAG='N'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 5======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='T'
        UPLO='U'
        DIAG='U'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 6======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='T'
        UPLO='L'
        DIAG='N'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)
       
        PRINT * , "==CASE 7======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='T'
        UPLO='L'
        DIAG='U'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

          PRINT * , "==CASE 8== >> INCX =-1 << ====="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='C'
        UPLO='U'
        DIAG='N'
c
c       NOTE INCX CHANGED TO MINUS 1
c            
        INCX=-1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 9======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='C'
        UPLO='U'
        DIAG='U'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==CASE 10======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='C'
        UPLO='L'
        DIAG='N'
        INCX=1

        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

        PRINT * , "==CASE 11======="
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='C'
        UPLO='L'
        DIAG='U'
        INCX=1
        
        CALL ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

         PRINT * , "==>CASE 12 TRIVIAL N=0  "
       
        CALL CCPLX(V6,X,6)
        CALL CRCMPLX(AC, A, 6, 6)
        X(4) =  (0,0)
        X(2) =  (0,0)
       
        TRANS='C'
        UPLO='L'
        DIAG='U'
        INCX=1
        
        CALL ZTBMV(UPLO,TRANS,DIAG,0,K,A,LDA,X,INCX)
        CALL PRNVECC(X,6)

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
