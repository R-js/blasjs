      PROGRAM test
c      solves system of equations
c         SUBROUTINE  ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
c     C := alpha*op( A )*op( B ) + beta*C,

      EXTERNAL ZGEMM, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC
      EXTERNAL COPYCC

      INTEGER N, M, K, LDA, LDB, LDC
      INTEGER NN, MM, KK

      PARAMETER (N=8)
      PARAMETER (K=4)
      PARAMETER (LDA=8)
      PARAMETER (LDB=8)
      PARAMETER (LDC=8)
                 
      CHARACTER TRANSA, TRANSB

      COMPLEX*16 ALPHA, BETA, AC(LDA, N),A(LDA, N)
      COMPLEX*16 BC(LDB, N), C(LDC,N), B(LDB,N)
                

        DATA ((AC(I,J),I=1,LDA),J=1,N)/
     +(1.2629542848807933,0.9921603654457979),
     +(-0.3262333607056494,-0.42951310949188126),
     +(1.3297992629225006,1.2383041008533804),
     +(1.2724293214294047,-0.2793462818542693),
     +(0.4146414344564082,1.7579030898107073),
     +(-1.5399500419037095,0.5607460908880562),
     +(-0.9285670347135381,-0.4527839725531578),
     +(-0.2947204467905602,-0.8320432961178319),
     +(-0.005767172747536955,-1.166570547084707),
     +(2.404653388857951,-1.0655905803882961),
     +(0.7635934611404596,-1.563782051071005),
     +(-0.7990092489893682,1.1565369971501793),
     +(-1.1476570092363514,0.8320471285723897),
     +(-0.28946157368822334,-0.22732869142475534),
     +(-0.29921511789731614,0.2661373616721048),
     +(-0.411510832795067,-0.3767027185836281),
     +(0.2522234481561323,2.4413646288945894),
     +(-0.8919211272845686,-0.7953391172553718),
     +(0.43568329935571865,-0.054877473711578625),
     +(-1.237538421929958,0.2501413228541527),
     +(-0.22426788527830935,0.6182432935662469),
     +(0.37739564598170106,-0.17262350264585732),
     +(0.1333363608148414,-2.2239002740099374),
     +(0.8041895097449078,-1.263614384970583),
     +(-0.057106774383808755,0.3587288959713519),
     +(0.5036079722337261,-0.011045478465663564),
     +(1.085769362145687,-0.9406491626186084),
     +(-0.6909538396968303,-0.11582532215695436),
     +(-1.2845993538721883,-0.8149687088699175),
     +(0.04672617218835198,0.24226348085968588),
     +(-0.23570655643950122,-1.4250983947324998),
     +(-0.5428882550102544,0.36594112304921983),
     +(-0.4333103174567822,0.2484126488725964),
     +(-0.6494716467962331,0.06528818167162072),
     +(0.726750747385451,0.01915639166027384),
     +(1.1519117540872,0.2573383771555333),
     +(0.9921603654457979,1.2629542848807933),
     +(-0.42951310949188126,-0.3262333607056494),
     +(1.2383041008533804,1.3297992629225006),
     +(-0.2793462818542693,1.2724293214294047),
     +(1.7579030898107073,0.4146414344564082),
     +(0.5607460908880562,-1.5399500419037095),
     +(-0.4527839725531578,-0.9285670347135381),
     +(-0.8320432961178319,-0.2947204467905602),
     +(-1.166570547084707,-0.005767172747536955),
     +(-1.0655905803882961,2.404653388857951),
     +(-1.563782051071005,0.7635934611404596),
     +(1.1565369971501793,-0.7990092489893682),
     +(0.8320471285723897,-1.1476570092363514),
     +(-0.22732869142475534,-0.28946157368822334),
     +(0.2661373616721048,-0.29921511789731614),
     +(-0.3767027185836281,-0.411510832795067),
     +(2.4413646288945894,0.2522234481561323),
     +(-0.7953391172553718,-0.8919211272845686),
     +(-0.054877473711578625,0.43568329935571865),
     +(0.2501413228541527,-1.237538421929958),
     +(0.6182432935662469,-0.22426788527830935),
     +(-0.17262350264585732,0.37739564598170106),
     +(-2.2239002740099374,0.1333363608148414),
     +(-1.263614384970583,0.8041895097449078),
     +(0.3587288959713519,-0.057106774383808755),
     +(-0.011045478465663564,0.5036079722337261),
     +(-0.9406491626186084,1.085769362145687),
     +(-0.11582532215695436,-0.6909538396968303)/


          DATA ((BC(I,J),I=1,N),J=1,N)/
     +(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),
     +(0,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),
     +(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),
     +(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),(0,0),
     +(0,0),(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),
     +(0,0),(0,0),(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),
     +(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0),(0,0),
     +(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0)/
     
     
c    
c      C := alpha*op( A )*op( B ) + beta*C,
c       
      PRINT * , "==CASE 0======="
       
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

        PRINT * , "==CASE 1======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)


         PRINT * , "==CASE 2======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(1,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

          PRINT * , "==CASE 3======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)
        return
         PRINT * , "==CASE 4======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 5======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(1,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

           PRINT * , "==CASE 6======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)
     
           PRINT * , "==CASE 7======="
         CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)
            
         PRINT * , "==CASE 8======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 9======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

           PRINT * , "==CASE 10======="
       
        
       CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(1,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 11======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 12======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='N'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

          PRINT * , "==CASE 13======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 14======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='N'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

          PRINT * , "==CASE 15======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

          PRINT * , "==CASE 16======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='T'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 17======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)


         PRINT * , "==CASE 18======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='C'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 19======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0.2,0.8)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

        PRINT * , "==CASE 20======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0.3,-0.7)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 21======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0,0)
        BETA=(1,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 22======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0,0)
        BETA=(0.2,0.2)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)

         PRINT * , "==CASE 23======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0,0)
        BETA=(0,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)
             
             PRINT * , "==CASE 24======="
       
        
        CALL CRCMPLX(AC,A, LDA,N)
        CALL CRCMPLX(AC,B, LDA,N)
        CALL CRCMPLX(AC,C, LDA,N)

             
        TRANSA='C'
        TRANSB='T'
        ALPHA=(0,0)
        BETA=(0.2,0)

        NN=8
        MM=6
        KK=4     

        CALL ZGEMM(TRANSA,TRANSB,MM,NN,KK,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        CALL PRNMATRC(C, 8, 8, 8)
    
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

      SUBROUTINE FILLM(A,V,M,N, LD)

       INTEGER N,M,K
       COMPLEX*16 A(LD,*),V
       
       K=1
       DO 20 J=1,N
         DO I=1,M
          A(I,J)=V+K
          K = K + 1
40       END DO            
20     CONTINUE

      END SUBROUTINE

      SUBROUTINE PRNMATR(A, N, M, LDA)

      INTEGER N,M
      DOUBLE PRECISION A(LDA,*)
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

      SUBROUTINE PRNMATRC(A, N, M, LDA)

      INTEGER N,M
      COMPLEX *16 A(LDA,*)
      PRINT *, '['
      DO J=1,N
         DO I=1,M
            PRINT *, "complex", A(I,J), ","
         END DO
      END DO
     
      PRINT *, ']'

      END SUBROUTINE
