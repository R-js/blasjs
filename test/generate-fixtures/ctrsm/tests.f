      PROGRAM test
c     solves system of equations 
c     ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
c     PRNMATRC, CCOPYMA
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,

      EXTERNAL ZTRSM, PRNMATRC, CCOPYMA
      
      INTEGER M,N,LDA,LDB
           
      CHARACTER SIDE, UPLO, TRANSA, DIAG 

      COMPLEX*16 ALPHA, A(6,6), B(6,6)
      COMPLEX*16 M6X6(6,6)
                

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
     

      PRINT *, "==CASE 2=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

        CALL CCOPYMA(M6X6, A, 6,6)
        CALL CCOPYMA(M6X6, B, 6,6)
  
        SIDE="L"
        TRANSA="N"
        UPLO="U"
        DIAG="N"
      
        ALPHA = (0,0)
        
        N=6
        M=4

        LDA=6
        LDB=6
       
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

         PRINT *, "==CASE 3=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

        CALL CCOPYMA(M6X6, A, 6,6)
        CALL CCOPYMA(M6X6, B, 6,6)
 
        SIDE="L"
        TRANSA="N"
        UPLO="U"
        DIAG="N"
      
        ALPHA = (0.2,0.6)
        
        N=6
        M=4

        LDA=6
        LDB=6
       
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

         PRINT *, "==CASE 4=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

        CALL CCOPYMA(M6X6, A, 6,6)
        CALL CCOPYMA(M6X6, B, 6,6)
        B(1,1)=0
        SIDE="L"
        TRANSA="N"
        UPLO="U"
        DIAG="U"
      
        ALPHA = (1,0)
        
        N=6
        M=4

        LDA=6
        LDB=6
       
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

  
                 PRINT *, "==CASE 5=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

        CALL CCOPYMA(M6X6, A, 6,6)
        CALL CCOPYMA(M6X6, B, 6,6)
        B(1,1)=0
        
        SIDE="L"
        UPLO="L"
        TRANSA="N"
        DIAG="N"
      
        ALPHA = (1,0)
        
        N=6
        M=4

        LDA=6
        LDB=6
        
         CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
         CALL PRNMATRC(B, 6, 6)

         PRINT *, "==CASE 6=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="L"
         TRANSA="N"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 7=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="U"
         TRANSA="T"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 8=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="U"
         TRANSA="C"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

           PRINT *, "==CASE 9=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="L"
         TRANSA="C"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

         PRINT *, "==CASE 10=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="L"
         TRANSA="C"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 11=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         B(1,1)=0
        
         SIDE="L"
         UPLO="L"
         TRANSA="T"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 12=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,2)=0
        
         SIDE="R"
         UPLO="U"
         TRANSA="N"
         DIAG="N"
      
         ALPHA = (1,0)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

          PRINT *, "==CASE 13=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,2)=0
        
         SIDE="R"
         UPLO="U"
         TRANSA="N"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)
         PRINT *, "==CASE 14=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(5,4)=0
        
         SIDE="R"
         UPLO="L"
         TRANSA="N"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

         PRINT *, "==CASE 15=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(5,4)=0
        
         SIDE="R"
         UPLO="L"
         TRANSA="N"
         DIAG="U"
      
         ALPHA = (1,0)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 16=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,6)=0
        
         SIDE="R"
         UPLO="U"
         TRANSA="T"
         DIAG="N"
      
         ALPHA = (1,0)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

        PRINT *, "==CASE 17=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,6)=0
        
         SIDE="R"
         UPLO="U"
         TRANSA="C"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)

          PRINT *, "==CASE 18=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,6)=0
        
         SIDE="R"
         UPLO="U"
         TRANSA="C"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)
          PRINT *, "==CASE 19=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(1,6)=0
        
         SIDE="R"
         UPLO="L"
         TRANSA="C"
         DIAG="U"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4  

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)
           PRINT *, "==CASE 20=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(2,1)=0
        
         SIDE="R"
         UPLO="L"
         TRANSA="C"
         DIAG="N"
      
         ALPHA = (0.5,0.5)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)
           PRINT *, "==CASE 21=="

*    SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)

         CALL CCOPYMA(M6X6, A, 6,6)
         CALL CCOPYMA(M6X6, B, 6,6)
         A(2,1)=0
        
         SIDE="R"
         UPLO="L"
         TRANSA="T"
         DIAG="N"
      
         ALPHA = (1,0)
        
         N=6
         M=4

        LDA=6
        LDB=6
        
        CALL ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        CALL PRNMATRC(B, 6, 6)
            
      end

         SUBROUTINE CCOPYMA(AC,A,N,LDA)

        INTEGER N, LDA
        COMPLEX *16 AC(LDA, *), A(LDA,*)

       DO 400 J = 1, N
         DO 500 I = 1, LDA
              A(I,J)=AC(I,J);
500       CONTINUE
400     CONTINUE
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

    

