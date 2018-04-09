      PROGRAM dsymvtest
c      
c    A*x = b,   or   A**T*x = b,
c SUBROUTINE ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
c 
*> ZGBMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
*>
*>    y := alpha*A**H*x + beta*y,

      EXTERNAL ZGBMV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX



      CHARACTER TRANS
      COMPLEX*16 ALPHA,BETA
      
      INTEGER INCX,INCY,KL,KU,LDA,M,N
      PARAMETER (KL=3)
      PARAMETER (KU=2)
      PARAMETER (LDA=6)
      PARAMETER (N=8)
      PARAMETER (M=6)
      COMPLEX*16  A(LDA,N), XC(N), X(N), YC(M), Y(M)

      DOUBLE PRECISION ARE(LDA,N), AIM(LDA,N)
         
         
      DATA ((ARE(I,J),I=1,LDA),J=1,N)/
     +  1.2629542848807933,
     +   -0.3262333607056494,
     +   1.3297992629225006,
     +   1.2724293214294047,
     +   0,
     +   0,
     +   0.4146414344564082,
     +   -1.5399500419037095,
     +   -0.9285670347135381,
     +   -0.2947204467905602,
     +   -0.005767172747536955,
     +   0,
     +   2.404653388857951,
     +   0.7635934611404596,
     +   -0.7990092489893682,
     +   -1.1476570092363514,
     +   -0.28946157368822334,
     +   -0.29921511789731614,
     +   0,
     +   -0.411510832795067,
     +   0.2522234481561323,
     +   -0.8919211272845686,
     +   0.43568329935571865,
     +   -1.237538421929958,
     +   0,
     +   0,
     +   -0.22426788527830935,
     +   0.37739564598170106,
     +   0.1333363608148414,
     +   0.8041895097449078,
     +   0,
     +   0,
     +   0,
     +   -0.057106774383808755,
     +   0.5036079722337261,
     +   1.085769362145687,
     +   0,
     +   0,
     +   0,
     +   0,
     +   -0.6909538396968303,
     +   -1.2845993538721883,
     +   0,
     +   0,
     +   0,
     +   0,
     +   0,
     +   0.04672617218835198/   

      DATA ((AIM(I,J),I=1,LDA),J=1,N)/
     +   0.9921603654457979,
     +   -0.42951310949188126,
     +   1.2383041008533804,
     +   -0.2793462818542693,
     +   0,
     +   0,
     +   1.7579030898107073,
     +   0.5607460908880562,
     +   -0.4527839725531578,
     +   -0.8320432961178319,
     +   -1.166570547084707,
     +   0,
     +   -1.0655905803882961,
     +   -1.563782051071005,
     +   1.1565369971501793,
     +   0.8320471285723897,
     +   -0.22732869142475534,
     +   0.2661373616721048,
     +   0,
     +   -0.3767027185836281,
     +   2.4413646288945894,
     +   -0.7953391172553718,
     +   -0.054877473711578625,
     +   0.2501413228541527,
     +   0,
     +   0,
     +   0.6182432935662469,
     +   -0.17262350264585732,
     +   -2.2239002740099374,
     +   -1.263614384970583,
     +   0,
     +   0,
     +   0,
     +   0.3587288959713519,
     +   -0.011045478465663564,
     +   -0.9406491626186084,
     +   0,
     +   0,
     +   0,
     +   0,
     +   -0.11582532215695436,
     +   -0.8149687088699175,
     +   0,
     +   0,
     +   0,
     +   0,
     +   0,
     +   0.24226348085968588/


       DATA (XC(I),I=1,N)/
     +   (-0.6490100777088978,  0.7721421858045301),
     +   (-0.11916876241803812,-0.21951562675343952),
     +   (0.6641356998941105,-0.4248102833772871),
     +   (1.100969102194087,-0.418980099421959),
     +   (0.14377148075806995, 0.9969868609091059),
     +   (-0.11775359816595128,-0.27577802908802723),
     +   (-0.9120683669483379, 1.2560188173061),
     +   (-1.4375862408299789, 0.6466743904953449)/

       DATA (YC(I),I=1,M)/
     + (-0.6490100777088978, 0.7721421858045301),
     +  ( -0.11916876241803812, -0.21951562675343952),
     +  ( 0.6641356998941105,-0.4248102833772871),
     +  ( 1.100969102194087,-0.418980099421959),
     +  ( 0.14377148075806995,  0.9969868609091059),
     +  ( -0.11775359816595128 ,-0.27577802908802723)/
     
    

      PRINT * , "==CASE 0======="

        CALL CRCMPLX(ARE,AIM, A,N, LDA)
        CALL CCPLX(XC,X,N)
        CALL CCPLX(YC,Y,M)

        TRANS='N'
        beta= complex(2.5, +0.5)
        alpha= complex(0, 0)
       
        INCX=1
        INCY=1

        CALL ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y=',Y
     
      
    
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

