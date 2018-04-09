      PROGRAM test


      EXTERNAL ZHBMV, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC


      INTEGER INCX,INCY,LDA,N, K
      PARAMETER (LDA=4)
      PARAMETER (N=6)
      PARAMETER (K=3)
    
      CHARACTER UPLO

      COMPLEX*16 ALPHA, AC(LDA,N), ACLO(LDA,N)
      COMPLEX*16 BETA
           
      COMPLEX*16  V6(6), Y(6), X(6) 


      DATA ((AC(I,J),I=1,LDA),J=1,N)/
     + (0,0),
     + (0,0),
     + (0,0),
     + (1.2629542350769043,-0.4295130968093872),
     + (0,0),
     + (0,0),
     + (-0.9285670518875122,-0.8320432901382446),
     + (-0.29472044110298157,-1.1665705442428589),
     + (0,0),
     + (-1.147657036781311,-0.22732868790626526),
     + (-0.28946158289909363,0.26613736152648926),
     + (-0.2992151081562042,-0.3767027258872986),
     + (0.43568331003189087,0.25014132261276245),
     + (-1.237538456916809,0.6182432770729065),
     + (-0.2242678850889206,-0.17262350022792816),
     + (0.3773956596851349,-2.223900318145752),
     + (0.503607988357544,-0.940649151802063),
     + (1.0857694149017334,-0.11582532525062561),
     + (-0.6909538507461548,-0.8149687051773071),
     + (-1.2845993041992188,0.24226348102092743),
     + (-0.433310329914093,0.06528817862272263),
     + (-0.649471640586853,0.01915639080107212),
     + (0.7267507314682007,0.25733837485313416),
     + (1.151911735534668,0)/

            DATA ((ACLO(I,J),I=1,LDA),J=1,N)/
     + (1.2629542848807933,-0.42951310949188126),
     + (-0.3262333607056494,1.2383041008533804),
     + (1.3297992629225006,-0.2793462818542693),
     + (1.2724293214294047,1.7579030898107073),
     + (-0.2947204467905602,-1.166570547084707),
     + (-0.005767172747536955,-1.0655905803882961),
     + (2.404653388857951,-1.563782051071005),
     + (0.7635934611404596,1.1565369971501793),
     + (-0.29921511789731614,-0.3767027185836281),
     + (-0.411510832795067,2.4413646288945894),
     + (0.2522234481561323,-0.7953391172553718),
     + (-0.8919211272845686,-0.054877473711578625),
     + (0.37739564598170106,-2.2239002740099374),
     + (0.1333363608148414,-1.263614384970583),
     + (0.8041895097449078,0.3587288959713519),
     + (0,0),
     + (-1.2845993538721883,0.24226348085968588),
     + (0.04672617218835198,-1.4250983947324998),
     + (0,0),
     + (0,0),
     + (1.1519117540872,0),
     + (0,0),
     + (0,0),
     + (0,0)/
     
       DATA (V6(I),I=1,6)/
     + (-0.6490100777088978, 0.7721421858045301),
     +  ( -0.11916876241803812, -0.21951562675343952),
     +  ( 0.6641356998941105,-0.4248102833772871),
     +  ( 1.100969102194087,-0.418980099421959),
     +  ( 0.14377148075806995,  0.9969868609091059),
     +  ( -0.11775359816595128 ,-0.27577802908802723)/
     
    
      PRINT * , "==CASE 0======="


        CALL CCPLX(V6,X,6)
        CALL CCPLX(V6,Y,6)

        UPLO='U'
        INCX=1
        INCY=1

        ALPHA = (0.2, 0.8)
        BETA = (0.3, -0.7)
        Y(4) = (0,0)
        Y(2) = (0,0)

        CALL ZHBMV(UPLO,N,K,ALPHA,AC,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y='
        CALL PRNVECC(Y,6)

        PRINT * , "==CASE 1======="


        CALL CCPLX(V6,X,6)
        CALL CCPLX(V6,Y,6)

        UPLO='U'
        INCX=-1
        INCY=-1

        ALPHA = (0.2, 0.8)
        BETA = (0, 0)
        Y(4) = (0,0)
        Y(2) = (0,0)

        CALL ZHBMV(UPLO,N,K,ALPHA,AC,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y='
        CALL PRNVECC(Y,6)

          PRINT * , "==CASE 1======="


        CALL CCPLX(V6,X,6)
        CALL CCPLX(V6,Y,6)

        UPLO='U'
        INCX=-1
        INCY=-1

        ALPHA = (0, 0)
        BETA = (0, 0)
        Y(4) = (0,0)
        Y(2) = (0,0)

        CALL ZHBMV(UPLO,N,K,ALPHA,AC,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y='
        CALL PRNVECC(Y,6)

           PRINT * , "==CASE 4======="


        CALL CCPLX(V6,X,6)
        CALL CCPLX(V6,Y,6)

        UPLO='U'
        INCX=-1
        INCY=-1

        ALPHA = (0.2, 0.8)
        BETA = (1, 0)
        Y(4) = (0,0)
        Y(2) = (0,0)

        CALL ZHBMV(UPLO,N,K,ALPHA,AC,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y='
        CALL PRNVECC(Y,6)

         PRINT * , "==CASE 5======="


        CALL CCPLX(V6,X,6)
        CALL CCPLX(V6,Y,6)

        UPLO='L'
        INCX=1
        INCY=1

        

        ALPHA = (0.2, 0.8)
        BETA = (1, 0)
        Y(4) = (0,0)
        Y(2) = (0,0)

        CALL ZHBMV(UPLO,N,K,ALPHA,ACLO,LDA,X,INCX,BETA,Y,INCY)

        PRINT *,'Y='
        CALL PRNVECC(Y,6)

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
