      PROGRAM rotmgtest
      
      EXTERNAL DGBMV
c     DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      CHARACTER TRANS
      INTEGER M,N,KL,KU,LDA,INCX,INCY 
      DOUBLE PRECISION ALPHA, BETA, Y(6), X(9)
      DOUBLE PRECISION A(3,9)

      DATA              ((A(I,J),I=1,3),J=1,9)/-1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0,
     +                 -1.0,2.0,-1.0/

      DATA              (Y(I),I=1,6)/1,1,1,1,1,1/
     
      DATA              (X(I),I=1,9)/1,1,2,2,3,3,4,4,5/
     

      ALPHA =1.5
      TRANS = 'N'
      M = 6
      N = 9
      KL = 1
      KU = 1
      BETA = 0
      INCX=1
      INCY=2
      LDA=3
    
      call DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
c           call DROTMG(DD1(i), DD2(I), DX1(I), DY1(I), DPARAM)
      print *,I,"|", Y
     
      end

     