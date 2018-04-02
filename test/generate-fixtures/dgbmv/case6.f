      PROGRAM rotmgtest
      
      EXTERNAL DGBMV
c     DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

      CHARACTER TRANS
      INTEGER M,N,KL,KU,LDA,INCX,INCY 
      DOUBLE PRECISION ALPHA, BETA, Y(9), X(6)
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

      DATA              (Y(I),I=1,9)/1,1,1,1,1,1,1,1,1/
     
      DATA              (X(I),I=1,6)/1,1,2,2,3,3/
     
     
  
      ALPHA =1.5
      TRANS = 'T'
      M = 6
      N = 9
      KL = 1
      KU = 1
      BETA = 1
      INCX=1
      INCY=1
      LDA=3
    
      call DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
c           call DROTMG(DD1(i), DD2(I), DX1(I), DY1(I), DPARAM)
      print *,I,"|", Y
     
      end

     