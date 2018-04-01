      PROGRAM rotmgtest
      
      EXTERNAL DROTM
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      INTEGER N, INCX, INCY

      DOUBLE PRECISION DX(4),DY(4), DPARAM(5)

c     0
      N = 4
      incx =1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8

      DPARAM(1) =-2
      DPARAM(2) = NaN
      DPARAM(3) = NaN
      DPARAM(4) = NaN
      DPARAM(5) = NaN
c      DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"0|", DX, DY, DPARAM

c     1
      N = 0
      incx =1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8

      DPARAM(1) =-2
      DPARAM(2) = NaN
      DPARAM(3) = NaN
      DPARAM(4) = NaN
      DPARAM(5) = NaN
c      DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"1|", DX, DY, DPARAM      

c     2
      N = 4
      incx =1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8

      DPARAM(1) =-1
      DPARAM(2) = 2
      DPARAM(3) = 3
      DPARAM(4) = 4
      DPARAM(5) = 5
c      DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"2|", DX, DY, DPARAM 

c     3
      N = 4
      incx =1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c 0, NaN, 3, 4, NaN
      DPARAM(1) =0
      DPARAM(2) = NaN
      DPARAM(3) = 3
      DPARAM(4) = 4
      DPARAM(5) = NaN
c      DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"3|", DX, DY, DPARAM       

c     4
      N = 4
      incx =1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c 0, [1, 2, NaN, NaN, 3]
      DPARAM(1) =1
      DPARAM(2) = 2
      DPARAM(3) = NaN
      DPARAM(4) = NaN
      DPARAM(5) = 3
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"4|", DX, DY, DPARAM 

c     5
      N = 4
      incx =-1
      incy =1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c 0, [1, 2, NaN, NaN, 3]
      DPARAM(1) =1
      DPARAM(2) = 2
      DPARAM(3) = NaN
      DPARAM(4) = NaN
      DPARAM(5) = 3
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"5|", DX, DY, DPARAM       
     
c     6
      N = 4
      incx =1
      incy =-1
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c 0, [1, 2, NaN, NaN, 3]
      DPARAM(1) =1
      DPARAM(2) = 2
      DPARAM(3) = NaN
      DPARAM(4) = NaN
      DPARAM(5) = 3
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"6|", DX, DY, DPARAM  

c     7
      N = 2
      incx =2
      incy =-2
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c    sparam: [0, NaN, 4, 5, NaN]
      DPARAM(1) =0
      DPARAM(2) = NaN
      DPARAM(3) = 4
      DPARAM(4) = 5
      DPARAM(5) = NaN
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"7|", DX, DY, DPARAM   

c     8
      N = 2
      incx =2
      incy =-2
c      
      DX(1)=1
      DX(2)=2
      DX(3)=3
      DX(4)=4
c      
      DY(1)=5
      DY(2)=6
      DY(3)=7
      DY(4)=8
c    sparam: [0, NaN, 4, 5, NaN]
      DPARAM(1) =-1
      DPARAM(2) = 1
      DPARAM(3) = 2
      DPARAM(4) = 3
      DPARAM(5) = 4
c     DROTM(N,DX,INCX,DY,INCY,DPARAM)
      call DROTM(N, DX, INCX, DY, INCY, DPARAM)
      print *,"8|", DX, DY, DPARAM   
      end

     