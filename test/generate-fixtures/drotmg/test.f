      PROGRAM rotmgtest
      
      EXTERNAL DROTMG
c      DROTMG(DD1,DD2,DX1,DY1,DPARAM)

      DOUBLE PRECISION DD1(11),DD2(11),DX1(11),DY1(11)
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5)
*     ..
*     .. Array Arguments ..
*     
      DD1(1)=-4
      DD2(1)=2
      DX1(1)=3
      DY1(1)=9

      DD1(2) = 4
      DD2(2) = 0
      DX1(2)=3
      DY1(2)=9

      DD1(3)=1
      DD2(3)=2
      DX1(3)=3
      DY1(3)=1

      DD1(4)=2
      DD2(4)=-1
      DX1(4)=3
      DY1(4)=8

      DD1(5)=2
      DD2(5)=1
      DX1(5)=3
      DY1(5)=8

      DD1(6)=1.0/(2**24)
      DD2(6)=1.0/(2**25)
      DX1(6)=3
      DY1(6)=8

      DD1(7)=1.0/(2**24)
      DD2(7)=1.0/(2**25)
      DX1(7)=3
      DY1(7)=2


      DD1(8)=2.0/(2**24)
      DD2(8)=1.0/(2**24)
      DX1(8)=3
      DY1(8)=2

      DD1(9) = (2 ** 24)
      DD2(9) = 2 * (2 ** 24)
      DX1(9) = 3
      DY1(9) = 2

      DD1(10) = 2 * (2 ** 24)
      DD2(10) = 1 * (2 ** 24)
      DX1(10) = 3
      DY1(10) = 2

      DO 100, I=1,10
            DPARAM(1)=0
            DPARAM(2)=0
            DPARAM(3)=0
            DPARAM(4)=0
            DPARAM(5)=0
            call DROTMG(DD1(i), DD2(I), DX1(I), DY1(I), DPARAM)
            print *,I,"|", DD1(I), DD2(I), DX1(I), DY1(I),"|",DPARAM
100   CONTINUE
     
      end

     