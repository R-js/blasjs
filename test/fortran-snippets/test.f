      PROGRAM INOUT
C
C This program reads in and prints out a name
C
      INTEGER ADD
      EXTERNAL ADD

      INTEGER I, LB
      REAL VALUE
      CHARACTER NAME*20
      CHARACTER POST*20, LAST*20
      REAL A, B, C, D, E
      DIMENSION LB(2)
      DIMENSION VALUE(10, -5:5)

      LB = LBOUND(VALUE)

      PRINT *, "LOWERBOUND", LB(2)
      PRINT *, "UPPERBOUND", UBOUND(VALUE)
c      INTEGER X1,X2
c      ADD2(X1,X2)
      RETURN
    
      
      PRINT *, " Type in your name, up to 20 characters"
      PRINT *, ' enclosed in quotes'
      READ *, NAME
      PRINT *, NAME
    
      A = 0.12345678901234567
      B = A
      C = A
      D = A 
      E = A 
c      DO 456,  I=0,100
c         VALUE(I)= ADD(I*144.0,4*1000.0);
c        VALUE(I)=VALUE(I)/1000.0
c 456   CONTINUE         

C      OPEN (UNIT=123,FILE='DATA.TXT')
      POST = "hello"
      LAST = "dingbats"
c      WRITE(*, FMT=101) (VALUE(I), POST, I=0,100)
c 101   FORMAT(F10.4,/, A20)

C      CLOSE(123)
      END 