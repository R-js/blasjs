      PROGRAM INOUT
C
C This program reads in and prints out a name
C
      INTEGER ADD
      EXTERNAL ADD

      REAL VALUE
      CHARACTER NAME*20
      CHARACTER POST*20, LAST*20
      REAL A, B, C, D, E
      DIMENSION VALUE(0:100)

c        INTEGER X1,X2
      ADD2(X1,X2) = X1+X2
    
      
      PRINT *, " Type in your name, up to 20 characters"
      PRINT *, ' enclosed in quotes'
      READ *, NAME
      PRINT *, NAME
    
      A = 0.12345678901234567
      B = A
      C = A
      D = A 
      E = A 
      DO 456,  I=0,100
         VALUE(I)= ADD2(I*144.0,4*1000.0);
         VALUE(I)=VALUE(I)/1000.0
456   CONTINUE         

C      OPEN (UNIT=123,FILE='DATA.TXT')
      POST = "hello"
      LAST = "dingbats"
      WRITE(*, FMT=101) (VALUE(I), POST, I=0,100)
101   FORMAT(F10.4,/, A20)

C      CLOSE(123)
      END 