      PROGRAM IZAMAXtest

    
      
      EXTERNAL IZAMAX
c
c     INTEGER FUNCTION IZAMAX(N,SX,INCX)
c
c     SUBROUTINE COPYCC(SRC, DST, N) 
      EXTERNAL COPYCC

c     FILL VECTOR WITH A VALUE      
      EXTERNAL FILL
c
c     PRINT A VECTOR      
      EXTERNAL PRNVECC
c
c 
c
c
      INTEGER RESULT
      INTEGER N, INCX, J

c      original data

      COMPLEX*16  OX(6)

c     working array     
      COMPLEX*16 SX(6) 
      DATA (OX(J),J=1,6)/(1,7),(2,8),(3,9),(4,10),(5,11),(6,12)/
c    == case 0  
      PRINT * , "==CASE 0=======" 
      N=6
      INCX=1
      CALL COPYCC(OX,SX,6)
      SX(2) = (0,0)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVECC(SX, 6)
      RESULT = IZAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT
      PRINT * , "==CASE 1======="    
      N=3
      INCX=2
      CALL COPYCC(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVECC(SX, 6)
      RESULT = IZAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT
      PRINT * , "==CASE 2======="    
      N=6
      INCX=-1
      CALL COPYCC(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVECC(SX, 6)
      RESULT = IZAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      PRINT * , "==CASE 3======="    
      N=0
      INCX=1
      CALL COPYCC(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVECC(SX, 6)
      RESULT = IZAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      PRINT * , "==CASE 4======="    
      N=1
      INCX=1
      CALL COPYCC(OX,SX,6)
      SX(2) = (0,1)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVECC(SX, 6)
      RESULT = IZAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      end
      

      SUBROUTINE COPYCC(X,Y, N)

       INTEGER N
       COMPLEX*16 X(*),Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE
      END SUBROUTINE

      SUBROUTINE PRNVECC(X, N)

      INTEGER N
      COMPLEX*16 X(*)

      PRINT *, '['
      DO J=1,N
            PRINT *, X(J), ", "
      END DO
     
      PRINT *, ']'

      END SUBROUTINE

