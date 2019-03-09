      PROGRAM IDAMAXtest

    
      
      EXTERNAL IDAMAX
c
c     INTEGER FUNCTION IDAMAX(N,SX,INCX)
c
c     SUBROUTINE COPY(SRC, DST, N) 
      EXTERNAL COPY
c     FILL AN ARRAY WITH ZEROS 
      EXTERNAL ZEROA
c     FILL VECTOR WITH A VALUE      
      EXTERNAL FILL
c
c     PRINT A VECTOR      
      EXTERNAL PRNVEC
c
c 
c
c
      INTEGER RESULT
      INTEGER N, INCX, J

c      original data
      DOUBLE PRECISION  OX(6)

c     working array     
      DOUBLE PRECISION SX(6) 
      DATA (OX(J),J=1,6)/1,2,3,4,5,6/
c    == case 0  
      PRINT * , "==CASE 0======="    
      N=6
      INCX=1
      CALL COPY(OX,SX,6)
      SX(2) = 0
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVEC(SX, 6)
      RESULT = IDAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      PRINT * , "==CASE 1======="    
      N=3
      INCX=2
      CALL COPY(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVEC(SX, 6)
      RESULT = IDAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT
      STOP
      PRINT * , "==CASE 2======="    
      N=6
      INCX=-1
      CALL COPY(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVEC(SX, 6)
      RESULT = IDAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      PRINT * , "==CASE 3======="    
      N=0
      INCX=1
      CALL COPY(OX,SX,6)
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVEC(SX, 6)
      RESULT = IDAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      PRINT * , "==CASE 4======="    
      N=1
      INCX=1
      CALL COPY(OX,SX,6)
      SX(2) = 0
      PRINT *, "N=", N, "INCX=", INCX, "SX="
      CALL PRNVEC(SX, 6)
      RESULT = IDAMAX(N, SX, INCX)
      PRINT *, "RESULT=", RESULT

      end
      
      



      SUBROUTINE COPY(X,Y,N)

       INTEGER N
       DOUBLE PRECISION X(*), Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE
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

      SUBROUTINE PRNVEC(X, N)

      INTEGER N
      DOUBLE PRECISION X(*)

      PRINT *, '['
      DO J=1,N
            PRINT *, X(J), ", "
      END DO
     
      PRINT *, ']'

      END SUBROUTINE

