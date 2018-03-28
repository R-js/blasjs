      PROGRAM complextest
      INTRINSIC CONJG

      COMPLEX C1, C2, C3,c4,c5, beta,X(2)
      LOGICAL T
       INTEGER I,N, INCX
       COMPLEX ZERO,  CV(8,5,2)

      DATA              ((CV(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (0.3E0,-0.4E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (0.1E0,-0.3E0), (0.5E0,-0.1E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (0.1E0,0.1E0),
     +                  (-0.6E0,0.1E0), (0.1E0,-0.3E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (0.3E0,0.1E0), (0.5E0,0.0E0),
     +                  (0.0E0,0.5E0), (0.0E0,0.2E0), (2.0E0,3.0E0),
     +                  (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CV(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (0.3E0,-0.4E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (0.1E0,-0.3E0), (8.0E0,9.0E0), (0.5E0,-0.1E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (0.1E0,0.1E0),
     +                  (3.0E0,6.0E0), (-0.6E0,0.1E0), (4.0E0,7.0E0),
     +                  (0.1E0,-0.3E0), (7.0E0,2.0E0), (7.0E0,2.0E0),
     +                  (7.0E0,2.0E0), (0.3E0,0.1E0), (5.0E0,8.0E0),
     +                  (0.5E0,0.0E0), (6.0E0,9.0E0), (0.0E0,0.5E0),
     +                  (8.0E0,3.0E0), (0.0E0,0.2E0), (9.0E0,4.0E0)/

       PARAMETER (ZERO= (0.0E+0,0.0E+0))
      X(1) =(0.1E0,0.1E0)
      X(2) =(1.0E0,2.0E0)
      N=1
      INCX=1
      


      C1 = CMPLX(1,2)
      C2 = CMPLX(1,1)
      C3 = C1*C2;
      c4 = cmplx(3,2)
      c5 = real(c4);
      beta = cmplx(0,0);

      
     

     
      

      T= BETA.EQ.REAL(ZERO)

      IF (C1.NE.1) THEN
         PRINT *, "c1 looks like", c1
      END IF

      print *, "c1=", c1
      print *, "c2=", c2
      print *, "c3=", C3
      print *, "c4", c4
      print *, "conj(c3)", CONJG(c3)
      print *, "REAL(c3)", REAL(c3)
      print *, "c4*REAL(c3)", c4*REAL(c3)
      print *, "c5", c5
      print *, "BETA.EQ.REAL(ZERO)", T

      DO 99, I=1,4 
            print *, 'I',I
99    CONTINUE
      print *, 'I',I
      print *, 'I',I
      write (6,100) "hello" 
100   FORMAT(1X,A6)      
      print *, "cv(1,2,1)=",cv(1,2,1)
      print *, "SCNRM2(N=1,X,INCX)", SCNRM2(N,X,INCX)
      print *, "SCNRM2(N=2,X,INCX)", SCNRM2(2,X,INCX)
      end

      REAL FUNCTION SCNRM2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
*     ..
*     .. Local Scalars ..
      REAL NORM,SCALE,SSQ,TEMP
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
c         X=(0.1E0,0.1E0), (1.0E0,2.0E0), LEN=2, N=1, STRUE2=0.5       
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (REAL(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(REAL(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
c             SSQ = 1, SCALE= 0.1 
              IF (AIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(AIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
c             SSQ = 1+0.1/0.1 = 2, SCALE = 0.1
   10     CONTINUE
c              0.1*SQRT(2)   
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      SCNRM2 = NORM
      RETURN
*
*     End of SCNRM2.
*
      END
