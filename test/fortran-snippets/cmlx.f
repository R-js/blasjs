      PROGRAM complextest
      INTRINSIC CONJG
      EXTERNAL CROTG, ZROTG

      COMPLEX C1, C2, C3,c4,c5, beta,X(2), CA, CB,S
      REAL C
      LOGICAL T
       INTEGER I,N, INCX
       COMPLEX ZERO,  CV(8,5,2)
      DOUBLE COMPLEX DCA, DCB, DS
      REAL*8 DC


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

      CB = CMPLX(34,23)
      CA = CMPLX(11,19)

      DCB = CMPLX(34,23)
      DCA = CMPLX(11,19)
      DS = CMPLX(0,0)
      DC = 0


      print *, " before  CA,CB,C,S", CA,CB,C,S
      CALL CROTG(CA,CB,C,S)
      print *, " after  CA,CB,C,S", CA,CB,C,S

c      print *, " before  DCA,DCB,DC,DS", DCA,DCB,DC,DS
c      CALL ZROTG(DCA,DCB,DC,DS)
c      print *, " after  DCA,DCB,DC,DS", DCA,DCB,DC,DS

     

     
      

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
     
   
     
      end

     