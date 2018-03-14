      PROGRAM complextest
      INTRINSIC CONJG

      COMPLEX C1, C2, C3
      C1 = CMPLX(1,2)
      C2 = CMPLX(1,1)
      C3 = C1*C2;
      IF (C1.NE.1) THEN
         PRINT *, "c1 looks like", c1
      END IF

      print *, "c1=",c1
      print *, "c2=",c2
      print *, "c3=",C3
      print *, "conj(c3)", CONJG(c3)

      end program complextest
