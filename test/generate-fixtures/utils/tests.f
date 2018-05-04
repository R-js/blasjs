      PROGRAM test

      EXTERNAL COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC, CRCMPLX
      EXTERNAL CCPLX, PRNMATRC,PRNVECC, CCOPYMA
      
      COMPLEX*16 A,B,C


      PRINT *, "=====COMPLEX MULTIPLICATION"
      A=(0.25, -.9)
      B=(0.1,0.5)
      PRINT *, "C=  ( 0.25 - 0.9i )*( 0.1 + 0.5i )", A*B     

      return
   
      end

      SUBROUTINE COPY(X,Y,N)

       INTEGER N
       DOUBLE PRECISION X(*), Y(*)

       DO 20 I=1,N
            Y(I)=X(I)
20     CONTINUE
      END SUBROUTINE

      SUBROUTINE COPYMA(AC,A,N,LDA)

        INTEGER N, LDA
        DOUBLE PRECISION AC(LDA, *), A(LDA,*)

       DO 40 J = 1, N
         DO 50 I = 1, LDA
              A(I,J)=AC(I,J);
50       CONTINUE
40     CONTINUE
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

      SUBROUTINE FILLM(MM,N,M, V)

       INTEGER N,M
       DOUBLE PRECISION MM(M,*),V

       DO 20 J=1,N
         DO I=1,M
            MM(I,J)=V
40       END DO            
20     CONTINUE

      END SUBROUTINE

      SUBROUTINE PRNMATR(A, N, M)

      INTEGER N,M
      DOUBLE PRECISION A(M,*)
      PRINT *, '['
      DO J=1,N
         DO I=1,M
            PRINT *, A(I,J), ","
         END DO
      END DO
     
      PRINT *, ']'

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

