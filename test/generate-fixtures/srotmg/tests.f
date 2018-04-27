      PROGRAM dsymvtest
c      
c     x := A*x,   or   x := A**T*x,
c
c     HELPERS
c
c
*       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DD1,DD2,DX1,DY1
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DPARAM(5)
*       ..

*> \param[in,out] DD1
*> \verbatim
*>          DD1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in,out] DD2
*> \verbatim
*>          DD2 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in,out] DX1

      EXTERNAL DROTMG, COPY, ZEROA, COPYMA,FILL
      EXTERNAL FILLM, PRNMATR, PRNVEC

      DOUBLE PRECISION DD1,DD2,DX1,DY1
      DOUBLE PRECISION DPARAM(5), DCPARAM(5)


        DATA (DCPARAM(I),I=1,5)/.0,.0,.0,.0,.0/
      
    
      PRINT * , "==CASE 6======="

     
       CALL COPY(DCPARAM,DPARAM,5)

    
        DD1 = 5.960464477539063e-8
        DD2 = 2.9802322387695312e-8
        DX1 = 3.0
        DY1 = 2.0

        CALL DROTMG(DD1,DD2,DX1,DY1,DPARAM)
     
        PRINT *,"X="
        PRINT *, "DD1",DD1, "DD2",DD2, "DX1",DX1, "DY1", DY1
        CALL PRNVEC(DPARAM,5)


        PRINT * , "==CASE 8======="

     
        CALL COPY(DCPARAM,DPARAM,5)

    
        DD1 = 16777216
        DD2 = 33554432
        DX1 = 3.0
        DY1 = 2.0

        CALL DROTMG(DD1,DD2,DX1,DY1,DPARAM)
     
        PRINT *,"X="
        PRINT *, "DD1",DD1, "DD2",DD2, "DX1",DX1, "DY1", DY1
        CALL PRNVEC(DPARAM,5)

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

