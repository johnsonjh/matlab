      SUBROUTINE PUTID(X,Y)
C     STORE A NAME
      INTEGER X(4),Y(4)
      DO 10 I = 1, 4
   10 X(I) = Y(I)
      RETURN
      END
