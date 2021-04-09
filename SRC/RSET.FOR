      SUBROUTINE RSET(N,DX,DY,INCY)
C 
C     COPIES A SCALAR, X, TO A SCALAR, Y.
      DOUBLE PRECISION DX,DY(1)
C 
      IF (N.LE.0) RETURN
      IY = 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
