      SUBROUTINE RSWAP(N,X,INCX,Y,INCY)
      DOUBLE PRECISION X(1),Y(1),T
      IF (N .LE. 0) RETURN
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I = 1, N
         T = X(IX)
         X(IX) = Y(IY)
         Y(IY) = T
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
      END
