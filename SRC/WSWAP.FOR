      SUBROUTINE WSWAP(N,XR,XI,INCX,YR,YI,INCY)
      DOUBLE PRECISION XR(1),XI(1),YR(1),YI(1),T
      IF (N .LE. 0) RETURN
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         T = XR(IX)
         XR(IX) = YR(IY)
         YR(IY) = T
         T = XI(IX)
         XI(IX) = YI(IY)
         YI(IY) = T
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
      END
