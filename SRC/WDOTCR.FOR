      DOUBLE PRECISION FUNCTION WDOTCR(N,XR,XI,INCX,YR,YI,INCY)
      DOUBLE PRECISION XR(1),XI(1),YR(1),YI(1),S,FLOP
      S = 0.0D0
      IF (N .LE. 0) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         S = FLOP(S + XR(IX)*YR(IY) + XI(IX)*YI(IY))
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
   20 WDOTCR = S
      RETURN
      END
