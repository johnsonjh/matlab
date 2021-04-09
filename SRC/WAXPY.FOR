      SUBROUTINE WAXPY(N,SR,SI,XR,XI,INCX,YR,YI,INCY)
      DOUBLE PRECISION SR,SI,XR(1),XI(1),YR(1),YI(1),FLOP
      IF (N .LE. 0) RETURN
      IF (SR .EQ. 0.0D0 .AND. SI .EQ. 0.0D0) RETURN
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         YR(IY) = FLOP(YR(IY) + SR*XR(IX) - SI*XI(IX))
         YI(IY) = YI(IY) + SR*XI(IX) + SI*XR(IX)
         IF (YI(IY) .NE. 0.0D0) YI(IY) = FLOP(YI(IY))
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
      END
