      SUBROUTINE WSET(N,XR,XI,YR,YI,INCY)
      INTEGER N,INCY
      DOUBLE PRECISION XR,XI,YR(1),YI(1)
      IY = 1
      IF (N .LE. 0 ) RETURN
      DO 10 I = 1,N
         YR(IY) = XR
         YI(IY) = XI
         IY = IY + INCY
   10 CONTINUE
      RETURN
      END
