      SUBROUTINE RROT(N,DX,INCX,DY,INCY,C,S)
C 
C     APPLIES A PLANE ROTATION.
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S,FLOP
      INTEGER I,INCX,INCY,IX,IY,N
C 
      IF (N.LE.0) RETURN
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = FLOP(C*DX(IX) + S*DY(IY))
        DY(IY) = FLOP(C*DY(IY) - S*DX(IX))
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
