      SUBROUTINE WSQRT(XR,XI,YR,YI)
      DOUBLE PRECISION XR,XI,YR,YI,S,TR,TI,PYTHAG,FLOP
C     Y = SQRT(X) WITH YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C 
      TR = XR
      TI = XI
      S = DSQRT(0.5D0*(PYTHAG(TR,TI) + DABS(TR)))
      IF (TR .GE. 0.0D0) YR = FLOP(S)
      IF (TI .LT. 0.0D0) S = -S
      IF (TR .LE. 0.0D0) YI = FLOP(S)
      IF (TR .LT. 0.0D0) YR = FLOP(0.5D0*(TI/YI))
      IF (TR .GT. 0.0D0) YI = FLOP(0.5D0*(TI/YR))
      RETURN
      END
