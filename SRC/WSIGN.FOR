      SUBROUTINE WSIGN(XR,XI,YR,YI,ZR,ZI)
      DOUBLE PRECISION XR,XI,YR,YI,ZR,ZI,PYTHAG,T
C     IF Y .NE. 0, Z = X*Y/ABS(Y)
C     IF Y .EQ. 0, Z = X
      T = PYTHAG(YR,YI)
      ZR = XR
      ZI = XI
      IF (T .NE. 0.0D0) CALL WMUL(YR/T,YI/T,ZR,ZI,ZR,ZI)
      RETURN
      END
