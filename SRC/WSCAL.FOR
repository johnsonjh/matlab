      SUBROUTINE WSCAL(N,SR,SI,XR,XI,INCX)
      DOUBLE PRECISION SR,SI,XR(1),XI(1)
      IF (N .LE. 0) RETURN
      IX = 1
      DO 10 I = 1, N
         CALL WMUL(SR,SI,XR(IX),XI(IX),XR(IX),XI(IX))
         IX = IX + INCX
   10 CONTINUE
      RETURN
      END
