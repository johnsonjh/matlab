      SUBROUTINE WATAN(XR,XI,YR,YI)
C     Y = ATAN(X) = (I/2)*LOG((I+X)/(I-X))
      DOUBLE PRECISION XR,XI,YR,YI,TR,TI
      IF (XI .NE. 0.0D0) GO TO 10
         YR = DATAN2(XR,1.0D0)
         YI = 0.0D0
         RETURN
   10 IF (XR.NE.0.0D0 .OR. DABS(XI).NE.1.0D0) GO TO 20
         CALL ERROR(32)
         RETURN
   20 CALL WDIV(XR,1.0D0+XI,-XR,1.0D0-XI,TR,TI)
      CALL WLOG(TR,TI,TR,TI)
      YR = -TI/2.0D0
      YI = TR/2.0D0
      RETURN
      END
