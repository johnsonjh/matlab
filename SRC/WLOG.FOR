      SUBROUTINE WLOG(XR,XI,YR,YI)
      DOUBLE PRECISION XR,XI,YR,YI,T,R,PYTHAG
C     Y = LOG(X)
      R = PYTHAG(XR,XI)
      IF (R .EQ. 0.0D0) CALL ERROR(32)
      IF (R .EQ. 0.0D0) RETURN
      T = DATAN2(XI,XR)
      IF (XI.EQ.0.0D0 .AND. XR.LT.0.0D0) T = DABS(T)
      YR = DLOG(R)
      YI = T
      RETURN
      END
