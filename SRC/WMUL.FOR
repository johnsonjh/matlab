      SUBROUTINE WMUL(AR,AI,BR,BI,CR,CI)
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI,T,FLOP
C     C = A*B
      T = AR*BI + AI*BR
      IF (T .NE. 0.0D0) T = FLOP(T)
      CR = FLOP(AR*BR - AI*BI)
      CI = T
      RETURN
      END
