      LOGICAL FUNCTION EQID(X,Y)
C     CHECK FOR EQUALITY OF TWO NAMES
      INTEGER X(4),Y(4)
      EQID = .TRUE.
      DO 10 I = 1, 4
   10 EQID = EQID .AND. (X(I).EQ.Y(I))
      RETURN
      END
