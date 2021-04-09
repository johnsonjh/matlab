      SUBROUTINE RROTG(DA,DB,C,S)
C 
C     CONSTRUCT GIVENS PLANE ROTATION.
C 
      DOUBLE PRECISION DA,DB,C,S,RHO,PYTHAG,FLOP,R,Z
C 
      RHO = DB
      IF ( DABS(DA) .GT. DABS(DB) ) RHO = DA
      C = 1.0D0
      S = 0.0D0
      Z = 1.0D0
      R = FLOP(DSIGN(PYTHAG(DA,DB),RHO))
      IF (R .NE. 0.0D0) C = FLOP(DA/R)
      IF (R .NE. 0.0D0) S = FLOP(DB/R)
      IF ( DABS(DA) .GT. DABS(DB) ) Z = S
      IF ( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = FLOP(1.0D0/C)
      DA = R
      DB = Z
      RETURN
      END
