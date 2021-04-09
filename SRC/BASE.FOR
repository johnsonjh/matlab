      SUBROUTINE BASE(X,B,EPS,S,N)
      DOUBLE PRECISION X,B,EPS,S(1),T
C 
C     STORE BASE B REPRESENTATION OF X IN S(1:N)
C 
      INTEGER PLUS,MINUS,DOT,ZERO,COMMA
      DATA PLUS/41/,MINUS/42/,DOT/47/,ZERO/0/,COMMA/48/
      L = 1
      IF (X .GE. 0.0D0) S(L) = PLUS
      IF (X .LT. 0.0D0) S(L) = MINUS
      S(L+1) = ZERO
      S(L+2) = DOT
      X = DABS(X)
      IF (X .NE. 0.0D0) K = DLOG(X)/DLOG(B)
      IF (X .EQ. 0.0D0) K = 0
      IF (X .GT. 1.0D0) K = K + 1
      X = X/B**K
      IF (B*X .GE. B) K = K + 1
      IF (B*X .GE. B) X = X/B
      IF (EPS .NE. 0.0D0) M = -DLOG(EPS)/DLOG(B) + 4
      IF (EPS .EQ. 0.0D0) M = 54
      DO 10 L = 4, M
      X = B*X
      J = IDINT(X)
      S(L) = DBLE(J)
      X = X - S(L)
   10 CONTINUE
      S(M+1) = COMMA
      IF (K .GE. 0) S(M+2) = PLUS
      IF (K .LT. 0) S(M+2) = MINUS
      T = DABS(DBLE(K))
      N = M + 3
      IF (T .GE. B) N = N + IDINT(DLOG(T)/DLOG(B))
      L = N
   20 J = IDINT(DMOD(T,B))
      S(L) = DBLE(J)
      L = L - 1
      T = T/B
      IF (L .GE. M+3) GO TO 20
      RETURN
      END
