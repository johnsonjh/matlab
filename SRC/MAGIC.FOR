      SUBROUTINE MAGIC(A,LDA,N)
C 
C     ALGORITHMS FOR MAGIC SQUARES TAKEN FROM
C        MATHEMATICAL RECREATIONS AND ESSAYS, 12TH ED.,
C        BY W. W. ROUSE BALL AND H. S. M. COXETER
C 
      DOUBLE PRECISION A(LDA,N),T
C 
      IF (MOD(N,4) .EQ. 0) GO TO 100
      IF (MOD(N,2) .EQ. 0) M = N/2
      IF (MOD(N,2) .NE. 0) M = N
C 
C     ODD ORDER OR UPPER CORNER OF EVEN ORDER
C 
      DO 20 J = 1,M
         DO 10 I = 1,M
            A(I,J) = 0
   10    CONTINUE
   20 CONTINUE
      I = 1
      J = (M+1)/2
      MM = M*M
      DO 40 K = 1, MM
         A(I,J) = K
         I1 = I-1
         J1 = J+1
         IF(I1.LT.1) I1 = M
         IF(J1.GT.M) J1 = 1
         IF(IDINT(A(I1,J1)).EQ.0) GO TO 30
            I1 = I+1
            J1 = J
   30    I = I1
         J = J1
   40 CONTINUE
      IF (MOD(N,2) .NE. 0) RETURN
C 
C     REST OF EVEN ORDER
C 
      T = M*M
      DO 60 I = 1, M
         DO 50 J = 1, M
            IM = I+M
            JM = J+M
            A(I,JM) = A(I,J) + 2*T
            A(IM,J) = A(I,J) + 3*T
            A(IM,JM) = A(I,J) + T
   50    CONTINUE
   60 CONTINUE
      M1 = (M-1)/2
      IF (M1.EQ.0) RETURN
      DO 70 J = 1, M1
         CALL RSWAP(M,A(1,J),1,A(M+1,J),1)
   70 CONTINUE
      M1 = (M+1)/2
      M2 = M1 + M
      CALL RSWAP(1,A(M1,1),1,A(M2,1),1)
      CALL RSWAP(1,A(M1,M1),1,A(M2,M1),1)
      M1 = N+1-(M-3)/2
      IF(M1.GT.N) RETURN
      DO 80 J = M1, N
         CALL RSWAP(M,A(1,J),1,A(M+1,J),1)
   80 CONTINUE
      RETURN
C 
C     DOUBLE EVEN ORDER
C 
  100 K = 1
      DO 120 I = 1, N
         DO 110 J = 1, N
            A(I,J) = K
            IF (MOD(I,4)/2 .EQ. MOD(J,4)/2) A(I,J) = N*N+1 - K
            K = K+1
  110    CONTINUE
  120 CONTINUE
      RETURN
      END
