      SUBROUTINE HILBER(A,LDA,N)
      DOUBLE PRECISION A(LDA,N)
C     GENERATE INVERSE HILBERT MATRIX
      DOUBLE PRECISION P,R
      P = DBLE(N)
      DO 20 I = 1, N
        IF (I.NE.1) P = (DBLE(N-I+1)*P*DBLE(N+I-1))/DBLE(I-1)**2
        R = P*P
        A(I,I) = R/DBLE(2*I-1)
        IF (I.EQ.N) GO TO 20
        IP1 = I+1
        DO 10 J = IP1, N
          R = -(DBLE(N-J+1)*R*(N+J-1))/DBLE(J-1)**2
          A(I,J) = R/DBLE(I+J-1)
          A(J,I) = A(I,J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
