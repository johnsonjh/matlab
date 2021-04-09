      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR,JOB)
C 
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,F,G,P,R,S
      DOUBLE PRECISION FLOP
C 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C 
C     ON INPUT.
C 
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C 
C        N IS THE ORDER OF THE MATRIX.
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C 
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C 
C      ON OUTPUT.
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C 
C        E HAS BEEN DESTROYED.
C 
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C 
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C 
C     ------------------------------------------------------------------
C 
C 
C*****
C     MODIFIED BY C. MOLER TO ELIMINATE MACHEP 11/22/78
C     MODIFIED TO ADD JOB PARAMETER 08/27/79
C*****
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C 
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C 
      E(N) = 0.0D0
C 
      DO 240 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
C*****
            P = FLOP(DABS(D(M)) + DABS(D(M+1)))
            S = FLOP(P + DABS(E(M)))
            IF (P .EQ. S) GO TO 120
C*****
  110    CONTINUE
C 
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = FLOP((D(L+1) - P)/(2.0D0*E(L)))
         R = FLOP(DSQRT(G*G+1.0D0))
         G = FLOP(D(M) - P + E(L)/(G + DSIGN(R,G)))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = FLOP(S*E(I))
            B = FLOP(C*E(I))
            IF (DABS(F) .LT. DABS(G)) GO TO 150
            C = FLOP(G/F)
            R = FLOP(DSQRT(C*C+1.0D0))
            E(I+1) = FLOP(F*R)
            S = FLOP(1.0D0/R)
            C = FLOP(C*S)
            GO TO 160
  150       S = FLOP(F/G)
            R = FLOP(DSQRT(S*S+1.0D0))
            E(I+1) = FLOP(G*R)
            C = FLOP(1.0D0/R)
            S = FLOP(S*C)
  160       G = FLOP(D(I+1) - P)
            R = FLOP((D(I) - G)*S + 2.0D0*C*B)
            P = FLOP(S*R)
            D(I+1) = G + P
            G = FLOP(C*R - B)
            IF (JOB .EQ. 0) GO TO 185
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               F = Z(K,I+1)
               Z(K,I+1) = FLOP(S*Z(K,I) + C*F)
               Z(K,I) = FLOP(C*Z(K,I) - S*F)
  180       CONTINUE
  185       CONTINUE
C 
  200    CONTINUE
C 
         D(L) = FLOP(D(L) - P)
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C 
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C 
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C 
         IF (JOB .EQ. 0) GO TO 285
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  285    CONTINUE
C 
  300 CONTINUE
C 
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
