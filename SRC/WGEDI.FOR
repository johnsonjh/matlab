      SUBROUTINE WGEDI(AR,AI,LDA,N,IPVT,DETR,DETI,WORKR,WORKI,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION AR(LDA,1),AI(LDA,1),DETR(2),DETI(2),WORKR(1),
     *                 WORKI(1)
C 
C     WGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY WGECO OR WGEFA.
C 
C     ON ENTRY
C 
C        A       DOUBLE-COMPLEX(LDA, N)
C                THE OUTPUT FROM WGECO OR WGEFA.
C 
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C 
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C 
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM WGECO OR WGEFA.
C 
C        WORK    DOUBLE-COMPLEX(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C 
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C 
C     ON RETURN
C 
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C 
C        DET     DOUBLE-COMPLEX(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. CABS1(DET(1) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C 
C     ERROR CONDITION
C 
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF WGECO HAS SET RCOND .GT. 0.0 OR WGEFA HAS SET
C        INFO .EQ. 0 .
C 
C     LINPACK. THIS VERSION DATED 07/01/79 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C 
C     SUBROUTINES AND FUNCTIONS
C 
C     BLAS WAXPY,WSCAL,WSWAP
C     FORTRAN DABS,MOD
C 
C     INTERNAL VARIABLES
C 
      DOUBLE PRECISION TR,TI
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C 
      DOUBLE PRECISION ZDUMR,ZDUMI
      DOUBLE PRECISION CABS1
      CABS1(ZDUMR,ZDUMI) = DABS(ZDUMR) + DABS(ZDUMI)
C 
C     COMPUTE DETERMINANT
C 
      IF (JOB/10 .EQ. 0) GO TO 80
         DETR(1) = 1.0D0
         DETI(1) = 0.0D0
         DETR(2) = 0.0D0
         DETI(2) = 0.0D0
         TEN = 10.0D0
         DO 60 I = 1, N
            IF (IPVT(I) .EQ. I) GO TO 10
               DETR(1) = -DETR(1)
               DETI(1) = -DETI(1)
   10       CONTINUE
            CALL WMUL(AR(I,I),AI(I,I),DETR(1),DETI(1),DETR(1),DETI(1))
C           ...EXIT
C        ...EXIT
            IF (CABS1(DETR(1),DETI(1)) .EQ. 0.0D0) GO TO 70
   20       IF (CABS1(DETR(1),DETI(1)) .GE. 1.0D0) GO TO 30
               DETR(1) = TEN*DETR(1)
               DETI(1) = TEN*DETI(1)
               DETR(2) = DETR(2) - 1.0D0
               DETI(2) = DETI(2) - 0.0D0
            GO TO 20
   30       CONTINUE
   40       IF (CABS1(DETR(1),DETI(1)) .LT. TEN) GO TO 50
               DETR(1) = DETR(1)/TEN
               DETI(1) = DETI(1)/TEN
               DETR(2) = DETR(2) + 1.0D0
               DETI(2) = DETI(2) + 0.0D0
            GO TO 40
   50       CONTINUE
   60    CONTINUE
   70    CONTINUE
   80 CONTINUE
C 
C     COMPUTE INVERSE(U)
C 
      IF (MOD(JOB,10) .EQ. 0) GO TO 160
         DO 110 K = 1, N
            CALL WDIV(1.0D0,0.0D0,AR(K,K),AI(K,K),AR(K,K),AI(K,K))
            TR = -AR(K,K)
            TI = -AI(K,K)
            CALL WSCAL(K-1,TR,TI,AR(1,K),AI(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 100
            DO 90 J = KP1, N
               TR = AR(K,J)
               TI = AI(K,J)
               AR(K,J) = 0.0D0
               AI(K,J) = 0.0D0
               CALL WAXPY(K,TR,TI,AR(1,K),AI(1,K),1,AR(1,J),AI(1,J),1)
   90       CONTINUE
  100       CONTINUE
  110    CONTINUE
C 
C        FORM INVERSE(U)*INVERSE(L)
C 
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 150
         DO 140 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 120 I = KP1, N
               WORKR(I) = AR(I,K)
               WORKI(I) = AI(I,K)
               AR(I,K) = 0.0D0
               AI(I,K) = 0.0D0
  120       CONTINUE
            DO 130 J = KP1, N
               TR = WORKR(J)
               TI = WORKI(J)
               CALL WAXPY(N,TR,TI,AR(1,J),AI(1,J),1,AR(1,K),AI(1,K),1)
  130       CONTINUE
            L = IPVT(K)
            IF (L .NE. K)
     *         CALL WSWAP(N,AR(1,K),AI(1,K),1,AR(1,L),AI(1,L),1)
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
      RETURN
      END
