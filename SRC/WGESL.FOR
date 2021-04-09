      SUBROUTINE WGESL(AR,AI,LDA,N,IPVT,BR,BI,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION AR(LDA,1),AI(LDA,1),BR(1),BI(1)
C 
C     WGESL SOLVES THE DOUBLE-COMPLEX SYSTEM
C     A * X = B  OR  CTRANS(A) * X = B
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
C        B       DOUBLE-COMPLEX(N)
C                THE RIGHT HAND SIDE VECTOR.
C 
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE
C                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.
C 
C     ON RETURN
C 
C        B       THE SOLUTION VECTOR  X .
C 
C     ERROR CONDITION
C 
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF WGECO HAS SET RCOND .GT. 0.0
C        OR WGEFA HAS SET INFO .EQ. 0 .
C 
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL WGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL WGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C 
C     LINPACK. THIS VERSION DATED 07/01/79 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C 
C     SUBROUTINES AND FUNCTIONS
C 
C     BLAS WAXPY,WDOTC
C 
C     INTERNAL VARIABLES
C 
      DOUBLE PRECISION WDOTCR,WDOTCI,TR,TI
      INTEGER K,KB,L,NM1
C 
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C 
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C 
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            TR = BR(L)
            TI = BI(L)
            IF (L .EQ. K) GO TO 10
               BR(L) = BR(K)
               BI(L) = BI(K)
               BR(K) = TR
               BI(K) = TI
   10       CONTINUE
            CALL WAXPY(N-K,TR,TI,AR(K+1,K),AI(K+1,K),1,BR(K+1),BI(K+1),
     *                 1)
   20    CONTINUE
   30    CONTINUE
C 
C        NOW SOLVE  U*X = Y
C 
         DO 40 KB = 1, N
            K = N + 1 - KB
            CALL WDIV(BR(K),BI(K),AR(K,K),AI(K,K),BR(K),BI(K))
            TR = -BR(K)
            TI = -BI(K)
            CALL WAXPY(K-1,TR,TI,AR(1,K),AI(1,K),1,BR(1),BI(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C 
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
C        FIRST SOLVE  CTRANS(U)*Y = B
C 
         DO 60 K = 1, N
            TR = BR(K) - WDOTCR(K-1,AR(1,K),AI(1,K),1,BR(1),BI(1),1)
            TI = BI(K) - WDOTCI(K-1,AR(1,K),AI(1,K),1,BR(1),BI(1),1)
            CALL WDIV(TR,TI,AR(K,K),-AI(K,K),BR(K),BI(K))
   60    CONTINUE
C 
C        NOW SOLVE CTRANS(L)*X = Y
C 
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            BR(K) = BR(K)
     *            + WDOTCR(N-K,AR(K+1,K),AI(K+1,K),1,BR(K+1),BI(K+1),1)
            BI(K) = BI(K)
     *            + WDOTCI(N-K,AR(K+1,K),AI(K+1,K),1,BR(K+1),BI(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               TR = BR(L)
               TI = BI(L)
               BR(L) = BR(K)
               BI(L) = BI(K)
               BR(K) = TR
               BI(K) = TI
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
