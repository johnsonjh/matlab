      SUBROUTINE WGEFA(AR,AI,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION AR(LDA,1),AI(LDA,1)
C 
C     WGEFA FACTORS A DOUBLE-COMPLEX MATRIX BY GAUSSIAN ELIMINATION.
C 
C     WGEFA IS USUALLY CALLED BY WGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR WGECO) = (1 + 9/N)*(TIME FOR WGEFA) .
C 
C     ON ENTRY
C 
C        A       DOUBLE-COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C 
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C 
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C 
C     ON RETURN
C 
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C 
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C 
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT WGESL OR WGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN WGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C 
C     LINPACK. THIS VERSION DATED 07/01/79 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C 
C     SUBROUTINES AND FUNCTIONS
C 
C     BLAS WAXPY,WSCAL,IWAMAX
C     FORTRAN DABS
C 
C     INTERNAL VARIABLES
C 
      DOUBLE PRECISION TR,TI
      INTEGER IWAMAX,J,K,KP1,L,NM1
C 
      DOUBLE PRECISION ZDUMR,ZDUMI
      DOUBLE PRECISION CABS1
      CABS1(ZDUMR,ZDUMI) = DABS(ZDUMR) + DABS(ZDUMI)
C 
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C 
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C 
C        FIND L = PIVOT INDEX
C 
         L = IWAMAX(N-K+1,AR(K,K),AI(K,K),1) + K - 1
         IPVT(K) = L
C 
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C 
         IF (CABS1(AR(L,K),AI(L,K)) .EQ. 0.0D0) GO TO 40
C 
C           INTERCHANGE IF NECESSARY
C 
            IF (L .EQ. K) GO TO 10
               TR = AR(L,K)
               TI = AI(L,K)
               AR(L,K) = AR(K,K)
               AI(L,K) = AI(K,K)
               AR(K,K) = TR
               AI(K,K) = TI
   10       CONTINUE
C 
C           COMPUTE MULTIPLIERS
C 
            CALL WDIV(-1.0D0,0.0D0,AR(K,K),AI(K,K),TR,TI)
            CALL WSCAL(N-K,TR,TI,AR(K+1,K),AI(K+1,K),1)
C 
C           ROW ELIMINATION WITH COLUMN INDEXING
C 
            DO 30 J = KP1, N
               TR = AR(L,J)
               TI = AI(L,J)
               IF (L .EQ. K) GO TO 20
                  AR(L,J) = AR(K,J)
                  AI(L,J) = AI(K,J)
                  AR(K,J) = TR
                  AI(K,J) = TI
   20          CONTINUE
               CALL WAXPY(N-K,TR,TI,AR(K+1,K),AI(K+1,K),1,AR(K+1,J),
     *                    AI(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (CABS1(AR(N,N),AI(N,N)) .EQ. 0.0D0) INFO = N
      RETURN
      END
