      SUBROUTINE WGECO(AR,AI,LDA,N,IPVT,RCOND,ZR,ZI)
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION AR(LDA,1),AI(LDA,1),ZR(1),ZI(1)
      DOUBLE PRECISION RCOND
C 
C     WGECO FACTORS A DOUBLE-COMPLEX MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C 
C     IF  RCOND  IS NOT NEEDED, WGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW WGECO BY WGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW WGECO BY WGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW WGECO BY WGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW WGECO BY WGEDI.
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
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C 
C        Z       DOUBLE-COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C 
C     LINPACK. THIS VERSION DATED 07/01/79 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C 
C     SUBROUTINES AND FUNCTIONS
C 
C     LINPACK WGEFA
C     BLAS WAXPY,WDOTC,WASUM
C     FORTRAN DABS,DMAX1
C 
C     INTERNAL VARIABLES
C 
      DOUBLE PRECISION WDOTCR,WDOTCI,EKR,EKI,TR,TI,WKR,WKI,WKMR,WKMI
      DOUBLE PRECISION ANORM,S,WASUM,SM,YNORM,FLOP
      INTEGER INFO,J,K,KB,KP1,L
C 
      DOUBLE PRECISION ZDUMR,ZDUMI
      DOUBLE PRECISION CABS1
      CABS1(ZDUMR,ZDUMI) = DABS(ZDUMR) + DABS(ZDUMI)
C 
C     COMPUTE 1-NORM OF A
C 
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,WASUM(N,AR(1,J),AI(1,J),1))
   10 CONTINUE
C 
C     FACTOR
C 
      CALL WGEFA(AR,AI,LDA,N,IPVT,INFO)
C 
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C 
C     SOLVE CTRANS(U)*W = E
C 
      EKR = 1.0D0
      EKI = 0.0D0
      DO 20 J = 1, N
         ZR(J) = 0.0D0
         ZI(J) = 0.0D0
   20 CONTINUE
      DO 110 K = 1, N
         CALL WSIGN(EKR,EKI,-ZR(K),-ZI(K),EKR,EKI)
         IF (CABS1(EKR-ZR(K),EKI-ZI(K))
     *       .LE. CABS1(AR(K,K),AI(K,K))) GO TO 40
            S = CABS1(AR(K,K),AI(K,K))
     *          /CABS1(EKR-ZR(K),EKI-ZI(K))
            CALL WRSCAL(N,S,ZR,ZI,1)
            EKR = S*EKR
            EKI = S*EKI
   40    CONTINUE
         WKR = EKR - ZR(K)
         WKI = EKI - ZI(K)
         WKMR = -EKR - ZR(K)
         WKMI = -EKI - ZI(K)
         S = CABS1(WKR,WKI)
         SM = CABS1(WKMR,WKMI)
         IF (CABS1(AR(K,K),AI(K,K)) .EQ. 0.0D0) GO TO 50
            CALL WDIV(WKR,WKI,AR(K,K),-AI(K,K),WKR,WKI)
            CALL WDIV(WKMR,WKMI,AR(K,K),-AI(K,K),WKMR,WKMI)
         GO TO 60
   50    CONTINUE
            WKR = 1.0D0
            WKI = 0.0D0
            WKMR = 1.0D0
            WKMI = 0.0D0
   60    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 100
            DO 70 J = KP1, N
               CALL WMUL(WKMR,WKMI,AR(K,J),-AI(K,J),TR,TI)
               SM = FLOP(SM + CABS1(ZR(J)+TR,ZI(J)+TI))
               CALL WAXPY(1,WKR,WKI,AR(K,J),-AI(K,J),1,
     $                    ZR(J),ZI(J),1)
               S = FLOP(S + CABS1(ZR(J),ZI(J)))
   70       CONTINUE
            IF (S .GE. SM) GO TO 90
               TR = WKMR - WKR
               TI = WKMI - WKI
               WKR = WKMR
               WKI = WKMI
               DO 80 J = KP1, N
                  CALL WAXPY(1,TR,TI,AR(K,J),-AI(K,J),1,
     $                       ZR(J),ZI(J),1)
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
         ZR(K) = WKR
         ZI(K) = WKI
  110 CONTINUE
      S = 1.0D0/WASUM(N,ZR,ZI,1)
      CALL WRSCAL(N,S,ZR,ZI,1)
C 
C     SOLVE CTRANS(L)*Y = W
C 
      DO 140 KB = 1, N
         K = N + 1 - KB
         IF (K .GE. N) GO TO 120
            ZR(K) = ZR(K)
     *            + WDOTCR(N-K,AR(K+1,K),AI(K+1,K),1,ZR(K+1),ZI(K+1),1)
            ZI(K) = ZI(K)
     *            + WDOTCI(N-K,AR(K+1,K),AI(K+1,K),1,ZR(K+1),ZI(K+1),1)
  120    CONTINUE
         IF (CABS1(ZR(K),ZI(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/CABS1(ZR(K),ZI(K))
            CALL WRSCAL(N,S,ZR,ZI,1)
  130    CONTINUE
         L = IPVT(K)
         TR = ZR(L)
         TI = ZI(L)
         ZR(L) = ZR(K)
         ZI(L) = ZI(K)
         ZR(K) = TR
         ZI(K) = TI
  140 CONTINUE
      S = 1.0D0/WASUM(N,ZR,ZI,1)
      CALL WRSCAL(N,S,ZR,ZI,1)
C 
      YNORM = 1.0D0
C 
C     SOLVE L*V = Y
C 
      DO 160 K = 1, N
         L = IPVT(K)
         TR = ZR(L)
         TI = ZI(L)
         ZR(L) = ZR(K)
         ZI(L) = ZI(K)
         ZR(K) = TR
         ZI(K) = TI
         IF (K .LT. N)
     *      CALL WAXPY(N-K,TR,TI,AR(K+1,K),AI(K+1,K),1,ZR(K+1),ZI(K+1),
     *                 1)
         IF (CABS1(ZR(K),ZI(K)) .LE. 1.0D0) GO TO 150
            S = 1.0D0/CABS1(ZR(K),ZI(K))
            CALL WRSCAL(N,S,ZR,ZI,1)
            YNORM = S*YNORM
  150    CONTINUE
  160 CONTINUE
      S = 1.0D0/WASUM(N,ZR,ZI,1)
      CALL WRSCAL(N,S,ZR,ZI,1)
      YNORM = S*YNORM
C 
C     SOLVE  U*Z = V
C 
      DO 200 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(ZR(K),ZI(K))
     *       .LE. CABS1(AR(K,K),AI(K,K))) GO TO 170
            S = CABS1(AR(K,K),AI(K,K))
     *          /CABS1(ZR(K),ZI(K))
            CALL WRSCAL(N,S,ZR,ZI,1)
            YNORM = S*YNORM
  170    CONTINUE
         IF (CABS1(AR(K,K),AI(K,K)) .EQ. 0.0D0) GO TO 180
            CALL WDIV(ZR(K),ZI(K),AR(K,K),AI(K,K),ZR(K),ZI(K))
  180    CONTINUE
         IF (CABS1(AR(K,K),AI(K,K)) .NE. 0.0D0) GO TO 190
            ZR(K) = 1.0D0
            ZI(K) = 0.0D0
  190    CONTINUE
         TR = -ZR(K)
         TI = -ZI(K)
         CALL WAXPY(K-1,TR,TI,AR(1,K),AI(1,K),1,ZR(1),ZI(1),1)
  200 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/WASUM(N,ZR,ZI,1)
      CALL WRSCAL(N,S,ZR,ZI,1)
      YNORM = S*YNORM
C 
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
