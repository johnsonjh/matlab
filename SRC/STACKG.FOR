      SUBROUTINE STACKG(ID)
      INTEGER ID(4)
C
C     GET VARIABLES FROM STORAGE
C
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      LOGICAL EQID
      IF (DDT .EQ. 1) WRITE(WTE,100) ID
  100 FORMAT(1X,'STACKG',4I4)
      K = LSIZE+1
   10 K = K-1
      IF (K .EQ. BOT-1) GO TO 98
      IF (.NOT. EQID(IDSTK(1,K), ID(1))) GO TO 10
      IF (K .GE. LSIZE-1 .AND. RHS .GT. 0) GO TO 98
      LK = LSTK(K)
      IF (RHS .EQ. 1) GO TO 40
      IF (RHS .EQ. 2) GO TO 60
      IF (RHS .GT. 2) CALL ERROR(21)
      IF (ERR .GT. 0) RETURN
      L = 1
      IF (TOP .GT. 0) L = LSTK(TOP) + MSTK(TOP)*NSTK(TOP)
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
C
C     LOAD VARIABLE TO TOP OF STACK
      LSTK(TOP) = L
      MSTK(TOP) = MSTK(K)
      NSTK(TOP) = NSTK(K)
      MN = MSTK(K)*NSTK(K)
      ERR = L+MN - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
C     IF RAND, MATFN6 GENERATES RANDOM NUMBER
      IF (K .EQ. LSIZE) GO TO 97
      CALL WCOPY(MN,STKR(LK),STKI(LK),1,STKR(L),STKI(L),1)
      GO TO 99
C
C     VECT(ARG)
   40 IF (MSTK(TOP) .EQ. 0) GO TO 99
      L = LSTK(TOP)
      MN = MSTK(TOP)*NSTK(TOP)
      MNK = MSTK(K)*NSTK(K)
      IF (MSTK(TOP) .LT. 0) MN = MNK
      DO 50 I = 1, MN
        LL = L+I-1
        LS = LK+I-1
        IF (MSTK(TOP) .GT. 0) LS = LK + IDINT(STKR(LL)) - 1
        IF (LS .LT. LK .OR. LS .GE. LK+MNK) CALL ERROR(21)
        IF (ERR .GT. 0) RETURN
        STKR(LL) = STKR(LS)
        STKI(LL) = STKI(LS)
   50 CONTINUE
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      IF (MSTK(K) .GT. 1) MSTK(TOP) = MN
      IF (MSTK(K) .EQ. 1) NSTK(TOP) = MN
      GO TO 99
C 
C     MATRIX(ARG,ARG)
   60 TOP = TOP-1
      L = LSTK(TOP)
      IF (MSTK(TOP+1) .EQ. 0) MSTK(TOP) = 0
      IF (MSTK(TOP) .EQ. 0) GO TO 99
      L2 = LSTK(TOP+1)
      M = MSTK(TOP)*NSTK(TOP)
      IF (MSTK(TOP) .LT. 0) M = MSTK(K)
      N = MSTK(TOP+1)*NSTK(TOP+1)
      IF (MSTK(TOP+1) .LT. 0) N = NSTK(K)
      L3 = L2 + N
      MK = MSTK(K)
      MNK = MSTK(K)*NSTK(K)
      DO 70 J = 1, N
      DO 70 I = 1, M
        LI = L+I-1
        IF (MSTK(TOP) .GT. 0) LI = L + IDINT(STKR(LI)) - 1
        LJ = L2+J-1
        IF (MSTK(TOP+1) .GT. 0) LJ = L2 + IDINT(STKR(LJ)) - 1
        LS = LK + LI-L + (LJ-L2)*MK
        IF (LS.LT.LK .OR. LS.GE.LK+MNK) CALL ERROR(21)
        IF (ERR .GT. 0) RETURN
        LL = L3 + I-1 + (J-1)*M
        STKR(LL) = STKR(LS)
        STKI(LL) = STKI(LS)
   70 CONTINUE
      MN = M*N
      CALL WCOPY(MN,STKR(L3),STKI(L3),1,STKR(L),STKI(L),1)
      MSTK(TOP) = M
      NSTK(TOP) = N
      GO TO 99
   97 FIN = 7
      FUN = 6
      RETURN
   98 FIN = 0
      RETURN
   99 FIN = -1
      FUN = 0
      RETURN
      END
