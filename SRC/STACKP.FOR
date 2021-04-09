      SUBROUTINE STACKP(ID)
      INTEGER ID(4)
C 
C     PUT VARIABLES INTO STORAGE
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      LOGICAL EQID
      INTEGER SEMI
      DATA SEMI/39/
      IF (DDT .EQ. 1) WRITE(WTE,100) ID
  100 FORMAT(1X,'STACKP',4I4)
      IF (TOP .LE. 0) CALL ERROR(1)
      IF (ERR .GT. 0) RETURN
      CALL FUNS(ID)
      IF (FIN .NE. 0) CALL ERROR(25)
      IF (ERR .GT. 0) RETURN
      M = MSTK(TOP)
      N = NSTK(TOP)
      IF (M .GT. 0) L = LSTK(TOP)
      IF (M .LT. 0) CALL ERROR(14)
      IF (ERR .GT. 0) RETURN
      IF (M .EQ. 0 .AND. N .NE. 0) GO TO 99
      MN = M*N
      LK = 0
      MK = 1
      NK = 0
      LT = 0
      MT = 0
      NT = 0
C 
C     DOES VARIABLE ALREADY EXIST
      CALL PUTID(IDSTK(1,BOT-1),ID)
      K = LSIZE+1
   05 K = K-1
      IF (.NOT.EQID(IDSTK(1,K),ID)) GO TO 05
      IF (K .EQ. BOT-1) GO TO 30
      LK = LSTK(K)
      MK = MSTK(K)
      NK = NSTK(K)
      MNK = MK*NK
      IF (RHS .EQ. 0) GO TO 20
      IF (RHS .GT. 2) CALL ERROR(15)
      IF (ERR .GT. 0) RETURN
      MT = MK
      NT = NK
      LT = L + MN
      ERR = LT + MNK - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WCOPY(MNK,STKR(LK),STKI(LK),1,STKR(LT),STKI(LT),1)
C 
C     DOES IT FIT
   20 IF (RHS.EQ.0 .AND. MN.EQ.MNK) GO TO 40
      IF (K .GE. LSIZE-3) CALL ERROR(13)
      IF (ERR .GT. 0) RETURN
C 
C     SHIFT STORAGE
      IF (K .EQ. BOT) GO TO 25
      LS = LSTK(BOT)
      LL = LS + MNK
      CALL WCOPY(LK-LS,STKR(LS),STKI(LS),-1,STKR(LL),STKI(LL),-1)
      KM1 = K-1
      DO 24 IB = BOT, KM1
        I = BOT+KM1-IB
        CALL PUTID(IDSTK(1,I+1),IDSTK(1,I))
        MSTK(I+1) = MSTK(I)
        NSTK(I+1) = NSTK(I)
        LSTK(I+1) = LSTK(I)+MNK
   24 CONTINUE
C 
C     DESTROY OLD VARIABLE
   25 BOT = BOT+1
C 
C     CREATE NEW VARIABLE
   30 IF (MN .EQ. 0) GO TO 99
      IF (BOT-2 .LE. TOP) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      K = BOT-1
      CALL PUTID(IDSTK(1,K), ID)
      IF (RHS .EQ. 1) GO TO 50
      IF (RHS .EQ. 2) GO TO 55
C 
C     STORE
   40 IF (K .LT. LSIZE) LSTK(K) = LSTK(K+1) - MN
      MSTK(K) = M
      NSTK(K) = N
      LK = LSTK(K)
      CALL WCOPY(MN,STKR(L),STKI(L),-1,STKR(LK),STKI(LK),-1)
      GO TO 90
C 
C     VECT(ARG)
   50 IF (MSTK(TOP-1) .LT. 0) GO TO 59
      MN1 = 1
      MN2 = 1
      L1 = 0
      L2 = 0
      IF (N.NE.1 .OR. NK.NE.1) GO TO 52
      L1 = LSTK(TOP-1)
      M1 = MSTK(TOP-1)
      MN1 = M1*NSTK(TOP-1)
      M2 = -1
      GO TO 60
   52 IF (M.NE.1 .OR. MK.NE.1) CALL ERROR(15)
      IF (ERR .GT. 0) RETURN
      L2 = LSTK(TOP-1)
      M2 = MSTK(TOP-1)
      MN2 = M2*NSTK(TOP-1)
      M1 = -1
      GO TO 60
C 
C     MATRIX(ARG,ARG)
   55 IF (MSTK(TOP-1).LT.0 .AND. MSTK(TOP-2).LT.0) GO TO 59
      L2 = LSTK(TOP-1)
      M2 = MSTK(TOP-1)
      MN2 = M2*NSTK(TOP-1)
      IF (M2 .LT. 0) MN2 = N
      L1 = LSTK(TOP-2)
      M1 = MSTK(TOP-2)
      MN1 = M1*NSTK(TOP-2)
      IF (M1 .LT. 0) MN1 = M
      GO TO 60
C 
   59 IF (MN .NE. MNK) CALL ERROR(15)
      IF (ERR .GT. 0) RETURN
      LK = LSTK(K)
      CALL WCOPY(MN,STKR(L),STKI(L),-1,STKR(LK),STKI(LK),-1)
      GO TO 90
C 
   60 IF (MN1.NE.M .OR. MN2.NE.N) CALL ERROR(15)
      IF (ERR .GT. 0) RETURN
      LL = 1
      IF (M1 .LT. 0) GO TO 62
      DO 61 I = 1, MN1
         LS = L1+I-1
         MK = MAX0(MK,IDINT(STKR(LS)))
         LL = MIN0(LL,IDINT(STKR(LS)))
   61 CONTINUE
   62 MK = MAX0(MK,M)
      IF (M2 .LT. 0) GO TO 64
      DO 63 I = 1, MN2
         LS = L2+I-1
         NK = MAX0(NK,IDINT(STKR(LS)))
         LL = MIN0(LL,IDINT(STKR(LS)))
   63 CONTINUE
   64 NK = MAX0(NK,N)
      IF (LL .LT. 1) CALL ERROR(21)
      IF (ERR .GT. 0) RETURN
      MNK = MK*NK
      LK = LSTK(K+1) - MNK
      ERR = LT + MT*NT - LK
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      LSTK(K) = LK
      MSTK(K) = MK
      NSTK(K) = NK
      CALL WSET(MNK,0.0D0,0.0D0,STKR(LK),STKI(LK),1)
      IF (NT .LT. 1) GO TO 67
      DO 66 J = 1, NT
         LS = LT+(J-1)*MT
         LL = LK+(J-1)*MK
         CALL WCOPY(MT,STKR(LS),STKI(LS),-1,STKR(LL),STKI(LL),-1)
   66 CONTINUE
   67 DO 68 J = 1, N
      DO 68 I = 1, M
        LI = L1+I-1
        IF (M1 .GT. 0) LI = L1 + IDINT(STKR(LI)) - 1
        LJ = L2+J-1
        IF (M2 .GT. 0) LJ = L2 + IDINT(STKR(LJ)) - 1
        LL = LK+LI-L1+(LJ-L2)*MK
        LS = L+I-1+(J-1)*M
        STKR(LL) = STKR(LS)
        STKI(LL) = STKI(LS)
   68 CONTINUE
      GO TO 90
C 
C     PRINT IF DESIRED AND POP STACK
   90 IF (SYM.NE.SEMI .AND. LCT(3).EQ.0) CALL PRINT(ID,K)
      IF (SYM.EQ.SEMI .AND. LCT(3).EQ.1) CALL PRINT(ID,K)
      IF (K .EQ. BOT-1) BOT = BOT-1
   99 IF (M .NE. 0) TOP = TOP - 1 - RHS
      IF (M .EQ. 0) TOP = TOP - 1
      RETURN
      END
