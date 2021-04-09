      SUBROUTINE MATFN6
C 
C     EVALUATE UTILITY FUNCTIONS
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      INTEGER SEMI,ID(4),UNIFOR(4),NORMAL(4),SEED(4)
      DOUBLE PRECISION EPS0,EPS,S,SR,SI,T
      DOUBLE PRECISION FLOP,URAND
      LOGICAL EQID
      DATA SEMI/39/
      DATA UNIFOR/30,23,18,15/,NORMAL/23,24,27,22/,SEED/28,14,14,13/
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN6',I4)
C     FUNCTIONS/FIN
C     MAGI DIAG SUM  PROD USER EYE  RAND ONES CHOP SIZE KRON  TRIL TRIU
C       1    2    3    4    5    6    7    8    9   10  11-13  14   15
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      GO TO (75,80,65,67,70,90,90,90,60,77,50,50,50,80,80),FIN
C 
C     KRONECKER PRODUCT
   50 IF (RHS .NE. 2) CALL ERROR(39)
      IF (ERR .GT. 0) RETURN
      TOP = TOP - 1
      L = LSTK(TOP)
      MA = MSTK(TOP)
      NA = NSTK(TOP)
      LA = L + MAX0(M*N*MA*NA,M*N+MA*NA)
      LB = LA + MA*NA
      ERR = LB + M*N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
C     MOVE A AND B ABOVE RESULT
      CALL WCOPY(MA*NA+M*N,STKR(L),STKI(L),1,STKR(LA),STKI(LA),1)
      DO 54 JA = 1, NA
        DO 53 J = 1, N
          LJ = LB + (J-1)*M
          DO 52 IA = 1, MA
C           GET J-TH COLUMN OF B
            CALL WCOPY(M,STKR(LJ),STKI(LJ),1,STKR(L),STKI(L),1)
C           ADDRESS OF A(IA,JA)
            LS = LA + IA-1 + (JA-1)*MA
            DO 51 I = 1, M
C             A(IA,JA) OP B(I,J)
              IF (FIN .EQ. 11) CALL WMUL(STKR(LS),STKI(LS),
     $           STKR(L),STKI(L),STKR(L),STKI(L))
              IF (FIN .EQ. 12) CALL WDIV(STKR(LS),STKI(LS),
     $           STKR(L),STKI(L),STKR(L),STKI(L))
              IF (FIN .EQ. 13) CALL WDIV(STKR(L),STKI(L),
     $           STKR(LS),STKI(LS),STKR(L),STKI(L))
              IF (ERR .GT. 0) RETURN
              L = L + 1
   51       CONTINUE
   52     CONTINUE
   53   CONTINUE
   54 CONTINUE
      MSTK(TOP) = M*MA
      NSTK(TOP) = N*NA
      GO TO 99
C 
C     CHOP
   60 EPS0 = 1.0D0
   61 EPS0 = EPS0/2.0D0
      T = FLOP(1.0D0 + EPS0)
      IF (T .GT. 1.0D0) GO TO 61
      EPS0 = 2.0D0*EPS0
      FLP(2) = IDINT(STKR(L))
      IF (SYM .NE. SEMI) WRITE(WTE,62) FLP(2)
   62 FORMAT(/1X,'CHOP ',I2,' PLACES.')
      EPS = 1.0D0
   63 EPS = EPS/2.0D0
      T = FLOP(1.0D0 + EPS)
      IF (T .GT. 1.0D0) GO TO 63
      EPS = 2.0D0*EPS
      T = STKR(VSIZE-4)
      IF (T.LT.EPS .OR. T.EQ.EPS0) STKR(VSIZE-4) = EPS
      MSTK(TOP) = 0
      GO TO 99
C 
C     SUM
   65 SR = 0.0D0
      SI = 0.0D0
      MN = M*N
      DO 66 I = 1, MN
         LS = L+I-1
         SR = FLOP(SR+STKR(LS))
         SI = FLOP(SI+STKI(LS))
   66 CONTINUE
      GO TO 69
C 
C     PROD
   67 SR = 1.0D0
      SI = 0.0D0
      MN = M*N
      DO 68 I = 1, MN
         LS = L+I-1
         CALL WMUL(STKR(LS),STKI(LS),SR,SI,SR,SI)
   68 CONTINUE
   69 STKR(L) = SR
      STKI(L) = SI
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
C 
C     USER
   70 S = 0.0D0
      T = 0.0D0
      IF (RHS .LT. 2) GO TO 72
      IF (RHS .LT. 3) GO TO 71
      T = STKR(L)
      TOP = TOP-1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
   71 S = STKR(L)
      TOP = TOP-1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
   72 CALL USER(STKR(L),M,N,S,T)
      CALL RSET(M*N,0.0D0,STKI(L),1)
      MSTK(TOP) = M
      NSTK(TOP) = N
      GO TO 99
C 
C     MAGIC
   75 N = MAX0(IDINT(STKR(L)),0)
      IF (N .EQ. 2) N = 0
      IF (N .GT. 0) CALL MAGIC(STKR(L),N,N)
      CALL RSET(N*N,0.0D0,STKI(L),1)
      MSTK(TOP) = N
      NSTK(TOP) = N
      GO TO 99
C 
C     SIZE
   77 STKR(L) = M
      STKR(L+1) = N
      STKI(L) = 0.0D0
      STKI(L+1) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 2
      IF (LHS .EQ. 1) GO TO 99
      NSTK(TOP) = 1
      TOP = TOP + 1
      LSTK(TOP) = L+1
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
C 
C     DIAG, TRIU, TRIL
   80 K = 0
      IF (RHS .NE. 2) GO TO 81
         K = IDINT(STKR(L))
         TOP = TOP-1
         L = LSTK(TOP)
         M = MSTK(TOP)
         N = NSTK(TOP)
   81 IF (FIN .GE. 14) GO TO 85
      IF (M .EQ. 1 .OR. N .EQ. 1) GO TO 83
      IF (K.GE.0) MN=MIN0(M,N-K)
      IF (K.LT.0) MN=MIN0(M+K,N)
      MSTK(TOP) = MAX0(MN,0)
      NSTK(TOP) = 1
      IF (MN .LE. 0) GO TO 99
      DO 82 I = 1, MN
         IF (K.GE.0) LS = L+(I-1)+(I+K-1)*M
         IF (K.LT.0) LS = L+(I-K-1)+(I-1)*M
         LL = L+I-1
         STKR(LL) = STKR(LS)
         STKI(LL) = STKI(LS)
   82 CONTINUE
      GO TO 99
   83 N = MAX0(M,N)+IABS(K)
      ERR = L+N*N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      MSTK(TOP) = N
      NSTK(TOP) = N
      DO 84 JB = 1, N
      DO 84 IB = 1, N
         J = N+1-JB
         I = N+1-IB
         SR = 0.0D0
         SI = 0.0D0
         IF (K.GE.0) LS = L+I-1
         IF (K.LT.0) LS = L+J-1
         LL = L+I-1+(J-1)*N
         IF (J-I .EQ. K) SR = STKR(LS)
         IF (J-I .EQ. K) SI = STKI(LS)
         STKR(LL) = SR
         STKI(LL) = SI
   84 CONTINUE
      GO TO 99
C 
C     TRIL, TRIU
   85 DO 87 J = 1, N
         LD = L + J - K - 1 + (J-1)*M
         IF (FIN .EQ. 14) LL = J - K - 1
         IF (FIN .EQ. 14) LS = LD - LL
         IF (FIN .EQ. 15) LL = M - J + K
         IF (FIN .EQ. 15) LS = LD + 1
         IF (LL .GT. 0) CALL WSET(LL,0.0D0,0.0D0,STKR(LS),STKI(LS),1)
   87 CONTINUE
      GO TO 99
C 
C     EYE, RAND, ONES
   90 IF (M.GT.1 .OR. RHS.EQ.0) GO TO 94
      IF (RHS .NE. 2) GO TO 91
        NN = IDINT(STKR(L))
        TOP = TOP-1
        L = LSTK(TOP)
        N = NSTK(TOP)
   91 IF (FIN.NE.7 .OR. N.LT.4) GO TO 93
      DO 92 I = 1, 4
        LS = L+I-1
        ID(I) = IDINT(STKR(LS))
   92 CONTINUE
      IF (EQID(ID,UNIFOR).OR.EQID(ID,NORMAL)) GO TO 97
      IF (EQID(ID,SEED)) GO TO 98
   93 IF (N .GT. 1) GO TO 94
      M = MAX0(IDINT(STKR(L)),0)
      IF (RHS .EQ. 2) N = MAX0(NN,0)
      IF (RHS .NE. 2) N = M
      ERR = L+M*N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      MSTK(TOP) = M
      NSTK(TOP) = N
      IF (M*N .EQ. 0) GO TO 99
   94 DO 96 J = 1, N
      DO 96 I = 1, M
        LL = L+I-1+(J-1)*M
        STKR(LL) = 0.0D0
        STKI(LL) = 0.0D0
        IF (I.EQ.J .OR. FIN.EQ.8) STKR(LL) = 1.0D0
        IF (FIN.EQ.7 .AND. RAN(2).EQ.0) STKR(LL) = FLOP(URAND(RAN(1)))
        IF (FIN.NE.7 .OR. RAN(2).EQ.0) GO TO 96
   95      SR = 2.0D0*URAND(RAN(1))-1.0D0
           SI = 2.0D0*URAND(RAN(1))-1.0D0
           T = SR*SR + SI*SI
           IF (T .GT. 1.0D0) GO TO 95
        STKR(LL) = FLOP(SR*DSQRT(-2.0D0*DLOG(T)/T))
   96 CONTINUE
      GO TO 99
C 
C     SWITCH UNIFORM AND NORMAL
   97 RAN(2) = ID(1) - UNIFOR(1)
      MSTK(TOP) = 0
      GO TO 99
C 
C     SEED
   98 IF (RHS .EQ. 2) RAN(1) = NN
      STKR(L) = RAN(1)
      MSTK(TOP) = 1
      IF (RHS .EQ. 2) MSTK(TOP) = 0
      NSTK(TOP) = 1
      GO TO 99
C 
   99 RETURN
      END
