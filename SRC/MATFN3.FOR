      SUBROUTINE MATFN3
C 
C     EVALUATE FUNCTIONS INVOLVING SINGULAR VALUE DECOMPOSITION
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      LOGICAL FRO,INF
      DOUBLE PRECISION P,S,T,TOL,EPS
      DOUBLE PRECISION WDOTCR,WDOTCI,PYTHAG,WNRM2,WASUM,FLOP
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN3',I4)
C 
      IF (FIN.EQ.1 .AND. RHS.EQ.2) TOP = TOP-1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      MN = M*N
      GO TO (50,70,10,30,70), FIN
C 
C     COND
C 
   10 LD = L + M*N
      L1 = LD + MIN0(M+1,N)
      L2 = L1 + N
      ERR = L2+MIN0(M,N) - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WSVDC(STKR(L),STKI(L),M,M,N,STKR(LD),STKI(LD),
     $           STKR(L1),STKI(L1),T,T,1,T,T,1,STKR(L2),STKI(L2),
     $           0,ERR)
      IF (ERR .NE. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
      S = STKR(LD)
      LD = LD + MIN0(M,N) - 1
      T = STKR(LD)
      IF (T .EQ. 0.0D0) GO TO 13
      STKR(L) = FLOP(S/T)
      STKI(L) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
   13 WRITE(WTE,14)
      IF (WIO .NE. 0) WRITE(WIO,14)
   14 FORMAT(1X,'CONDITION IS INFINITE')
      MSTK(TOP) = 0
      GO TO 99
C 
C     NORM
C 
   30 P = 2.0D0
      INF = .FALSE.
      IF (RHS .NE. 2) GO TO 31
      FRO = IDINT(STKR(L)).EQ.15 .AND. MN.GT.1
      INF = IDINT(STKR(L)).EQ.18 .AND. MN.GT.1
      IF (.NOT. FRO) P = STKR(L)
      TOP = TOP-1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      MN = M*N
      IF (FRO) M = MN
      IF (FRO) N = 1
   31 IF (M .GT. 1 .AND. N .GT. 1) GO TO 40
      IF (P .EQ. 1.0D0) GO TO 36
      IF (P .EQ. 2.0D0) GO TO 38
      I = IWAMAX(MN,STKR(L),STKI(L),1) + L - 1
      S = DABS(STKR(I)) + DABS(STKI(I))
      IF (INF .OR. S .EQ. 0.0D0) GO TO 49
      T = 0.0D0
      DO 33 I = 1, MN
         LS = L+I-1
         T = FLOP(T + (PYTHAG(STKR(LS),STKI(LS))/S)**P)
   33 CONTINUE
      IF (P .NE. 0.0D0) P = 1.0D0/P
      S = FLOP(S*T**P)
      GO TO 49
   36 S = WASUM(MN,STKR(L),STKI(L),1)
      GO TO 49
   38 S = WNRM2(MN,STKR(L),STKI(L),1)
      GO TO 49
C 
C     MATRIX NORM
C 
   40 IF (INF) GO TO 43
      IF (P .EQ. 1.0D0) GO TO 46
      IF (P .NE. 2.0D0) CALL ERROR(23)
      IF (ERR .GT. 0) RETURN
      LD = L + M*N
      L1 = LD + MIN0(M+1,N)
      L2 = L1 + N
      ERR = L2+MIN0(M,N) - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WSVDC(STKR(L),STKI(L),M,M,N,STKR(LD),STKI(LD),
     $           STKR(L1),STKI(L1),T,T,1,T,T,1,STKR(L2),STKI(L2),
     $           0,ERR)
      IF (ERR .NE. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
      S = STKR(LD)
      GO TO 49
   43 S = 0.0D0
      DO 45 I = 1, M
         LI = L+I-1
         T = WASUM(N,STKR(LI),STKI(LI),M)
         S = DMAX1(S,T)
   45 CONTINUE
      GO TO 49
   46 S = 0.0D0
      DO 48 J = 1, N
         LJ = L+(J-1)*M
         T = WASUM(M,STKR(LJ),STKI(LJ),1)
         S = DMAX1(S,T)
   48 CONTINUE
      GO TO 49
   49 STKR(L) = S
      STKI(L) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
C 
C     SVD
C 
   50 IF (LHS .NE. 3) GO TO 52
      K = M
      IF (RHS .EQ. 2) K = MIN0(M,N)
      LU = L + M*N
      LD = LU + M*K
      LV = LD + K*N
      L1 = LV + N*N
      L2 = L1 + N
      ERR = L2+MIN0(M,N) - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      JOB = 11
      IF (RHS .EQ. 2) JOB = 21
      CALL WSVDC(STKR(L),STKI(L),M,M,N,STKR(LD),STKI(LD),
     $        STKR(L1),STKI(L1),STKR(LU),STKI(LU),M,STKR(LV),STKI(LV),
     $        N,STKR(L2),STKI(L2),JOB,ERR)
      DO 51 JB = 1, N
      DO 51 I = 1, K
        J = N+1-JB
        LL = LD+I-1+(J-1)*K
        IF (I.NE.J) STKR(LL) = 0.0D0
        STKI(LL) = 0.0D0
        LS = LD+I-1
        IF (I.EQ.J) STKR(LL) = STKR(LS)
        LS = L1+I-1
        IF (ERR.NE.0 .AND. I.EQ.J-1) STKR(LL) = STKR(LS)
   51 CONTINUE
      IF (ERR .NE. 0) CALL ERROR(24)
      ERR = 0
      CALL WCOPY(M*K+K*N+N*N,STKR(LU),STKI(LU),1,STKR(L),STKI(L),1)
      MSTK(TOP) = M
      NSTK(TOP) = K
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = L + M*K
      MSTK(TOP) = K
      NSTK(TOP) = N
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = L + M*K + K*N
      MSTK(TOP) = N
      NSTK(TOP) = N
      GO TO 99
C 
   52 LD = L + M*N
      L1 = LD + MIN0(M+1,N)
      L2 = L1 + N
      ERR = L2+MIN0(M,N) - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WSVDC(STKR(L),STKI(L),M,M,N,STKR(LD),STKI(LD),
     $           STKR(L1),STKI(L1),T,T,1,T,T,1,STKR(L2),STKI(L2),
     $           0,ERR)
      IF (ERR .NE. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
      K = MIN0(M,N)
      CALL WCOPY(K,STKR(LD),STKI(LD),1,STKR(L),STKI(L),1)
      MSTK(TOP) = K
      NSTK(TOP) = 1
      GO TO 99
C 
C     PINV AND RANK
C 
   70 TOL = -1.0D0
      IF (RHS .NE. 2) GO TO 71
      TOL = STKR(L)
      TOP = TOP-1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
   71 LU = L + M*N
      LD = LU + M*M
      IF (FIN .EQ. 5) LD = L + M*N
      LV = LD + M*N
      L1 = LV + N*N
      IF (FIN .EQ. 5) L1 = LD + N
      L2 = L1 + N
      ERR = L2+MIN0(M,N) - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      IF (FIN .EQ. 2) JOB = 11
      IF (FIN .EQ. 5) JOB = 0
      CALL WSVDC(STKR(L),STKI(L),M,M,N,STKR(LD),STKI(LD),
     $        STKR(L1),STKI(L1),STKR(LU),STKI(LU),M,STKR(LV),STKI(LV),
     $        N,STKR(L2),STKI(L2),JOB,ERR)
      IF (ERR .NE. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
      EPS = STKR(VSIZE-4)
      IF (TOL .LT. 0.0D0) TOL = FLOP(DBLE(MAX0(M,N))*EPS*STKR(LD))
      MN = MIN0(M,N)
      K = 0
      DO 72 J = 1, MN
        LS = LD+J-1
        S = STKR(LS)
        IF (S .LE. TOL) GO TO 73
        K = J
        LL = LV+(J-1)*N
        IF (FIN .EQ. 2) CALL WRSCAL(N,1.0D0/S,STKR(LL),STKI(LL),1)
   72 CONTINUE
   73 IF (FIN .EQ. 5) GO TO 78
      DO 76 J = 1, M
      DO 76 I = 1, N
        LL = L+I-1+(J-1)*N
        L1 = LV+I-1
        L2 = LU+J-1
        STKR(LL) = WDOTCR(K,STKR(L2),STKI(L2),M,STKR(L1),STKI(L1),N)
        STKI(LL) = WDOTCI(K,STKR(L2),STKI(L2),M,STKR(L1),STKI(L1),N)
   76 CONTINUE
      MSTK(TOP) = N
      NSTK(TOP) = M
      GO TO 99
   78 STKR(L) = DBLE(K)
      STKI(L) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
C 
   99 RETURN
      END
