      SUBROUTINE MATFN4
C 
C     EVALUATE FUNCTIONS INVOLVING QR DECOMPOSITION (LEAST SQUARES)
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      DOUBLE PRECISION T,TOL,EPS,FLOP
      INTEGER QUOTE
      DATA QUOTE/49/
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN4',I4)
C 
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      IF (FIN .EQ. -1) GO TO 10
      IF (FIN .EQ. -2) GO TO 20
      GO TO 40
C 
C     RECTANGULAR MATRIX RIGHT DIVISION, A/A2
   10 L2 = LSTK(TOP+1)
      M2 = MSTK(TOP+1)
      N2 = NSTK(TOP+1)
      TOP = TOP + 1
      IF (N.GT.1 .AND. N.NE.N2) CALL ERROR(11)
      IF (ERR .GT. 0) RETURN
      CALL STACK1(QUOTE)
      IF (ERR .GT. 0) RETURN
      LL = L2+M2*N2
      CALL WCOPY(M*N,STKR(L),STKI(L),1,STKR(LL),STKI(LL),1)
      CALL WCOPY(M*N+M2*N2,STKR(L2),STKI(L2),1,STKR(L),STKI(L),1)
      LSTK(TOP) = L+M2*N2
      MSTK(TOP) = M
      NSTK(TOP) = N
      CALL STACK1(QUOTE)
      IF (ERR .GT. 0) RETURN
      TOP = TOP - 1
      M = N2
      N = M2
      GO TO 20
C 
C     RECTANGULAR MATRIX LEFT DIVISION A BACKSLASH A2
C 
   20 L2 = LSTK(TOP+1)
      M2 = MSTK(TOP+1)
      N2 = NSTK(TOP+1)
      IF (M2*N2 .GT. 1) GO TO 21
        M2 = M
        N2 = M
        ERR = L2+M*M - LSTK(BOT)
        IF (ERR .GT. 0) CALL ERROR(17)
        IF (ERR .GT. 0) RETURN
        CALL WSET(M*M-1,0.0D0,0.0D0,STKR(L2+1),STKI(L2+1),1)
        CALL WCOPY(M,STKR(L2),STKI(L2),0,STKR(L2),STKI(L2),M+1)
   21 IF (M2 .NE. M) CALL ERROR(12)
      IF (ERR .GT. 0) RETURN
      L3 = L2 + MAX0(M,N)*N2
      L4 = L3 + N
      ERR = L4 + N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      IF (M .GT. N) GO TO 23
      DO 22 JB = 1, N2
        J = N+1-JB
        LS = L2 + (J-1)*M
        LL = L2 + (J-1)*N
        CALL WCOPY(M,STKR(LS),STKI(LS),-1,STKR(LL),STKI(LL),-1)
   22 CONTINUE
   23 DO 24 J = 1, N
        BUF(J) = 0
   24 CONTINUE
      CALL WQRDC(STKR(L),STKI(L),M,M,N,STKR(L4),STKI(L4),
     $           BUF,STKR(L3),STKI(L3),1)
      K = 0
      EPS = STKR(VSIZE-4)
      T = DABS(STKR(L))+DABS(STKI(L))
      TOL = FLOP(DBLE(MAX0(M,N))*EPS*T)
      MN = MIN0(M,N)
      DO 27 J = 1, MN
        LS = L+J-1+(J-1)*M
        T = DABS(STKR(LS)) + DABS(STKI(LS))
        IF (T .GT. TOL) K = J
   27 CONTINUE
      IF (K .LT. MN) WRITE(WTE,28) K,TOL
      IF (K.LT.MN .AND. WIO.NE.0) WRITE(WIO,28) K,TOL
   28 FORMAT(1X,'RANK DEFICIENT,  RANK =',I4,',  TOL =',1PD13.4)
      MN = MAX0(M,N)
      DO 29 J = 1, N2
        LS = L2+(J-1)*MN
        CALL WQRSL(STKR(L),STKI(L),M,M,K,STKR(L4),STKI(L4),
     $             STKR(LS),STKI(LS),T,T,STKR(LS),STKI(LS),
     $             STKR(LS),STKI(LS),T,T,T,T,100,INFO)
        LL = LS+K
        CALL WSET(N-K,0.0D0,0.0D0,STKR(LL),STKI(LL),1)
   29 CONTINUE
      DO 31 J = 1, N
        BUF(J) = -BUF(J)
   31 CONTINUE
      DO 35 J = 1, N
        IF (BUF(J) .GT. 0) GO TO 35
        K = -BUF(J)
        BUF(J) = K
   33   CONTINUE
          IF (K .EQ. J) GO TO 34
          LS = L2+J-1
          LL = L2+K-1
          CALL WSWAP(N2,STKR(LS),STKI(LS),MN,STKR(LL),STKI(LL),MN)
          BUF(K) = -BUF(K)
          K = BUF(K)
          GO TO 33
   34   CONTINUE
   35 CONTINUE
      DO 36 J = 1, N2
        LS = L2+(J-1)*MN
        LL = L+(J-1)*N
        CALL WCOPY(N,STKR(LS),STKI(LS),1,STKR(LL),STKI(LL),1)
   36 CONTINUE
      MSTK(TOP) = N
      NSTK(TOP) = N2
      IF (FIN .EQ. -1) CALL STACK1(QUOTE)
      IF (ERR .GT. 0) RETURN
      GO TO 99
C 
C     QR
C 
   40 MM = MAX0(M,N)
      LS = L + MM*MM
      IF (LHS.EQ.1 .AND. FIN.EQ.1) LS = L
      LE = LS + M*N
      L4 = LE + MM
      ERR = L4+MM - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      IF (LS.NE.L) CALL WCOPY(M*N,STKR(L),STKI(L),1,STKR(LS),STKI(LS),1)
      JOB = 1
      IF (LHS.LT.3) JOB = 0
      DO 42 J = 1, N
        BUF(J) = 0
   42 CONTINUE
      CALL WQRDC(STKR(LS),STKI(LS),M,M,N,STKR(L4),STKI(L4),
     $            BUF,STKR(LE),STKI(LE),JOB)
      IF (LHS.EQ.1 .AND. FIN.EQ.1) GO TO 99
      CALL WSET(M*M,0.0D0,0.0D0,STKR(L),STKI(L),1)
      CALL WSET(M,1.0D0,0.0D0,STKR(L),STKI(L),M+1)
      DO 43 J = 1, M
        LL = L+(J-1)*M
        CALL WQRSL(STKR(LS),STKI(LS),M,M,N,STKR(L4),STKI(L4),
     $             STKR(LL),STKI(LL),STKR(LL),STKI(LL),T,T,
     $             T,T,T,T,T,T,10000,INFO)
   43 CONTINUE
      IF (FIN .EQ. 2) GO TO 99
      NSTK(TOP) = M
      DO 45 J = 1, N
        LL = LS+J+(J-1)*M
        CALL WSET(M-J,0.0D0,0.0D0,STKR(LL),STKI(LL),1)
   45 CONTINUE
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = LS
      MSTK(TOP) = M
      NSTK(TOP) = N
      IF (LHS .EQ. 2) GO TO 99
      CALL WSET(N*N,0.0D0,0.0D0,STKR(LE),STKI(LE),1)
      DO 47 J = 1, N
        LL = LE+BUF(J)-1+(J-1)*N
        STKR(LL) = 1.0D0
   47 CONTINUE
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = LE
      MSTK(TOP) = N
      NSTK(TOP) = N
      GO TO 99
C 
   99 RETURN
      END
