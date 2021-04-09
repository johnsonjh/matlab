      SUBROUTINE MATFN1
C 
C     EVALUATE FUNCTIONS INVOLVING GAUSSIAN ELIMINATION
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      DOUBLE PRECISION DTR(2),DTI(2),SR,SI,RCOND,T,T0,T1,FLOP,EPS,WASUM
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN1',I4)
C 
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      IF (FIN .EQ. -1) GO TO 10
      IF (FIN .EQ. -2) GO TO 20
      GO TO (30,40,50,60,70,80,85),FIN
C 
C     MATRIX RIGHT DIVISION, A/A2
   10 L2 = LSTK(TOP+1)
      M2 = MSTK(TOP+1)
      N2 = NSTK(TOP+1)
      IF (M2 .NE. N2) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      IF (M*N .EQ. 1) GO TO 16
      IF (N .NE. N2) CALL ERROR(11)
      IF (ERR .GT. 0) RETURN
      L3 = L2 + M2*N2
      ERR = L3+N2 - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WGECO(STKR(L2),STKI(L2),M2,N2,BUF,RCOND,STKR(L3),STKI(L3))
      IF (RCOND .EQ. 0.0D0) CALL ERROR(19)
      IF (ERR .GT. 0) RETURN
      T = FLOP(1.0D0 + RCOND)
      IF (T.EQ.1.0D0 .AND. FUN.NE.21) WRITE(WTE,11) RCOND
      IF (T.EQ.1.0D0 .AND. FUN.NE.21 .AND. WIO.NE.0) WRITE(WIO,11) RCOND
   11 FORMAT(1X,'WARNING.'
     $  /1X,'MATRIX IS CLOSE TO SINGULAR OR BADLY SCALED.'
     $  /1X,'RESULTS MAY BE INACCURATE.  RCOND =', 1PD13.4/)
      IF (T.EQ.1.0D0 .AND. FUN.EQ.21) WRITE(WTE,12) RCOND
      IF (T.EQ.1.0D0 .AND. FUN.EQ.21 .AND. WIO.NE.0) WRITE(WIO,12) RCOND
   12 FORMAT(1X,'WARNING.'
     $  /1X,'EIGENVECTORS ARE BADLY CONDITIONED.'
     $  /1X,'RESULTS MAY BE INACCURATE.  RCOND =', 1PD13.4/)
      DO 15 I = 1, M
         DO 13 J = 1, N
            LS = L+I-1+(J-1)*M
            LL = L3+J-1
            STKR(LL) = STKR(LS)
            STKI(LL) = -STKI(LS)
   13    CONTINUE
         CALL WGESL(STKR(L2),STKI(L2),M2,N2,BUF,STKR(L3),STKI(L3),1)
         DO 14 J = 1, N
            LL = L+I-1+(J-1)*M
            LS = L3+J-1
            STKR(LL) = STKR(LS)
            STKI(LL) = -STKI(LS)
   14    CONTINUE
   15 CONTINUE
      IF (FUN .NE. 21) GO TO 99
C 
C     CHECK FOR IMAGINARY ROUNDOFF IN MATRIX FUNCTIONS
      SR = WASUM(N*N,STKR(L),STKR(L),1)
      SI = WASUM(N*N,STKI(L),STKI(L),1)
      EPS = STKR(VSIZE-4)
      T = EPS*SR
      IF (DDT .EQ. 18) WRITE(WTE,115) SR,SI,EPS,T
  115 FORMAT(1X,'SR,SI,EPS,T',1P4D13.4)
      IF (SI .LE. EPS*SR) CALL RSET(N*N,0.0D0,STKI(L),1)
      GO TO 99
C 
   16 SR = STKR(L)
      SI = STKI(L)
      N = N2
      M = N
      MSTK(TOP) = N
      NSTK(TOP) = N
      CALL WCOPY(N*N,STKR(L2),STKI(L2),1,STKR(L),STKI(L),1)
      GO TO 30
C 
C     MATRIX LEFT DIVISION A BACKSLASH A2
   20 L2 = LSTK(TOP+1)
      M2 = MSTK(TOP+1)
      N2 = NSTK(TOP+1)
      IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      IF (M2*N2 .EQ. 1) GO TO 26
      L3 = L2 + M2*N2
      ERR = L3+N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WGECO(STKR(L),STKI(L),M,N,BUF,RCOND,STKR(L3),STKI(L3))
      IF (RCOND .EQ. 0.0D0) CALL ERROR(19)
      IF (ERR .GT. 0) RETURN
      T = FLOP(1.0D0 + RCOND)
      IF (T .EQ. 1.0D0) WRITE(WTE,11) RCOND
      IF (T.EQ.1.0D0 .AND. WIO.NE.0) WRITE(WIO,11) RCOND
      IF (M2 .NE. N) CALL ERROR(12)
      IF (ERR .GT. 0) RETURN
      DO 23 J = 1, N2
         LJ = L2+(J-1)*M2
         CALL WGESL(STKR(L),STKI(L),M,N,BUF,STKR(LJ),STKI(LJ),0)
   23 CONTINUE
      NSTK(TOP) = N2
      CALL WCOPY(M2*N2,STKR(L2),STKI(L2),1,STKR(L),STKI(L),1)
      GO TO 99
   26 SR = STKR(L2)
      SI = STKI(L2)
      GO TO 30
C 
C     INV
C 
   30 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      IF (DDT .EQ. 17) GO TO 32
      DO 31 J = 1, N
      DO 31 I = 1, N
        LS = L+I-1+(J-1)*N
        T0 = STKR(LS)
        T1 = FLOP(1.0D0/(DBLE(I+J-1)))
        IF (T0 .NE. T1) GO TO 32
   31 CONTINUE
      GO TO 72
   32 L3 = L + N*N
      ERR = L3+N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WGECO(STKR(L),STKI(L),M,N,BUF,RCOND,STKR(L3),STKI(L3))
      IF (RCOND .EQ. 0.0D0) CALL ERROR(19)
      IF (ERR .GT. 0) RETURN
      T = FLOP(1.0D0 + RCOND)
      IF (T .EQ. 1.0D0) WRITE(WTE,11) RCOND
      IF (T.EQ.1.0D0 .AND. WIO.NE.0) WRITE(WIO,11) RCOND
      CALL WGEDI(STKR(L),STKI(L),M,N,BUF,DTR,DTI,STKR(L3),STKI(L3),1)
      IF (FIN .LT. 0) CALL WSCAL(N*N,SR,SI,STKR(L),STKI(L),1)
      GO TO 99
C 
C     DET
C 
   40 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      CALL WGEFA(STKR(L),STKI(L),M,N,BUF,INFO)
      CALL WGEDI(STKR(L),STKI(L),M,N,BUF,DTR,DTI,SR,SI,10)
      K = IDINT(DTR(2))
      KA = IABS(K)+2
      T = 1.0D0
      DO 41 I = 1, KA
         T = T/10.0D0
         IF (T .EQ. 0.0D0) GO TO 42
   41 CONTINUE
      STKR(L) = DTR(1)*10.D0**K
      STKI(L) = DTI(1)*10.D0**K
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      GO TO 99
   42 IF (DTI(1) .EQ. 0.0D0) WRITE(WTE,43) DTR(1),K
      IF (DTI(1) .NE. 0.0D0) WRITE(WTE,44) DTR(1),DTI(1),K
   43 FORMAT(1X,'DET =  ',F7.4,7H * 10**,I4)
   44 FORMAT(1X,'DET =  ',F7.4,' + ',F7.4,' i ',7H * 10**,I4)
      STKR(L) = DTR(1)
      STKI(L) = DTI(1)
      STKR(L+1) = DTR(2)
      STKI(L+1) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 2
      GO TO 99
C 
C     RCOND
C 
   50 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      L3 = L + N*N
      ERR = L3+N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WGECO(STKR(L),STKI(L),M,N,BUF,RCOND,STKR(L3),STKI(L3))
      STKR(L) = RCOND
      STKI(L) = 0.0D0
      MSTK(TOP) = 1
      NSTK(TOP) = 1
      IF (LHS .EQ. 1) GO TO 99
      L = L + 1
      CALL WCOPY(N,STKR(L3),STKI(L3),1,STKR(L),STKI(L),1)
      TOP = TOP + 1
      LSTK(TOP) = L
      MSTK(TOP) = N
      NSTK(TOP) = 1
      GO TO 99
C 
C     LU
C 
   60 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      CALL WGEFA(STKR(L),STKI(L),M,N,BUF,INFO)
      IF (LHS .NE. 2) GO TO 99
      NN = N*N
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = L + NN
      MSTK(TOP) = N
      NSTK(TOP) = N
      ERR = L+NN+NN - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      DO 64 KB = 1, N
        K = N+1-KB
        DO 61 I = 1, N
          LL = L+I-1+(K-1)*N
          LU = LL + NN
          IF (I .LE. K) STKR(LU) = STKR(LL)
          IF (I .LE. K) STKI(LU) = STKI(LL)
          IF (I .GT. K) STKR(LU) = 0.0D0
          IF (I .GT. K) STKI(LU) = 0.0D0
          IF (I .LT. K) STKR(LL) = 0.0D0
          IF (I .LT. K) STKI(LL) = 0.0D0
          IF (I .EQ. K) STKR(LL) = 1.0D0
          IF (I .EQ. K) STKI(LL) = 0.0D0
          IF (I .GT. K) STKR(LL) = -STKR(LL)
          IF (I .GT. K) STKI(LL) = -STKI(LL)
   61   CONTINUE
        I = BUF(K)
        IF (I .EQ. K) GO TO 64
        LI = L+I-1+(K-1)*N
        LK = L+K-1+(K-1)*N
        CALL WSWAP(N-K+1,STKR(LI),STKI(LI),N,STKR(LK),STKI(LK),N)
   64 CONTINUE
      GO TO 99
C 
C     HILBERT
   70 N = IDINT(STKR(L))
      MSTK(TOP) = N
      NSTK(TOP) = N
   72 CALL HILBER(STKR(L),N,N)
      CALL RSET(N*N,0.0D0,STKI(L),1)
      IF (FIN .LT. 0) CALL WSCAL(N*N,SR,SI,STKR(L),STKI(L),1)
      GO TO 99
C 
C     CHOLESKY
   80 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      CALL WPOFA(STKR(L),STKI(L),M,N,ERR)
      IF (ERR .NE. 0) CALL ERROR(29)
      IF (ERR .GT. 0) RETURN
      DO 81 J = 1, N
        LL = L+J+(J-1)*M
        CALL WSET(M-J,0.0D0,0.0D0,STKR(LL),STKI(LL),1)
   81 CONTINUE
      GO TO 99
C 
C     RREF
   85 IF (RHS .LT. 2) GO TO 86
        TOP = TOP-1
        L = LSTK(TOP)
        IF (MSTK(TOP) .NE. M) CALL ERROR(5)
        IF (ERR .GT. 0) RETURN
        N = N + NSTK(TOP)
   86 CALL RREF(STKR(L),STKI(L),M,M,N,STKR(VSIZE-4))
      NSTK(TOP) = N
      GO TO 99
C 
   99 RETURN
      END
