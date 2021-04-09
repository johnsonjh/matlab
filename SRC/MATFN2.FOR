      SUBROUTINE MATFN2
C 
C     EVALUATE ELEMENTARY FUNCTIONS AND FUNCTIONS INVOLVING
C     EIGENVALUES AND EIGENVECTORS
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      DOUBLE PRECISION PYTHAG,ROUNDM,TR,TI,SR,SI,POWR,POWI,FLOP
      LOGICAL HERM,SCHUR,VECT,HESS
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN2',I4)
C 
C     FUNCTIONS/FIN
C     **   SIN  COS ATAN  EXP  SQRT LOG
C      0    1    2    3    4    5    6
C    EIG  SCHU HESS POLY ROOT
C     11   12   13   14   15
C    ABS  ROUN REAL IMAG CONJ
C     21   22   23   24   25
      IF (FIN .NE. 0) GO TO 05
         L = LSTK(TOP+1)
         POWR = STKR(L)
         POWI = STKI(L)
   05 L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      IF (FIN .GE. 11 .AND. FIN .LE. 13) GO TO 10
      IF (FIN .EQ. 14 .AND. (M.EQ.1 .OR. N.EQ.1)) GO TO 50
      IF (FIN .EQ. 14) GO TO 10
      IF (FIN .EQ. 15) GO TO 60
      IF (FIN .GT. 20) GO TO 40
      IF (M .EQ. 1 .OR. N .EQ. 1) GO TO 40
C 
C     EIGENVALUES AND VECTORS
   10 IF (M .NE. N) CALL ERROR(20)
      IF (ERR .GT. 0) RETURN
      SCHUR = FIN .EQ. 12
      HESS = FIN .EQ. 13
      VECT = LHS.EQ.2 .OR. FIN.LT.10
      NN = N*N
      L2 = L + NN
      LD = L2 + NN
      LE = LD + N
      LW = LE + N
      ERR = LW+N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WCOPY(NN,STKR(L),STKI(L),1,STKR(L2),STKI(L2),1)
C 
C     CHECK IF HERMITIAN
      DO 15 J = 1, N
      DO 15 I = 1, J
         LS = L+I-1+(J-1)*N
         LL = L+(I-1)*N+J-1
         HERM = STKR(LL).EQ.STKR(LS) .AND. STKI(LL).EQ.-STKI(LS)
         IF (.NOT. HERM) GO TO 30
   15 CONTINUE
C 
C     HERMITIAN EIGENVALUE PROBLEM
      CALL WSET(NN,0.0D0,0.0D0,STKR(L),STKI(L),1)
      CALL WSET(N,1.0D0,0.0D0,STKR(L),STKI(L),N+1)
      CALL WSET(N,0.0D0,0.0D0,STKI(LD),STKI(LE),1)
      JOB = 0
      IF (VECT) JOB = 1
      CALL HTRIDI(N,N,STKR(L2),STKI(L2),STKR(LD),STKR(LE),
     $            STKR(LE),STKR(LW))
      IF (.NOT.HESS) CALL IMTQL2(N,N,STKR(LD),STKR(LE),STKR(L),ERR,JOB)
      IF (ERR .GT. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
      IF (JOB .NE. 0)
     $  CALL HTRIBK(N,N,STKR(L2),STKI(L2),STKR(LW),N,STKR(L),STKI(L))
      GO TO 31
C 
C     NON-HERMITIAN EIGENVALUE PROBLEM
   30 CALL CORTH(N,N,1,N,STKR(L2),STKI(L2),STKR(LW),STKI(LW))
      IF (.NOT.VECT .AND. HESS) GO TO 31
      JOB = 0
      IF (VECT) JOB = 2
      IF (VECT .AND. SCHUR) JOB = 1
      IF (HESS) JOB = 3
      CALL COMQR3(N,N,1,N,STKR(LW),STKI(LW),STKR(L2),STKI(L2),
     $            STKR(LD),STKI(LD),STKR(L),STKI(L),ERR,JOB)
      IF (ERR .GT. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
C 
C     VECTORS
   31 IF (.NOT.VECT) GO TO 34
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = L2
      MSTK(TOP) = N
      NSTK(TOP) = N
C 
C     DIAGONAL OF VALUES OR CANONICAL FORMS
   34 IF (.NOT.VECT .AND. .NOT.SCHUR .AND. .NOT.HESS) GO TO 37
      DO 36 J = 1, N
         LJ = L2+(J-1)*N
         IF (SCHUR .AND. (.NOT.HERM)) LJ = LJ+J
         IF (HESS .AND. (.NOT.HERM)) LJ = LJ+J+1
         LL = L2+J*N-LJ
         CALL WSET(LL,0.0D0,0.0D0,STKR(LJ),STKI(LJ),1)
   36 CONTINUE
      IF (.NOT.HESS .OR. HERM)
     $   CALL WCOPY(N,STKR(LD),STKI(LD),1,STKR(L2),STKI(L2),N+1)
      LL = L2+1
      IF (HESS .AND. HERM)
     $   CALL WCOPY(N-1,STKR(LE+1),STKI(LE+1),1,STKR(LL),STKI(LL),N+1)
      LL = L2+N
      IF (HESS .AND. HERM)
     $   CALL WCOPY(N-1,STKR(LE+1),STKI(LE+1),1,STKR(LL),STKI(LL),N+1)
      IF (FIN .LT. 10) GO TO 42
      IF (VECT .OR. .NOT.(SCHUR.OR.HESS)) GO TO 99
      CALL WCOPY(NN,STKR(L2),STKI(L2),1,STKR(L),STKI(L),1)
      GO TO 99
C 
C     VECTOR OF EIGENVALUES
   37 IF (FIN .EQ. 14) GO TO 52
      CALL WCOPY(N,STKR(LD),STKI(LD),1,STKR(L),STKI(L),1)
      NSTK(TOP) = 1
      GO TO 99
C 
C     ELEMENTARY FUNCTIONS
C     FOR MATRICES.. X,D = EIG(A), FUN(A) = X*FUN(D)/X
   40 INC = 1
      N = M*N
      L2 = L
      GO TO 44
   42 INC = N+1
   44 DO 46 J = 1, N
        LS = L2+(J-1)*INC
        SR = STKR(LS)
        SI = STKI(LS)
        TI = 0.0D0
        IF (FIN .NE. 0) GO TO 45
          CALL WLOG(SR,SI,SR,SI)
          CALL WMUL(SR,SI,POWR,POWI,SR,SI)
          TR = DEXP(SR)*DCOS(SI)
          TI = DEXP(SR)*DSIN(SI)
   45   IF (FIN .EQ. 1) TR = DSIN(SR)*DCOSH(SI)
        IF (FIN .EQ. 1) TI = DCOS(SR)*DSINH(SI)
        IF (FIN .EQ. 2) TR = DCOS(SR)*DCOSH(SI)
        IF (FIN .EQ. 2) TI = -DSIN(SR)*DSINH(SI)
        IF (FIN .EQ. 3) CALL WATAN(SR,SI,TR,TI)
        IF (FIN .EQ. 4) TR = DEXP(SR)*DCOS(SI)
        IF (FIN .EQ. 4) TI = DEXP(SR)*DSIN(SI)
        IF (FIN .EQ. 5) CALL WSQRT(SR,SI,TR,TI)
        IF (FIN .EQ. 6) CALL WLOG(SR,SI,TR,TI)
        IF (FIN .EQ. 21) TR = PYTHAG(SR,SI)
        IF (FIN .EQ. 22) TR = ROUNDM(SR)
        IF (FIN .EQ. 23) TR = SR
        IF (FIN .EQ. 24) TR = SI
        IF (FIN .EQ. 25) TR = SR
        IF (FIN .EQ. 25) TI = -SI
        IF (ERR .GT. 0) RETURN
        STKR(LS) = FLOP(TR)
        STKI(LS) = 0.0D0
        IF (TI .NE. 0.0D0) STKI(LS) = FLOP(TI)
   46 CONTINUE
      IF (INC .EQ. 1) GO TO 99
      DO 48 J = 1, N
        LS = L2+(J-1)*INC
        SR = STKR(LS)
        SI = STKI(LS)
        LS = L+(J-1)*N
        LL = L2+(J-1)*N
        CALL WCOPY(N,STKR(LS),STKI(LS),1,STKR(LL),STKI(LL),1)
        CALL WSCAL(N,SR,SI,STKR(LS),STKI(LS),1)
   48 CONTINUE
C     SIGNAL MATFN1 TO DIVIDE BY EIGENVECTORS
      FUN = 21
      FIN = -1
      TOP = TOP-1
      GO TO 99
C 
C     POLY
C     FORM POLYNOMIAL WITH GIVEN VECTOR AS ROOTS
   50 N = MAX0(M,N)
      LD = L+N+1
      CALL WCOPY(N,STKR(L),STKI(L),1,STKR(LD),STKI(LD),1)
C 
C     FORM CHARACTERISTIC POLYNOMIAL
   52 CALL WSET(N+1,0.0D0,0.0D0,STKR(L),STKI(L),1)
      STKR(L) = 1.0D0
      DO 56 J = 1, N
         CALL WAXPY(J,-STKR(LD),-STKI(LD),STKR(L),STKI(L),-1,
     $              STKR(L+1),STKI(L+1),-1)
         LD = LD+1
   56 CONTINUE
      MSTK(TOP) = N+1
      NSTK(TOP) = 1
      GO TO 99
C 
C     ROOTS
   60 LL = L+M*N
      STKR(LL) = -1.0D0
      STKI(LL) = 0.0D0
      K = -1
   61 K = K+1
      L1 = L+K
      IF (DABS(STKR(L1))+DABS(STKI(L1)) .EQ. 0.0D0) GO TO 61
      N = MAX0(M*N - K-1, 0)
      IF (N .LE. 0) GO TO 65
      L2 = L1+N+1
      LW = L2+N*N
      ERR = LW+N - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WSET(N*N+N,0.0D0,0.0D0,STKR(L2),STKI(L2),1)
      DO 64 J = 1, N
         LL = L2+J+(J-1)*N
         STKR(LL) = 1.0D0
         LS = L1+J
         LL = L2+(J-1)*N
         CALL WDIV(-STKR(LS),-STKI(LS),STKR(L1),STKI(L1),
     $             STKR(LL),STKI(LL))
         IF (ERR .GT. 0) RETURN
   64 CONTINUE
      CALL COMQR3(N,N,1,N,STKR(LW),STKI(LW),STKR(L2),STKI(L2),
     $            STKR(L),STKI(L),TR,TI,ERR,0)
      IF (ERR .GT. 0) CALL ERROR(24)
      IF (ERR .GT. 0) RETURN
   65 MSTK(TOP) = N
      NSTK(TOP) = 1
      GO TO 99
   99 RETURN
      END
