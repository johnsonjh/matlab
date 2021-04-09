      SUBROUTINE MATFN5
C 
C     FILE HANDLING AND OTHER I/O
C 
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER ALFA(52),ALFB(52),ALFL,CASE
      INTEGER IDS(4,32),PSTK(32),RSTK(32),PSIZE,PT,PTZ
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /ALFS/ ALFA,ALFB,ALFL,CASE
      COMMON /RECU/ IDS,PSTK,RSTK,PSIZE,PT,PTZ
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      INTEGER EOL,CH,BLANK,FLAG,TOP2,PLUS,MINUS,QUOTE,SEMI,LRAT,MRAT
      INTEGER ID(4)
      DOUBLE PRECISION EPS,B,S,T,FLOP,WASUM
      LOGICAL TEXT
      DATA EOL/99/,BLANK/36/,PLUS/41/,MINUS/42/,QUOTE/49/,SEMI/39/
      DATA LRAT/5/,MRAT/100/
C 
      IF (DDT .EQ. 1) WRITE(WTE,100) FIN
  100 FORMAT(1X,'MATFN5',I4)
C     FUNCTIONS/FIN
C     EXEC SAVE LOAD PRIN DIAR DISP BASE LINE CHAR PLOT RAT  DEBU
C      1    2    3    4    5    6    7    8    9   10   11   12
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      IF (FIN .GT. 5) GO TO 15
C 
C     CONVERT FILE NAME
      MN = M*N
      FLAG = 3
      IF (SYM .EQ. SEMI) FLAG = 0
      IF (RHS .LT. 2) GO TO 12
         FLAG = IDINT(STKR(L))
         TOP2 = TOP
         TOP = TOP-1
         L = LSTK(TOP)
         MN = MSTK(TOP)*NSTK(TOP)
   12 LUN = -1
      IF (MN.EQ.1 .AND. STKR(L).LT.10.0D0) LUN = IDINT(STKR(L))
      IF (LUN .GE. 0) GO TO 15
      DO 14 J = 1, 32
         LS = L+J-1
         IF (J .LE. MN) CH = IDINT(STKR(LS))
         IF (J .GT. MN) CH = BLANK
         IF (CH.LT.0 .OR. CH.GE.ALFL) CALL ERROR(38)
         IF (ERR .GT. 0) RETURN
         IF (CASE .EQ. 0) BUF(J) = ALFA(CH+1)
         IF (CASE .EQ. 1) BUF(J) = ALFB(CH+1)
   14 CONTINUE
C 
   15 GO TO (20,30,35,25,27,60,65,70,50,80,40,95),FIN
C 
C     EXEC
   20 IF (LUN .EQ. 0) GO TO 23
      K = LPT(6)
      LIN(K+1) = LPT(1)
      LIN(K+2) = LPT(3)
      LIN(K+3) = LPT(6)
      LIN(K+4) = PTZ
      LIN(K+5) = RIO
      LIN(K+6) = LCT(4)
      LPT(1) = K + 7
      LCT(4) = FLAG
      PTZ = PT - 4
      IF (RIO .EQ. RTE) RIO = 10
      RIO = RIO + 1
      IF (LUN .GT. 0) RIO = LUN
      IF (LUN .LT. 0) CALL FILES(RIO,BUF)
      IF (FLAG .GE. 4) WRITE(WTE,22)
   22 FORMAT(1X,'PAUSE MODE. ENTER BLANK LINES.')
      SYM = EOL
      MSTK(TOP) = 0
      GO TO 99
C 
C     EXEC(0)
   23 RIO = RTE
      ERR = 99
      GO TO 99
C 
C     PRINT
   25 K = WTE
      WTE = LUN
      IF (LUN .LT. 0) WTE = 7
      IF (LUN .LT. 0) CALL FILES(WTE,BUF)
      L = LCT(2)
      LCT(2) = 9999
      IF (RHS .GT. 1) CALL PRINT(SYN,TOP2)
      LCT(2) = L
      WTE = K
      MSTK(TOP) = 0
      GO TO 99
C 
C     DIARY
   27 WIO = LUN
      IF (LUN .LT. 0) WIO = 8
      IF (LUN .LT. 0) CALL FILES(WIO,BUF)
      MSTK(TOP) = 0
      GO TO 99
C 
C     SAVE
   30 IF (LUN .LT. 0) LUNIT = 1
      IF (LUN .LT. 0) CALL FILES(LUNIT,BUF)
      IF (LUN .GT. 0) LUNIT = LUN
      K = LSIZE-4
      IF (K .LT. BOT) K = LSIZE
      IF (RHS .EQ. 2) K = TOP2
      IF (RHS .EQ. 2) CALL PUTID(IDSTK(1,K),SYN)
   32 L = LSTK(K)
      M = MSTK(K)
      N = NSTK(K)
      DO 34 I = 1, 4
         J = IDSTK(I,K)+1
         BUF(I) = ALFA(J)
   34 CONTINUE
      IMG = 0
      IF (WASUM(M*N,STKI(L),STKI(L),1) .NE. 0.0D0) IMG = 1
      CALL SAVLOD(LUNIT,BUF,M,N,IMG,0,STKR(L),STKI(L))
      K = K-1
      IF (K .GE. BOT) GO TO 32
      CALL FILES(-LUNIT,BUF)
      MSTK(TOP) = 0
      GO TO 99
C 
C     LOAD
   35 IF (LUN .LT. 0) LUNIT = 2
      IF (LUN .LT. 0) CALL FILES(LUNIT,BUF)
      IF (LUN .GT. 0) LUNIT = LUN
   36 JOB = LSTK(BOT) - L
      CALL SAVLOD(LUNIT,ID,MSTK(TOP),NSTK(TOP),IMG,JOB,STKR(L),STKI(L))
      MN = MSTK(TOP)*NSTK(TOP)
      IF (MN .EQ. 0) GO TO 39
      IF (IMG .EQ. 0) CALL RSET(MN,0.0D0,STKI(L),1)
      DO 38 I = 1, 4
         J = 0
   37    J = J+1
         IF (ID(I).NE.ALFA(J) .AND. J.LE.BLANK) GO TO 37
         ID(I) = J-1
   38 CONTINUE
      SYM = SEMI
      RHS = 0
      CALL STACKP(ID)
      TOP = TOP + 1
      GO TO 36
   39 CALL FILES(-LUNIT,BUF)
      MSTK(TOP) = 0
      GO TO 99
C 
C     RAT
   40 IF (RHS .EQ. 2) GO TO 44
      MN = M*N
      L2 = L
      IF (LHS .EQ. 2) L2 = L + MN
      LW = L2 + MN
      ERR = LW + LRAT - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      IF (LHS .EQ. 2) TOP = TOP + 1
      LSTK(TOP) = L2
      MSTK(TOP) = M
      NSTK(TOP) = N
      CALL RSET(LHS*MN,0.0D0,STKI(L),1)
      DO 42 I = 1, MN
         CALL RAT(STKR(L),LRAT,MRAT,S,T,STKR(LW))
         STKR(L) = S
         STKR(L2) = T
         IF (LHS .EQ. 1) STKR(L) = FLOP(S/T)
         L = L + 1
         L2 = L2 + 1
   42 CONTINUE
      GO TO 99
   44 MRAT = IDINT(STKR(L))
      LRAT = IDINT(STKR(L-1))
      TOP = TOP - 1
      MSTK(TOP) = 0
      GO TO 99
C 
C     CHAR
   50 K = IABS(IDINT(STKR(L)))
      IF (M*N.NE.1 .OR. K.GE.ALFL) CALL ERROR(36)
      IF (ERR .GT. 0) RETURN
      CH = ALFA(K+1)
      IF (STKR(L) .LT. 0.0D0) CH = ALFB(K+1)
      WRITE(WTE,51) CH
   51 FORMAT(1X,'REPLACE CHARACTER ',A1)
      READ(RTE,52) CH
   52 FORMAT(A1)
      IF (STKR(L) .GE. 0.0D0) ALFA(K+1) = CH
      IF (STKR(L) .LT. 0.0D0) ALFB(K+1) = CH
      MSTK(TOP) = 0
      GO TO 99
C 
C     DISP
   60 WRITE(WTE,61)
      IF (WIO .NE. 0) WRITE(WIO,61)
   61 FORMAT(1X,80A1)
      IF (RHS .EQ. 2) GO TO 65
      MN = M*N
      TEXT = .TRUE.
      DO 62 I = 1, MN
        LS = L+I-1
        CH = IDINT(STKR(LS))
        TEXT = TEXT .AND. (CH.GE.0) .AND. (CH.LT.ALFL)
        TEXT = TEXT .AND. (DBLE(CH).EQ.STKR(LS))
   62 CONTINUE
      DO 64 I = 1, M
      DO 63 J = 1, N
        LS = L+I-1+(J-1)*M
        IF (STKR(LS) .EQ. 0.0D0) CH = BLANK
        IF (STKR(LS) .GT. 0.0D0) CH = PLUS
        IF (STKR(LS) .LT. 0.0D0) CH = MINUS
        IF (TEXT) CH = IDINT(STKR(LS))
        BUF(J) = ALFA(CH+1)
   63 CONTINUE
      WRITE(WTE,61) (BUF(J),J=1,N)
      IF (WIO .NE. 0) WRITE(WIO,61) (BUF(J),J=1,N)
   64 CONTINUE
      MSTK(TOP) = 0
      GO TO 99
C 
C     BASE
   65 IF (RHS .NE. 2) CALL ERROR(39)
      IF (STKR(L) .LE. 1.0D0) CALL ERROR(36)
      IF (ERR .GT. 0) RETURN
      B = STKR(L)
      L2 = L
      TOP = TOP-1
      RHS = 1
      L = LSTK(TOP)
      M = MSTK(TOP)*NSTK(TOP)
      EPS = STKR(VSIZE-4)
      DO 66 I = 1, M
         LS = L2+(I-1)*N
         LL = L+I-1
         CALL BASE(STKR(LL),B,EPS,STKR(LS),N)
   66 CONTINUE
      CALL RSET(M*N,0.0D0,STKI(L2),1)
      CALL WCOPY(M*N,STKR(L2),STKI(L2),1,STKR(L),STKI(L),1)
      MSTK(TOP) = N
      NSTK(TOP) = M
      CALL STACK1(QUOTE)
      IF (FIN .EQ. 6) GO TO 60
      GO TO 99
C 
C     LINES
   70 LCT(2) = IDINT(STKR(L))
      MSTK(TOP) = 0
      GO TO 99
C 
C     PLOT
   80 IF (RHS .GE. 2) GO TO 82
      N = M*N
      DO 81 I = 1, N
         LL = L+I-1
         STKI(LL) = DBLE(I)
   81 CONTINUE
      CALL PLOT(WTE,STKI(L),STKR(L),N,T,0,BUF)
      IF (WIO .NE. 0) CALL PLOT(WIO,STKI(L),STKR(L),N,T,0,BUF)
      MSTK(TOP) = 0
      GO TO 99
   82 IF (RHS .EQ. 2) K = 0
      IF (RHS .EQ. 3) K = M*N
      IF (RHS .GT. 3) K = RHS - 2
      TOP = TOP - (RHS - 1)
      N = MSTK(TOP)*NSTK(TOP)
      IF (MSTK(TOP+1)*NSTK(TOP+1) .NE. N) CALL ERROR(5)
      IF (ERR .GT. 0) RETURN
      LX = LSTK(TOP)
      LY = LSTK(TOP+1)
      IF (RHS .GT. 3) L = LSTK(TOP+2)
      CALL PLOT(WTE,STKR(LX),STKR(LY),N,STKR(L),K,BUF)
      IF (WIO .NE. 0) CALL PLOT(WIO,STKR(LX),STKR(LY),N,STKR(L),K,BUF)
      MSTK(TOP) = 0
      GO TO 99
C 
C     DEBUG
   95 DDT = IDINT(STKR(L))
      WRITE(WTE,96) DDT
   96 FORMAT(1X,'DEBUG ',I4)
      MSTK(TOP) = 0
      GO TO 99
C 
   99 RETURN
      END
