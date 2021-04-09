      SUBROUTINE CLAUSE
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER IDS(4,32),PSTK(32),RSTK(32),PSIZE,PT,PTZ
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /RECU/ IDS,PSTK,RSTK,PSIZE,PT,PTZ
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      INTEGER FOR(4),WHILE(4),IFF(4),ELSE(4),ENND(4),DO(4),THENN(4)
      INTEGER SEMI,EQUAL,EOL,BLANK,R
      INTEGER OP,COMMA,LESS,GREAT,NAME
      LOGICAL EQID
      DOUBLE PRECISION E1,E2
      DATA SEMI/39/,EQUAL/46/,EOL/99/,BLANK/36/
      DATA COMMA/48/,LESS/50/,GREAT/51/,NAME/1/
      DATA FOR/15,24,27,36/,WHILE/32,17,18,21/,IFF/18,15,36,36/
      DATA ELSE/14,21,28,14/,ENND/14,23,13,36/
      DATA DO/13,24,36,36/,THENN/29,17,14,23/
      R = -FIN-10
      FIN = 0
      IF (DDT .EQ. 1) WRITE(WTE,100) PT,RSTK(PT),R
  100 FORMAT(1X,'CLAUSE',3I4)
      IF (R.LT.1 .OR. R.GT.6) GO TO 01
      GO TO (02,30,30,80,99,90),R
   01 R = RSTK(PT)
      GO TO (99,99,05,40,45,99,99,99,99,99,99,99,15,55,99,99,99),R
C 
C     FOR
C 
   02 CALL GETSYM
      IF (SYM .NE. NAME) CALL ERROR(34)
      IF (ERR .GT. 0) RETURN
      PT = PT+2
      CALL PUTID(IDS(1,PT),SYN)
      CALL GETSYM
      IF (SYM .NE. EQUAL) CALL ERROR(34)
      IF (ERR .GT. 0) RETURN
      CALL GETSYM
      RSTK(PT) = 3
C     *CALL* EXPR
      RETURN
   05 PSTK(PT-1) = 0
      PSTK(PT) = LPT(4) - 1
C
C     Comment:  Microsoft FORTRAN error labels EQID as an array.
C
C     IF (EQID(SYN,DO)) SYM = SEMI
C
      if (eqid(syn(1),do(1))) sym = semi
      IF (SYM .EQ. COMMA) SYM = SEMI
      IF (SYM .NE. SEMI) CALL ERROR(34)
      IF (ERR .GT. 0) RETURN
   10 J = PSTK(PT-1)
      LPT(4) = PSTK(PT)
      SYM = SEMI
      CHAR = BLANK
      J = J+1
      L = LSTK(TOP)
      M = MSTK(TOP)
      N = NSTK(TOP)
      LJ = L+(J-1)*M
      L2 = L + M*N
      IF (M .NE. -3) GO TO 12
      LJ = L+3
      L2 = LJ
      STKR(LJ) = STKR(L) + DBLE(J-1)*STKR(L+1)
      STKI(LJ) = 0.0
      IF (STKR(L+1).GT.0.0D0 .AND. STKR(LJ).GT.STKR(L+2)) GO TO 20
      IF (STKR(L+1).LT.0.0D0 .AND. STKR(LJ).LT.STKR(L+2)) GO TO 20
      M = 1
      N = J
   12 IF (J .GT. N) GO TO 20
      IF (TOP+1 .GE. BOT) CALL ERROR(18)
      IF (ERR .GT. 0) RETURN
      TOP = TOP+1
      LSTK(TOP) = L2
      MSTK(TOP) = M
      NSTK(TOP) = 1
      ERR = L2+M - LSTK(BOT)
      IF (ERR .GT. 0) CALL ERROR(17)
      IF (ERR .GT. 0) RETURN
      CALL WCOPY(M,STKR(LJ),STKI(LJ),1,STKR(L2),STKI(L2),1)
      RHS = 0
      CALL STACKP(IDS(1,PT))
      IF (ERR .GT. 0) RETURN
      PSTK(PT-1) = J
      PSTK(PT) = LPT(4)
      RSTK(PT) = 13
C     *CALL* PARSE
      RETURN
   15 GO TO 10
   20 MSTK(TOP) = 0
      NSTK(TOP) = 0
      RHS = 0
      CALL STACKP(IDS(1,PT))
      IF (ERR .GT. 0) RETURN
      PT = PT-2
      GO TO 80
C
C     WHILE OR IF
C
   30 PT = PT+1
      CALL PUTID(IDS(1,PT),SYN)
      PSTK(PT) = LPT(4)-1
   35 LPT(4) = PSTK(PT)
      CHAR = BLANK
      CALL GETSYM
      RSTK(PT) = 4
C     *CALL* EXPR
      RETURN
   40 IF (SYM.NE.EQUAL .AND. SYM.NE.LESS .AND. SYM.NE.GREAT)
     $    CALL ERROR(35)
      IF (ERR .GT. 0) RETURN
      OP = SYM
      CALL GETSYM
      IF (SYM.EQ.EQUAL .OR. SYM.EQ.GREAT) OP = OP + SYM
      IF (OP .GT. GREAT) CALL GETSYM
      PSTK(PT) = 256*PSTK(PT) + OP
      RSTK(PT) = 5
C     *CALL* EXPR
      RETURN
   45 OP = MOD(PSTK(PT),256)
      PSTK(PT) = PSTK(PT)/256
      L = LSTK(TOP-1)
      E1 = STKR(L)
      L = LSTK(TOP)
      E2 = STKR(L)
      TOP = TOP - 2
c
c     Microsoft FORTRAN Error does not correctly compile the following line
c
c     IF (EQID(SYN,DO) .OR. EQID(SYN,THENN)) SYM = SEMI
c
      if (eqid(syn(1),do(1)) .or. eqid(syn(1),thenn(1))) sym = semi
      IF (SYM .EQ. COMMA) SYM = SEMI
      IF (SYM .NE. SEMI) CALL ERROR(35)
      IF (ERR .GT. 0) RETURN
      IF (OP.EQ.EQUAL         .AND. E1.EQ.E2) GO TO 50
      IF (OP.EQ.LESS          .AND. E1.LT.E2) GO TO 50
      IF (OP.EQ.GREAT         .AND. E1.GT.E2) GO TO 50
      IF (OP.EQ.(LESS+EQUAL)  .AND. E1.LE.E2) GO TO 50
      IF (OP.EQ.(GREAT+EQUAL) .AND. E1.GE.E2) GO TO 50
      IF (OP.EQ.(LESS+GREAT)  .AND. E1.NE.E2) GO TO 50
      PT = PT-1
      GO TO 80
   50 RSTK(PT) = 14
C     *CALL* PARSE
      RETURN
c
c     Microsoft FORTRAN Error does not correctly compile the following line
c
c  55 IF (EQID(IDS(1,PT),WHILE)) GO TO 35
c
   55 if (eqid(ids(1,pt),while(1))) go to 35
      PT = PT-1
c     Microsoft FORTRAN Error does not correctly compile the following line
c
c
c     IF (EQID(SYN,ELSE)) GO TO 80
c
      if (eqid(syn(1), else(1))) go to 80
      RETURN
C
C     SEARCH FOR MATCHING END OR ELSE
   80 KOUNT = 0
      CALL GETSYM
   82 IF (SYM .EQ. EOL) RETURN
      IF (SYM .NE. NAME) GO TO 83
c
c     Microsoft FORTRAN Error does not correctly compile the following line
c
c     IF (EQID(SYN,ENND) .AND. KOUNT.EQ.0) RETURN
c     IF (EQID(SYN,ELSE) .AND. KOUNT.EQ.0) RETURN
c     IF (EQID(SYN,ENND) .OR. EQID(SYN,ELSE))
c    $       KOUNT = KOUNT-1
c     IF (EQID(SYN,FOR) .OR. EQID(SYN,WHILE)
c    $       .OR. EQID(SYN,IFF)) KOUNT = KOUNT+1
c
      if (eqid(syn(1),ennd(1)) .and. kount .eq. 0) return
      if (eqid(syn(1),else(1)) .and. kount .eq. 0) return
      if (eqid(syn(1),ennd(1)) .or. eqid(syn(1),else(1)))
     $       kount = kount-1
      if (eqid(syn(1),for(1)) .or. eqid(syn(1),while(1))
     $       .or. eqid(syn(1),iff(1))) kount = kount+1
   83 CALL GETSYM
      GO TO 82
C
C     EXIT FROM LOOP
   90 IF (DDT .EQ. 1) WRITE(WTE,190) (RSTK(I),I=1,PT)
  190 FORMAT(1X,'EXIT  ',10I4)
      IF (RSTK(PT) .EQ. 14) PT = PT-1
      IF (PT .LE. PTZ) RETURN
      IF (RSTK(PT) .EQ. 14) PT = PT-1
      IF (PT-1 .LE. PTZ) RETURN
      IF (RSTK(PT) .EQ. 13) TOP = TOP-1
      IF (RSTK(PT) .EQ. 13) PT = PT-2
      GO TO 80
C
   99 CALL ERROR(22)
      IF (ERR .GT. 0) RETURN
      RETURN
      END
