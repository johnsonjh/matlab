      SUBROUTINE COMANDM(ID)
      INTEGER ID(4)
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
      CHARACTER*1 CBUF(256)
      INTEGER CMD(4,17),CMDL,A,D,E,Z,LRECL,CH,BLANK,NAME,DOT,H(4)
      INTEGER SEMI,COMMA,EOL
      DOUBLE PRECISION URAND
      LOGICAL EQID
      DATA CMDL/17/,A/10/,D/13/,E/14/,Z/35/,EOL/99/,SEMI/39/,COMMA/48/
      DATA BLANK/36/,NAME/1/,DOT/47/
C
C       CLEAR ELSE  END   EXIT
C       FOR   HELP  IF    LONG
C       RETUR SEMI
C       SHORT WHAT  WHILE
C       WHO   WHY   LALA  FOO
C
      DATA CMD/
     $  12,21,14,10, 14,21,28,14, 14,23,13,36, 14,33,18,29,
     $  15,24,27,36, 17,14,21,25, 18,15,36,36, 21,24,23,16,
     $  27,14,29,30, 28,14,22,18,
     $  28,17,24,27, 32,17,10,29, 32,17,18,21,
     $  32,17,24,36, 32,17,34,36, 21,10,21,10, 15,30,12,20/
C
      DATA LRECL/80/
  101 FORMAT(80A1)
  102 FORMAT(1X,80A1)
C
      IF (DDT .EQ. 1) WRITE(WTE,100)
  100 FORMAT(1X,'COMAND')
      FUN = 0
      DO 10 K = 1, CMDL
C
C     Microsoft FORTRAN Error will not compile the following line correctly       C
C       IF (EQID(CMD(1,K), ID)) GO TO 20
C
        if (eqid(cmd(1,k), id(1))) go to 20
   10 CONTINUE
      FIN = 0
      RETURN
C
   20 IF (CHAR.EQ.COMMA .OR. CHAR.EQ.SEMI .OR. CHAR.EQ.EOL) GO TO 22
      IF (CHAR.LE.Z .OR. K.EQ.6) GO TO 22
      CALL ERROR(16)
      RETURN
C
   22 FIN = 1
      GO TO (25,36,38,40,30,80,34,52,44,55,50,65,32,60,70,46,48),K
C
C     CLEAR
   25 IF (CHAR.GE.A .AND. CHAR.LE.Z) GO TO 26
      BOT = LSIZE-3
      GO TO 98
   26 CALL GETSYM
      TOP = TOP+1
      MSTK(TOP) = 0
      NSTK(TOP) = 0
      RHS = 0
      CALL STACKP(SYN)
      IF (ERR .GT. 0) RETURN
      FIN = 1
      GO TO 98
C
C     FOR, WHILE, IF, ELSE, END
   30 FIN = -11
      GO TO 99
   32 FIN = -12
      GO TO 99
   34 FIN = -13
      GO TO 99
   36 FIN = -14
      GO TO 99
   38 FIN = -15
      GO TO 99
C
C     EXIT
   40 IF (PT .GT. PTZ) FIN = -16
      IF (PT .GT. PTZ) GO TO 98
      K = IDINT(STKR(VSIZE-2))
      WRITE(WTE,140) K
      IF (WIO .NE. 0) WRITE(WIO,140) K
  140 FORMAT(/1X,'total flops ',I9//1X,'ADIOS'/)
      FUN = 99
      GO TO 98
C
C     RETURN
   44 K = LPT(1) - 7
      IF (K .LE. 0) FUN = 99
      IF (K .LE. 0) GO TO 98
      CALL FILES(-RIO,BUF)
      LPT(1) = LIN(K+1)
      LPT(4) = LIN(K+2)
      LPT(6) = LIN(K+3)
      PTZ = LIN(K+4)
      RIO = LIN(K+5)
      LCT(4) = LIN(K+6)
      CHAR = BLANK
      SYM = COMMA
      GO TO 99
C
C     LALA
   46 WRITE(WTE,146)
  146 FORMAT(1X,'QUIT SINGING AND GET BACK TO WORK.')
      GO TO 98
C
C     FOO
   48 WRITE(WTE,148)
  148 FORMAT(1X,'YOUR PLACE OR MINE')
      GO TO 98
C
C     SHORT, LONG
   50 FMT = 1
      GO TO 54
   52 FMT = 2
   54 IF (CHAR.EQ.E .OR. CHAR.EQ.D) FMT = FMT+2
      IF (CHAR .EQ. Z) FMT = 5
      IF (CHAR.EQ.E .OR. CHAR.EQ.D .OR. CHAR.EQ.Z) CALL GETSYM
      GO TO 98
C
C     SEMI
   55 LCT(3) = 1 - LCT(3)
      GO TO 98
C
C     WHO
   60 WRITE(WTE,160)
      IF (WIO .NE. 0) WRITE(WIO,160)
  160 FORMAT(1X,'Your current variables are...')
      CALL PRNTID(IDSTK(1,BOT),LSIZE-BOT+1)
      L = VSIZE-LSTK(BOT)+1
      WRITE(WTE,161) L,VSIZE
      IF (WIO .NE. 0) WRITE(WIO,161) L,VSIZE
  161 FORMAT(1X,'using ',I7,' out of ',I7,' elements.')
C
C       PRINT OUT THE MATRICES AVAILIABLE IN TOTAL2
C
C      CALL MATTOTDIR
      GO TO 98
C
C     WHAT
   65 WRITE(WTE,165)
  165 FORMAT(1X,'The functions and commands are...')
      H(1) = 0
      CALL FUNS(H)
      CALL PRNTID(CMD,CMDL-2)
      GO TO 98
C
C     WHY
   70 K = IDINT(9.0D0*URAND(RAN(1))+1.0D0)
      GO TO (71,72,73,74,75,76,77,78,79),K
   71 WRITE(WTE,171)
  171 FORMAT(1X,'WHAT?')
      GO TO 98
   72 WRITE(WTE,172)
  172 FORMAT(1X,'R.T.F.M.')
      GO TO 98
   73 WRITE(WTE,173)
  173 FORMAT(1X,'HOW DO YOU EXPECT ME TO KNOW?')
      GO TO 98
   74 WRITE(WTE,174)
  174 FORMAT(1X,'PETE MADE ME DO IT.')
      GO TO 98
   75 WRITE(WTE,175)
  175 FORMAT(1X,'INSUFFICIENT DATA TO ANSWER.')
      GO TO 98
   76 WRITE(WTE,176)
  176 FORMAT(1X,'IT FEELS GOOD.')
      GO TO 98
   77 WRITE(WTE,177)
  177 FORMAT(1X,'WHY NOT?')
      GO TO 98
   78 WRITE(WTE,178)
  178 FORMAT(1X,'/--ERROR'/1X,'STUPID QUESTION.')
      GO TO 98
   79 WRITE(WTE,179)
  179 FORMAT(1X,'SYSTEM ERROR, RETRY')
      GO TO 98
C
C     HELP
   80 IF (CHAR .NE. EOL) GO TO 81
      WRITE(WTE,180)
      IF (WIO .NE. 0) WRITE(WIO,180)
  180 FORMAT(1X,'Type HELP followed by ...'
     $  /1X,'INTRO   (To get started)'
     $  /1X,'NEWS    (recent revisions)')
      H(1) = 0
      CALL FUNS(H)
      CALL PRNTID(CMD,CMDL-2)
      J = BLANK+2
      WRITE(WTE,181)
      IF (WIO .NE. 0) WRITE(WIO,181)
  181 FORMAT(1X,'ANS   EDIT  FILE  FUN   MACRO')
      WRITE(WTE,182) (ALFA(I),I=J,ALFL)
      IF (WIO .NE. 0) WRITE(WIO,182) (ALFA(I),I=J,ALFL)
  182 FORMAT(1X,17(A1,1X)/)
      GO TO 98
C
   81 CALL GETSYM
      IF (SYM .EQ. NAME) GO TO 82
      IF (SYM .EQ. 0) SYM = DOT
      H(1) = ALFA(SYM+1)
      H(2) = ALFA(BLANK+1)
      H(3) = ALFA(BLANK+1)
      H(4) = ALFA(BLANK+1)
      GO TO 84
   82 DO 83 I = 1, 4
        CH = SYN(I)
        H(I) = ALFA(CH+1)
   83 CONTINUE
   84 READ(HIO,101,END=89) (CBUF(I),I=1,LRECL)
      DO 1112 I = 1,LRECL
         BUF(I) = ICHAR(CBUF(I))
 1112 CONTINUE
CDC.. IF (EOF(HIO).NE.0) GO TO 89
      DO 85 I = 1, 4
        IF (H(I) .NE. BUF(I)) GO TO 84
   85 CONTINUE
      WRITE(WTE,102)
      IF (WIO .NE. 0) WRITE(WIO,102)
   86 K = LRECL + 1
   87 K = K - 1
      IF (BUF(K) .EQ. ALFA(BLANK+1)) GO TO 87
      WRITE(WTE,102) (CBUF(I),I=1,K)
      IF (WIO .NE. 0) WRITE(WIO,102) (CBUF(I),I=1,K)
      READ(HIO,101) (CBUF(I),I=1,LRECL)
      DO 88 I = 1, LRECL
         BUF(I) = ICHAR(CBUF(I))
   88 CONTINUE
      IF (BUF(1) .EQ. ALFA(BLANK+1)) GO TO 86
      CALL FILES(-HIO,BUF)
      GO TO 98
C
   89 WRITE(WTE,189) (H(I),I=1,4)
  189 FORMAT(1X,'SORRY, NO HELP ON ',4A1)
      CALL FILES(-HIO,BUF)
      GO TO 98
C
   98 CALL GETSYM
   99 RETURN
      END
