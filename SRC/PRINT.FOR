      SUBROUTINE PRINT(ID,K)
C     PRIMARY OUTPUT ROUTINE
      INTEGER ID(4),K
      DOUBLE PRECISION STKR(5005),STKI(5005)
      INTEGER IDSTK(4,48),LSTK(48),MSTK(48),NSTK(48),VSIZE,LSIZE,BOT,TOP
      INTEGER ALFA(52),ALFB(52),ALFL,CASE
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /VSTK/ STKR,STKI,IDSTK,LSTK,MSTK,NSTK,VSIZE,LSIZE,BOT,TOP
      COMMON /ALFS/ ALFA,ALFB,ALFL,CASE
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      DOUBLE PRECISION S,TR,TI,PR(12),PI(12),ROUNDM
      INTEGER FNO(11),FNL(11),SIG(12),PLUS,MINUS,BLANK,TYP,F
      DATA PLUS/41/,MINUS/42/,BLANK/36/
C     FORMAT NUMBERS AND LENGTHS
      DATA FNO /11,12,21,22,23,24,31,32,33,34,-1/
      DATA FNL /12, 6, 8, 4, 6, 3, 4, 2, 3, 1, 1/
C     FMT   1       2       3       4       5
C         SHORT   LONG   SHORT E  LONG E    Z
C     TYP   1       2       3
C         INTEGER  REAL   COMPLEX
      IF (LCT(1) .LT. 0) GO TO 99
      L = LSTK(K)
      M = MSTK(K)
      N = NSTK(K)
      MN = M*N
      TYP = 1
      S = 0.0D0
      DO 10 I = 1, MN
        LS = L+I-1
        TR = STKR(LS)
        TI = STKI(LS)
        S = DMAX1(S,DABS(TR),DABS(TI))
        IF (ROUNDM(TR) .NE. TR) TYP = MAX0(2,TYP)
        IF (TI .NE. 0.0D0) TYP = 3
   10 CONTINUE
      IF (S .NE. 0.0D0) S = DLOG10(S)
      KS = IDINT(S)
      IF (-2 .LE. KS .AND. KS .LE. 1) KS = 0
      IF (KS .EQ. 2 .AND. FMT .EQ. 1 .AND. TYP .EQ. 2) KS = 0
      IF (TYP .EQ. 1 .AND. KS .LE. 2) F = 1
      IF (TYP .EQ. 1 .AND. KS .GT. 2) F = 2
      IF (TYP .EQ. 1 .AND. KS .GT. 9) TYP = 2
      IF (TYP .EQ. 2) F = FMT + 2
      IF (TYP .EQ. 3) F = FMT + 6
      IF (MN.EQ.1 .AND. KS.NE.0 .AND. FMT.LT.3 .AND. TYP.NE.1) F = F+2
      IF (FMT .EQ. 5) F = 11
      JINC = FNL(F)
      F = FNO(F)
      S = 1.0D0
      IF (F.EQ.21 .OR. F.EQ.22 .OR. F.EQ.31 .OR. F.EQ.32) S = 10.0D0**KS
      LS = ((N-1)/JINC+1)*M + 2
      IF (LCT(1) + LS .LE. LCT(2)) GO TO 20
         LCT(1) = 0
         WRITE(WTE,43) LS
         READ(RTE,44,END=19) LS
CDC..    IF (EOF(RTE).NE.0) GO TO 19
         IF (LS .EQ. ALFA(BLANK+1)) GO TO 20
         LCT(1) = -1
         GO TO 99
   19    CALL FILES(-RTE,BUF)
   20 CONTINUE
      WRITE(WTE,44)
      IF (WIO .NE. 0) WRITE(WIO,44)
      CALL PRNTID(ID,-1)
      LCT(1) = LCT(1)+2
      LUNIT = WTE
   50 IF (S .NE. 1.0D0) WRITE(LUNIT,41) S
      DO 80 J1 = 1, N, JINC
        J2 = MIN0(N, J1+JINC-1)
        WRITE(LUNIT,44)
        IF (N .GT. JINC) WRITE(LUNIT,42) J1,J2
        DO 70 I = 1, M
          JM = J2-J1+1
          DO 60 J = 1, JM
             LS = L+I-1+(J+J1-2)*M
             PR(J) = STKR(LS)/S
             PI(J) = DABS(STKI(LS)/S)
             SIG(J) = ALFA(PLUS+1)
             IF (STKI(LS) .LT. 0.0D0) SIG(J) = ALFA(MINUS+1)
   60     CONTINUE
          IF (F .EQ. 11) WRITE(LUNIT,11)(PR(J),J=1,JM)
          IF (F .EQ. 12) WRITE(LUNIT,12)(PR(J),J=1,JM)
          IF (F .EQ. 21) WRITE(LUNIT,21)(PR(J),J=1,JM)
          IF (F .EQ. 22) WRITE(LUNIT,22)(PR(J),J=1,JM)
          IF (F .EQ. 23) WRITE(LUNIT,23)(PR(J),J=1,JM)
          IF (F .EQ. 24) WRITE(LUNIT,24)(PR(J),J=1,JM)
          IF (F .EQ. 31) WRITE(LUNIT,31)(PR(J),SIG(J),PI(J),J=1,JM)
          IF (F .EQ. 32) WRITE(LUNIT,32)(PR(J),SIG(J),PI(J),J=1,JM)
          IF (F .EQ. 33) WRITE(LUNIT,33)(PR(J),SIG(J),PI(J),J=1,JM)
          IF (F .EQ. 34) WRITE(LUNIT,34)(PR(J),SIG(J),PI(J),J=1,JM)
          IF (F .EQ. -1) CALL FORMZ(LUNIT,STKR(LS),STKI(LS))
          LCT(1) = LCT(1)+1
   70   CONTINUE
   80 CONTINUE
      IF (LUNIT.EQ.WIO .OR. WIO.EQ.0) GO TO 99
      LUNIT = WIO
      GO TO 50
   99 RETURN
C 
   11 FORMAT(1X,12F6.0)
   12 FORMAT(1X,6F12.0)
   21 FORMAT(1X,F9.4,7F10.4)
   22 FORMAT(1X,F19.15,3F20.15)
   23 FORMAT(1X,1P6D13.4)
   24 FORMAT(1X,1P3D24.15)
   31 FORMAT(1X,4(F9.4,' ',A1,F7.4,'i'))
   32 FORMAT(1X,F19.15,A1,F18.15,'i',F20.15,A1,F18.15,'i')
   33 FORMAT(1X,3(1PD13.4,' ',A1,1PD10.4,'i'))
   34 FORMAT(1X,1PD24.15,' ',A1,1PD21.15,'i')
   41 FORMAT(/1X,' ',1PD9.1,2H *)
   42 FORMAT(1X,'    COLUMNS',I3,' THRU',I3)
   43 FORMAT(/1X,'AT LEAST ',I5,' MORE LINES.',
     $       '  ENTER BLANK LINE TO CONTINUE OUTPUT.')
   44 FORMAT(A1)
C 
      END
