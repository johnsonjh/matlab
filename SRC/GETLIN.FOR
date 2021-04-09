      SUBROUTINE GETLIN
C     GET A NEW LINE
      INTEGER ALFA(52),ALFB(52),ALFL,CASE
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /ALFS/ ALFA,ALFB,ALFL,CASE
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
      CHARACTER*1 CBUF(256)
      INTEGER LRECL,EOL,SLASH,BSLASH,DOT,BLANK,RETU(4)
      DATA EOL/99/,DOT/47/,BLANK/36/,RETU/27,14,29,30/
      DATA SLASH/44/,BSLASH/45/,LRECL/80/
C
   10 L = LPT(1)
   11 DO 12 J = 1, LRECL
         BUF(J) = ALFA(BLANK+1)
   12 CONTINUE
      READ(RIO,101,END=50,ERR=15) (CBUF(J),J=1,LRECL)
CDC.. IF (EOF(RIO).NE.0) GO TO 50
  101 FORMAT(80A1)
      DO 13 J = 1, LRECL
         BUF(J) = ICHAR(CBUF(J))
   13 CONTINUE
      N = LRECL+1
   15 N = N-1
      IF (BUF(N) .EQ. ALFA(BLANK+1)) GO TO 15
      IF (MOD(LCT(4),2) .EQ. 1) WRITE(WTE,102) (BUF(J),J=1,N)
      IF (WIO .NE. 0) WRITE(WIO,102) (BUF(J),J=1,N)
  102 FORMAT(1X,80A1)
C
      DO 40 J = 1, N
         DO 20 K = 1, ALFL
           IF (BUF(J).EQ.ALFA(K) .OR. BUF(J).EQ.ALFB(K)) GO TO 30
   20    CONTINUE
         K = EOL+1
         CALL XCHAR(BUF(J),K)
         IF (K .GT. EOL) GO TO 10
         IF (K .EQ. EOL) GO TO 45
         IF (K .EQ. -1) L = L-1
         IF (K .LE. 0) GO TO 40
C
   30    K = K-1
         IF (K.EQ.SLASH .AND. BUF(J+1).EQ.BUF(J)) GO TO 45
         IF (K.EQ.DOT .AND. BUF(J+1).EQ.BUF(J)) GO TO 11
         IF (K.EQ.BSLASH .AND. N.EQ.1) GO TO 60
         LIN(L) = K
         IF (L.LT.1024) L = L+1
         IF (L.EQ.1024) WRITE(WTE,33) L
   33    FORMAT(1X,'INPUT BUFFER LIMIT IS ',I4,' CHARACTERS.')
   40 CONTINUE
   45 LIN(L) = EOL
      LPT(6) = L
      LPT(4) = LPT(1)
      LPT(3) = 0
      LPT(2) = 0
      LCT(1) = 0
      CALL GETCH
      RETURN
C
   50 IF (RIO .EQ. RTE) GO TO 52
      CALL PUTID(LIN(L),RETU)
      L = L + 4
      GO TO 45
   52 CALL FILES(-RTE,BUF)
      LIN(L) = EOL
      RETURN
C
   60 N = LPT(6) - LPT(1)
      DO 61 I = 1, N
         J = L+I-1
         K = LIN(J)
         BUF(I) = ALFA(K+1)
         IF (CASE.EQ.1 .AND. K.LT.36) BUF(I) = ALFB(K+1)
   61 CONTINUE
      CALL EDIT(BUF,N)
      N = N + 1
      GO TO 15
      END
