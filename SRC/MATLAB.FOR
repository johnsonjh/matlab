      SUBROUTINE MATLAB(INIT)
C
C     INIT = 0 FOR ORDINARY FIRST ENTRY
C          = POSITIVE FOR SUBSEQUENT ENTRIES
C          = NEGATIVE FOR SILENT INITIALIZATION (SEE MATZ)
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
C
      DOUBLE PRECISION S,T
      INTEGER EPS(4),FLOPS(4),EYE(4),RAND(4)
C
C     CHARACTER SET
C            0       10       20       30       40       50
C
C     0      0        A        K        U   COLON  :  LESS   <
C     1      1        B        L        V   PLUS   +  GREAT  >
C     2      2        C        M        W   MINUS  -
C     3      3        D        N        X   STAR   *
C     4      4        E        O        Y   SLASH  /
C     5      5        F        P        Z   BSLASH \
C     6      6        G        Q  BLANK     EQUAL  =
C     7      7        H        R  LPAREN (  DOT    .
C     8      8        I        S  RPAREN )  COMMA  ,
C     9      9        J        T  SEMI   ;  QUOTE  '
C
      CHARACTER ALPHA(52),ALPHB(52)
      DATA ALPHA /'0','1','2','3','4','5','6','7','8','9',
     $    'A','B','C','D','E','F','G','H','I','J',
     $    'K','L','M','N','O','P','Q','R','S','T',
     $    'U','V','W','X','Y','Z',' ','(',')',';',
     $    ':','+','-','*','/','\','=','.',',','''',
     $    '<','>'/
C
C     ALTERNATE CHARACTER SET
C
      DATA ALPHB /'0','1','2','3','4','5','6','7','8','9',
     $    'a','b','c','d','e','f','g','h','i','j',
     $    'k','l','m','n','o','p','q','r','s','t',
     $    'u','v','w','x','y','z',' ','(',')',';',
     $    '|','+','-','*','/','$','=','.',',','"',
     $    '[',']'/
C
      DATA EPS/14,25,28,36/,FLOPS/15,21,24,25/
      DATA EYE/14,34,14,36/,RAND/27,10,23,13/
C
      IF (INIT .GT. 0) GO TO 90
C
C     RTE = UNIT NUMBER FOR TERMINAL INPUT
      RTE = 0
      CALL FILES(RTE,BUF)
      RIO = RTE
C
C     WTE = UNIT NUMBER FOR TERMINAL OUTPUT
      WTE = 0
      CALL FILES(WTE,BUF)
      WIO = 0
C
      IF (INIT .GE. 0) WRITE(WTE,100)
  100 FORMAT(//1X,'     < M A T L A B >'
     $  /1X,'   Version of 05/25/82')
C
C     HIO = UNIT NUMBER FOR HELP FILE
      HIO = 9
      CALL FILES(HIO,BUF)
C
C     RANDOM NUMBER SEED
      RAN(1) = 0
C
C     INITIAL LINE LIMIT
      LCT(2) = 25
C
      ALFL = 52
      CASE = 0
C     CASE = 1 for file names in lower case
      DO 20 I = 1, ALFL
         ALFA(I) = ICHAR(ALPHA(I))
         ALFB(I) = ICHAR(ALPHB(I))
   20 CONTINUE
C
      VSIZE = 5005
      LSIZE = 48
      PSIZE = 32
      BOT = LSIZE-3
      CALL WSET(5,0.0D0,0.0D0,STKR(VSIZE-4),STKI(VSIZE-4),1)
      CALL PUTID(IDSTK(1,LSIZE-3),EPS)
      LSTK(LSIZE-3) = VSIZE-4
      MSTK(LSIZE-3) = 1
      NSTK(LSIZE-3) = 1
      S = 1.0D0
   30 S = S/2.0D0
      T = 1.0D0 + S
      IF (T .GT. 1.0D0) GO TO 30
      STKR(VSIZE-4) = 2.0D0*S
      CALL PUTID(IDSTK(1,LSIZE-2),FLOPS)
      LSTK(LSIZE-2) = VSIZE-3
      MSTK(LSIZE-2) = 1
      NSTK(LSIZE-2) = 2
      CALL PUTID(IDSTK(1,LSIZE-1), EYE)
      LSTK(LSIZE-1) = VSIZE-1
      MSTK(LSIZE-1) = -1
      NSTK(LSIZE-1) = -1
      STKR(VSIZE-1) = 1.0D0
      CALL PUTID(IDSTK(1,LSIZE), RAND)
      LSTK(LSIZE) = VSIZE
      MSTK(LSIZE) = 1
      NSTK(LSIZE) = 1
      FMT = 1
      FLP(1) = 0
      FLP(2) = 0
      DDT = 0
      RAN(2) = 0
      PTZ = 0
      PT = PTZ
      ERR = 0
      IF (INIT .LT. 0) RETURN
C
   90 CALL PARSE
      IF (FUN .EQ. 1) CALL MATFN1
      IF (FUN .EQ. 2) CALL MATFN2
      IF (FUN .EQ. 3) CALL MATFN3
      IF (FUN .EQ. 4) CALL MATFN4
      IF (FUN .EQ. 5) CALL MATFN5
      IF (FUN .EQ. 6) CALL MATFN6
      IF (FUN .EQ. 21) CALL MATFN1
      IF (FUN .NE. 99) GO TO 90
      RETURN
      END
