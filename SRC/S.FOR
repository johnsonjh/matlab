C     PROGRAM MAIN FOR VAX/VMS
C     For overflow handling, see
C     Volume 2A, Common Run-Time Library, page 6-40.
C     Also compile subroutine URAND with FORTRAN/NOCHECK.
      EXTERNAL UNTRAP
      CALL LIB$ESTABLISH(UNTRAP)
      CALL MATLAB(0)
      STOP
      END

 
      INTEGER FUNCTION UNTRAP(SIGARGS,MECHARGS)
C     VAX/VMS, Vol. 2A, Common Run-Time, p. 6-40.
      INTEGER SIGARGS(3),MECHARGS(5)
      INCLUDE 'SYS$LIBRARY:SIGDEF'
      INTEGER HUGE(2)
      DATA HUGE/'FFFF7FFF'X,'FFFFFFFF'X/
      IF (LIB$MATCH_COND (SIGARGS(2), SS$_FLTOVF)) THEN
         WRITE(6,10)
   10    FORMAT(' ---OVERFLOW---')
         CALL LIB$INSV(1,0,3,SIGARGS(2))
         UNTRAP = SS$_CONTINUE
      ELSEIF (LIB$MATCH_COND (SIGARGS(2), SS$_INTOVF)) THEN
         WRITE(6,12)
   12    FORMAT(' ---INTEGER OVERFLOW---')
         CALL LIB$INSV(2,0,3,SIGARGS(2))
         UNTRAP = SS$_CONTINUE
      ELSEIF (LIB$MATCH_COND (SIGARGS(2), SS$_ROPRAND)) THEN
         UNTRAP = LIB$FIXUP_FLT(SIGARGS,MECHARGS,HUGE)
      ELSE
         UNTRAP = SS$_RESIGNAL
      ENDIF
      RETURN
      END

 
      SUBROUTINE FILES(LUNIT,NAME)
      INTEGER LUNIT,NAME(32)
C 
C     SYSTEM DEPENDENT ROUTINE TO ALLOCATE FILES
C     LUNIT = LOGICAL UNIT NUMBER
C     NAME = FILE NAME, 1 CHARACTER PER WORD
C 
      character*32 nam
C 
      if (lunit .eq. 5) return
      if (lunit .eq. 6) return
      if (lunit .eq. 9) then
         open(unit=9,name='math$doc:mathelp.lis',
     $        type='old',shared,readonly)
         write(6,09)
   09    format(/1X,'HELP is available')
         return
      end if
      if (lunit .eq. -9) then
         rewind 9
         return
      end if
      if (lunit .lt. 0) then
         close(unit=-lunit)
         return
      end if
      do 10 i = 1,32
         nam(i:i) = char(name(i))
   10 continue
      open(unit=lunit,name=nam,type='unknown',recordsize=81)
      RETURN
      END

 
      SUBROUTINE SAVLOD(LUNIT,ID,M,N,IMG,JOB,XREAL,XIMAG)
      INTEGER LUNIT,ID(4),M,N,IMG,JOB
      DOUBLE PRECISION XREAL(1),XIMAG(1)
C 
C     IMPLEMENT SAVE AND LOAD
C     LUNIT = LOGICAL UNIT NUMBER
C     ID = NAME, FORMAT 4A1
C     M, N = DIMENSIONS
C     IMG = NONZERO IF XIMAG IS NONZERO
C     JOB = 0     FOR SAVE
C         = SPACE AVAILABLE FOR LOAD
C     XREAL, XIMAG = REAL AND OPTIONAL IMAGINARY PARTS
C 
C     SYSTEM DEPENDENT FORMATS
  101 FORMAT(4A1,3I4)
  102 FORMAT(4Z18)
C 
      IF (JOB .GT. 0) GO TO 20
C 
C     SAVE
   10 WRITE(LUNIT,101) ID,M,N,IMG
      DO 15 J = 1, N
         K = (J-1)*M+1
         L = J*M
         WRITE(LUNIT,102) (XREAL(I),I=K,L)
         IF (IMG .NE. 0) WRITE(LUNIT,102) (XIMAG(I),I=K,L)
   15 CONTINUE
      RETURN
C 
C     LOAD
   20 READ(LUNIT,101,END=30) ID,M,N,IMG
      IF (M*N .GT. JOB) GO TO 30
      DO 25 J = 1, N
         K = (J-1)*M+1
         L = J*M
         READ(LUNIT,102,END=30) (XREAL(I),I=K,L)
         IF (IMG .NE. 0) READ(LUNIT,102,END=30) (XIMAG(I),I=K,L)
   25 CONTINUE
      RETURN
C 
C     END OF FILE
   30 M = 0
      N = 0
      RETURN
      END

 
      SUBROUTINE FORMZ(LUNIT,X,Y)
      DOUBLE PRECISION X,Y
C 
C     SYSTEM DEPENDENT ROUTINE TO PRINT WITH Z FORMAT
C 
      IF (Y .NE. 0.0D0) WRITE(LUNIT,10) X,Y
      IF (Y .EQ. 0.0D0) WRITE(LUNIT,10) X
   10 FORMAT(2Z18)
      RETURN
      END

 
      DOUBLE PRECISION FUNCTION FLOP(X)
      DOUBLE PRECISION X
C     SYSTEM DEPENDENT FUNCTION
C     COUNT AND POSSIBLY CHOP EACH FLOATING POINT OPERATION
C     FLP(1) IS FLOP COUNTER
C     FLP(2) IS NUMBER OF PLACES TO BE CHOPPED
C 
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
C 
      DOUBLE PRECISION MASK(14),XX,MM
      real mas(2,14)
      LOGICAL LX(2),LM(2)
      EQUIVALENCE (LX(1),XX),(LM(1),MM)
      equivalence (mask(1),mas(1))
      data mas/
     $ 'ffffffff'x,'fff0ffff'x,
     $ 'ffffffff'x,'ff00ffff'x,
     $ 'ffffffff'x,'f000ffff'x,
     $ 'ffffffff'x,'0000ffff'x,
     $ 'ffffffff'x,'0000fff0'x,
     $ 'ffffffff'x,'0000ff00'x,
     $ 'ffffffff'x,'0000f000'x,
     $ 'ffffffff'x,'00000000'x,
     $ 'fff0ffff'x,'00000000'x,
     $ 'ff00ffff'x,'00000000'x,
     $ 'f000ffff'x,'00000000'x,
     $ '0000ffff'x,'00000000'x,
     $ '0000fff0'x,'00000000'x,
     $ '0000ff80'x,'00000000'x/
C 
      FLP(1) = FLP(1) + 1
      K = FLP(2)
      FLOP = X
      IF (K .LE. 0) RETURN
      FLOP = 0.0D0
      IF (K .GE. 15) RETURN
      XX = X
      MM = MASK(K)
      LX(1) = LX(1) .AND. LM(1)
      LX(2) = LX(2) .AND. LM(2)
      FLOP = XX
      RETURN
      END

 
      SUBROUTINE XCHAR(BUF,K)
      INTEGER BUF(1),K
C 
C     SYSTEM DEPENDENT ROUTINE TO HANDLE SPECIAL CHARACTERS
C 
C 
      INTEGER BACK,MASK
      DATA BACK/'20202008'X/,MASK/'000000FF'X/
C 
      IF (BUF(1) .EQ. BACK) K = -1
      L = BUF(1) .AND. MASK
      IF (K .NE. -1) WRITE(6,10) BUF(1),L
   10 FORMAT(1X,1H',A1,4H' = ,Z2,' hex is not a MATLAB character.')
      RETURN
      END

 
      SUBROUTINE USER(A,M,N,S,T)
      DOUBLE PRECISION A(M,N),S,T
C 
      INTEGER A3(9)
      DATA A3 /-149,537,-27,-50,180,-9,-154,546,-25/
      IF (A(1,1) .NE. 3.0D0) RETURN
      DO 10 I = 1, 9
         A(I,1) = DFLOAT(A3(I))
   10 CONTINUE
      M = 3
      N = 3
      RETURN
      END

 
      SUBROUTINE PROMPT(PAUSE)
      INTEGER PAUSE
C 
C     ISSUE MATLAB PROMPT WITH OPTIONAL PAUSE
C 
      INTEGER DDT,ERR,FMT,LCT(4),LIN(1024),LPT(6),HIO,RIO,WIO,RTE,WTE
      COMMON /IOP/ DDT,ERR,FMT,LCT,LIN,LPT,HIO,RIO,WIO,RTE,WTE
      WRITE(WTE,10)
      IF (WIO .NE. 0) WRITE(WIO,10)
   10 FORMAT(/1X,'<>',$)
      IF (PAUSE .EQ. 1) READ(RTE,20) DUMMY
   20 FORMAT(A1)
      RETURN
      END

 
      SUBROUTINE PLOT(LUNIT,X,Y,N,P,K,BUF)
      DOUBLE PRECISION X(N),Y(N),P(1)
      INTEGER BUF(79)
C 
C     PLOT X VS. Y ON LUNIT
C     IF K IS NONZERO, THEN P(1),...,P(K) ARE EXTRA PARAMETERS
C     BUF IS WORK SPACE
C 
      DOUBLE PRECISION XMIN,YMIN,XMAX,YMAX,DY,DX,Y1,Y0
      INTEGER AST,BLANK,H,W
      DATA AST/1H*/,BLANK/1H /,H/20/,W/79/
C 
C     H = HEIGHT, W = WIDTH
C 
      IF (K .GT. 0) WRITE(LUNIT,01) (P(I), I=1,K)
   01 FORMAT('Extra parameters',10f5.1)
      XMIN = X(1)
      XMAX = X(1)
      YMIN = Y(1)
      YMAX = Y(1)
      DO 10 I = 1, N
         XMIN = DMIN1(XMIN,X(I))
         XMAX = DMAX1(XMAX,X(I))
         YMIN = DMIN1(YMIN,Y(I))
         YMAX = DMAX1(YMAX,Y(I))
   10 CONTINUE
      DX = XMAX - XMIN
      IF (DX .EQ. 0.0D0) DX = 1.0D0
      DY = YMAX - YMIN
      WRITE(LUNIT,35)
      DO 40 L = 1, H
         DO 20 J = 1, W
            BUF(J) = BLANK
   20    CONTINUE
         Y1 = YMIN + (H-L+1)*DY/H
         Y0 = YMIN + (H-L)*DY/H
         JMAX = 1
         DO 30 I = 1, N
            IF (Y(I) .GT. Y1) GO TO 30
            IF (L.NE.H .AND. Y(I).LE.Y0) GO TO 30
            J = 1 + (W-1)*(X(I) - XMIN)/DX
            BUF(J) = AST
            JMAX = MAX0(JMAX,J)
   30    CONTINUE
         WRITE(LUNIT,35) (BUF(J),J=1,JMAX)
   35    FORMAT(79A1)
   40 CONTINUE
      RETURN
      END

 
      SUBROUTINE EDIT(BUF,N)
      INTEGER BUF(N)
C 
C     CALLED AFTER INPUT OF A SINGLE BACKSLASH
C     BUF CONTAINS PREVIOUS INPUT LINE, ONE CHAR PER WORD
C     ENTER LOCAL EDITOR IF AVAILABLE
C     OTHERWISE JUST
      RETURN
      END