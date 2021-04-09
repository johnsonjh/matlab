      PROGRAM MAIN
C
C     PROGRAM MAIN FOR MS-DOS
C
      CALL MATLAB(0)
      STOP
      END

      SUBROUTINE FILES(LUNIT,NAME)
      INTEGER LUNIT,NAME(32)
C
C     SYSTEM DEPENDENT ROUTINE TO ALLOCATE FILES
C     LUNIT = LOGICAL UNIT NUMBER
C     NAME = FILE NAME, 1 CHARACTER PER WORD
C
C
      logical fexist
      integer temp
      character*1 namec(32)
      character*32 names
      equivalence (names, namec)
C
      if (lunit .eq. 0) then
         return
      end if
      if (lunit .eq. 9) then
         open(9,file='matlab.hlp',status='old')
         write(*,09)
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
C
      DO 10 I = 1, 32
        TEMP = NAME(I)
        NAMEC(I) = CHAR(TEMP)
   10 CONTINUE
C
      inquire(file=names,exist=fexist)
      if (fexist) then
        open(lunit,file=names,status='old')
       else
        open(lunit,file=names,status='new')
      end if
C
   20 RETURN
      END

      SUBROUTINE SAVLOD(LUNIT,ID,M,N,IMG,JOB,XREAL,XIMAG)
      CHARACTER*1 IDC(4)
      INTEGER LUNIT,ID(4),M,N,IMG,JOB,TEMP
      DOUBLE PRECISION XREAL(1),XIMAG(1)
C
C     IMPLEMENT SAVE AND LOAD
C     LUNIT = LOGICAL UNIT NUMBER
C     ID = NAME, FORMAT 4I1
C     IDC = NAME, FORMAT 4A1
C     M, N = DIMENSIONS
C     IMG = NONZERO IF XIMAG IS NONZERO
C     JOB = 0     FOR SAVE
C         = SPACE AVAILABLE FOR LOAD
C     XREAL, XIMAG = REAL AND OPTIONAL IMAGINARY PARTS
C
C     SYSTEM DEPENDENT FORMATS
  101 FORMAT(4A1,3I4)
  102 FORMAT(G20.8E3)
C
      IF (JOB .GT. 0) GO TO 20
C
C     SAVE
C
C     CONVERT ID TO CHARACTER STRING
C
      DO 5 I = 1, 4
        TEMP = ID(I)
        IDC(I) = CHAR(TEMP)
    5 CONTINUE
C
   10 WRITE(LUNIT,101) IDC,M,N,IMG
      DO 15 J = 1, N
         K = (J-1)*M+1
         L = J*M
         WRITE(LUNIT,102) (XREAL(I),I=K,L)
         IF (IMG .NE. 0) WRITE(LUNIT,102) (XIMAG(I),I=K,L)
   15 CONTINUE
      RETURN
C
C     LOAD
   20 READ(LUNIT,101,END=30) IDC,M,N,IMG
C
C     CONVERT IDC TO INTEGER STRING
C
      DO 21 I = 1, 4
         ID(I) = ICHAR(IDC(I))
   21 CONTINUE
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
   10 FORMAT(2G18.10E3)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FLOP(X)
      DOUBLE PRECISION X
C     SYSTEM DEPENDENT FUNCTION
C     COUNT EACH FLOATING POINT OPERATION
C     FLP(1) IS FLOP COUNTER
C     FLP(2) IS NUMBER OF PLACES TO BE CHOPPED
C
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN
C
C
      FLP(1) = FLP(1) + 1
      K = FLP(2)
      FLOP = X
      RETURN
      END

      SUBROUTINE XCHAR(BUF,K)
      INTEGER BUF(1),K
C
C     SYSTEM DEPENDENT ROUTINE TO HANDLE SPECIAL CHARACTERS
C
C
      RETURN
      END

      SUBROUTINE USER(A,M,N,S,T)
      DOUBLE PRECISION A(M,N),S,T
C
      INTEGER A3(9)
      DATA A3 /-149,537,-27,-50,180,-9,-154,546,-25/
      IF (A(1,1) .NE. 3.0D0) RETURN
      DO 10 I = 1, 9
         A(I,1) = DBLE(A3(I))
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
      CHARACTER ASTC,BLANKC
      INTEGER AST,BLANK,H,W
      DATA ASTC/'*'/,BLANKC/' '/,H/20/,W/79/
C
      AST=ICHAR(ASTC)
      BLANK=ICHAR(BLANKC)
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
C
      RETURN
      END