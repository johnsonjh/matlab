      SUBROUTINE WQRSL(XR,XI,LDX,N,K,QRAUXR,QRAUXI,YR,YI,QYR,QYI,QTYR,
     *                 QTYI,BR,BI,RSDR,RSDI,XBR,XBI,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      DOUBLE PRECISION XR(LDX,1),XI(LDX,1),QRAUXR(1),QRAUXI(1),YR(1),
     *                 YI(1),QYR(1),QYI(1),QTYR(1),QTYI(1),BR(1),BI(1),
     *                 RSDR(1),RSDI(1),XBR(1),XBI(1)
C 
C     WQRSL APPLIES THE OUTPUT OF WQRDC TO COMPUTE COORDINATE
C     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
C     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
C 
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C 
C     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
C     N X P MATRIX X THAT WAS INPUT TO WQRDC (IF NO PIVOTING WAS
C     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
C     ORIGINAL ORDER).  WQRDC PRODUCES A FACTORED UNITARY MATRIX Q
C     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
C 
C              XK = Q * (R)
C                       (0)
C 
C     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
C     X AND QRAUX.
C 
C     ON ENTRY
C 
C        X      DOUBLE-COMPLEX(LDX,P).
C               X CONTAINS THE OUTPUT OF WQRDC.
C 
C        LDX    INTEGER.
C               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C 
C        N      INTEGER.
C               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
C               HAVE THE SAME VALUE AS N IN WQRDC.
C 
C        K      INTEGER.
C               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
C               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
C               SAME AS IN THE CALLING SEQUENCE TO WQRDC.
C 
C        QRAUX  DOUBLE-COMPLEX(P).
C               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM WQRDC.
C 
C        Y      DOUBLE-COMPLEX(N)
C               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
C               BY WQRSL.
C 
C        JOB    INTEGER.
C               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
C               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
C               MEANING.
C 
C                    IF A.NE.0, COMPUTE QY.
C                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
C                    IF C.NE.0, COMPUTE B.
C                    IF D.NE.0, COMPUTE RSD.
C                    IF E.NE.0, COMPUTE XB.
C 
C               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
C               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
C               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
C               SEQUENCE.
C 
C     ON RETURN
C 
C        QY     DOUBLE-COMPLEX(N).
C               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
C               REQUESTED.
C 
C        QTY    DOUBLE-COMPLEX(N).
C               QTY CONTAINS CTRANS(Q)*Y, IF ITS COMPUTATION HAS
C               BEEN REQUESTED.  HERE CTRANS(Q) IS THE CONJUGATE
C               TRANSPOSE OF THE MATRIX Q.
C 
C        B      DOUBLE-COMPLEX(K)
C               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
C 
C                    MINIMIZE NORM2(Y - XK*B),
C 
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
C               IF PIVOTING WAS REQUESTED IN WQRDC, THE J-TH
C               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
C               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO WQRDC.)
C 
C        RSD    DOUBLE-COMPLEX(N).
C               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
C               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
C               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
C 
C        XB     DOUBLE-COMPLEX(N).
C               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
C               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
C               OF X.
C 
C        INFO   INTEGER.
C               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
C               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
C               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
C               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
C 
C     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
C     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
C     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
C     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
C     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
C     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
C     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
C     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
C     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
C     COMPUTED.  THUS THE CALLING SEQUENCE
C 
C          CALL WQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C 
C     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
C     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
C     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
C     A SINGLE CALLINNG SEQUENCE.
C 
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C 
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C 
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C 
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C 
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C 
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C 
C     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
C     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
C 
C     LINPACK. THIS VERSION DATED 07/03/79 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C 
C     WQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C 
C     BLAS WAXPY,WCOPY,WDOTCR,WDOTCI
C     FORTRAN DABS,DIMAG,MIN0,MOD
C 
C     INTERNAL VARIABLES
C 
      INTEGER I,J,JJ,JU,KP1
      DOUBLE PRECISION WDOTCR,WDOTCI,TR,TI,TEMPR,TEMPI
      LOGICAL CB,CQY,CQTY,CR,CXB
C 
      DOUBLE PRECISION ZDUMR,ZDUMI
      DOUBLE PRECISION CABS1
      CABS1(ZDUMR,ZDUMI) = DABS(ZDUMR) + DABS(ZDUMI)
C 
C     SET INFO FLAG.
C 
      INFO = 0
C 
C     DETERMINE WHAT IS TO BE COMPUTED.
C 
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C 
C     SPECIAL ACTION WHEN N=1.
C 
      IF (JU .NE. 0) GO TO 80
         IF (.NOT.CQY) GO TO 10
            QYR(1) = YR(1)
            QYI(1) = YI(1)
   10    CONTINUE
         IF (.NOT.CQTY) GO TO 20
            QTYR(1) = YR(1)
            QTYI(1) = YI(1)
   20    CONTINUE
         IF (.NOT.CXB) GO TO 30
            XBR(1) = YR(1)
            XBI(1) = YI(1)
   30    CONTINUE
         IF (.NOT.CB) GO TO 60
            IF (CABS1(XR(1,1),XI(1,1)) .NE. 0.0D0) GO TO 40
               INFO = 1
            GO TO 50
   40       CONTINUE
               CALL WDIV(YR(1),YI(1),XR(1,1),XI(1,1),BR(1),BI(1))
   50       CONTINUE
   60    CONTINUE
         IF (.NOT.CR) GO TO 70
            RSDR(1) = 0.0D0
            RSDI(1) = 0.0D0
   70    CONTINUE
      GO TO 290
   80 CONTINUE
C 
C        SET UP TO COMPUTE QY OR QTY.
C 
         IF (CQY) CALL WCOPY(N,YR,YI,1,QYR,QYI,1)
         IF (CQTY) CALL WCOPY(N,YR,YI,1,QTYR,QTYI,1)
         IF (.NOT.CQY) GO TO 110
C 
C           COMPUTE QY.
C 
            DO 100 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUXR(J),QRAUXI(J)) .EQ. 0.0D0)
     *            GO TO 90
                  TEMPR = XR(J,J)
                  TEMPI = XI(J,J)
                  XR(J,J) = QRAUXR(J)
                  XI(J,J) = QRAUXI(J)
                  TR = -WDOTCR(N-J+1,XR(J,J),XI(J,J),1,QYR(J),QYI(J),1)
                  TI = -WDOTCI(N-J+1,XR(J,J),XI(J,J),1,QYR(J),QYI(J),1)
                  CALL WDIV(TR,TI,XR(J,J),XI(J,J),TR,TI)
                  CALL WAXPY(N-J+1,TR,TI,XR(J,J),XI(J,J),1,QYR(J),
     *                       QYI(J),1)
                  XR(J,J) = TEMPR
                  XI(J,J) = TEMPI
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
         IF (.NOT.CQTY) GO TO 140
C 
C           COMPUTE CTRANS(Q)*Y.
C 
            DO 130 J = 1, JU
               IF (CABS1(QRAUXR(J),QRAUXI(J)) .EQ. 0.0D0)
     *            GO TO 120
                  TEMPR = XR(J,J)
                  TEMPI = XI(J,J)
                  XR(J,J) = QRAUXR(J)
                  XI(J,J) = QRAUXI(J)
                  TR = -WDOTCR(N-J+1,XR(J,J),XI(J,J),1,QTYR(J),
     *                         QTYI(J),1)
                  TI = -WDOTCI(N-J+1,XR(J,J),XI(J,J),1,QTYR(J),
     *                         QTYI(J),1)
                  CALL WDIV(TR,TI,XR(J,J),XI(J,J),TR,TI)
                  CALL WAXPY(N-J+1,TR,TI,XR(J,J),XI(J,J),1,QTYR(J),
     *                       QTYI(J),1)
                  XR(J,J) = TEMPR
                  XI(J,J) = TEMPI
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
C 
C        SET UP TO COMPUTE B, RSD, OR XB.
C 
         IF (CB) CALL WCOPY(K,QTYR,QTYI,1,BR,BI,1)
         KP1 = K + 1
         IF (CXB) CALL WCOPY(K,QTYR,QTYI,1,XBR,XBI,1)
         IF (CR .AND. K .LT. N)
     *      CALL WCOPY(N-K,QTYR(KP1),QTYI(KP1),1,RSDR(KP1),RSDI(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 160
            DO 150 I = KP1, N
               XBR(I) = 0.0D0
               XBI(I) = 0.0D0
  150       CONTINUE
  160    CONTINUE
         IF (.NOT.CR) GO TO 180
            DO 170 I = 1, K
               RSDR(I) = 0.0D0
               RSDI(I) = 0.0D0
  170       CONTINUE
  180    CONTINUE
         IF (.NOT.CB) GO TO 230
C 
C           COMPUTE B.
C 
            DO 210 JJ = 1, K
               J = K - JJ + 1
               IF (CABS1(XR(J,J),XI(J,J)) .NE. 0.0D0) GO TO 190
                  INFO = J
C                 ......EXIT
C           ......EXIT
                  GO TO 220
  190          CONTINUE
               CALL WDIV(BR(J),BI(J),XR(J,J),XI(J,J),BR(J),BI(J))
               IF (J .EQ. 1) GO TO 200
                  TR = -BR(J)
                  TI = -BI(J)
                  CALL WAXPY(J-1,TR,TI,XR(1,J),XI(1,J),1,BR,BI,1)
  200          CONTINUE
  210       CONTINUE
  220       CONTINUE
  230    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 280
C 
C           COMPUTE RSD OR XB AS REQUIRED.
C 
            DO 270 JJ = 1, JU
               J = JU - JJ + 1
               IF (CABS1(QRAUXR(J),QRAUXI(J)) .EQ. 0.0D0)
     *            GO TO 260
                  TEMPR = XR(J,J)
                  TEMPI = XI(J,J)
                  XR(J,J) = QRAUXR(J)
                  XI(J,J) = QRAUXI(J)
                  IF (.NOT.CR) GO TO 240
                     TR = -WDOTCR(N-J+1,XR(J,J),XI(J,J),1,RSDR(J),
     *                            RSDI(J),1)
                     TI = -WDOTCI(N-J+1,XR(J,J),XI(J,J),1,RSDR(J),
     *                            RSDI(J),1)
                     CALL WDIV(TR,TI,XR(J,J),XI(J,J),TR,TI)
                     CALL WAXPY(N-J+1,TR,TI,XR(J,J),XI(J,J),1,RSDR(J),
     *                          RSDI(J),1)
  240             CONTINUE
                  IF (.NOT.CXB) GO TO 250
                     TR = -WDOTCR(N-J+1,XR(J,J),XI(J,J),1,XBR(J),
     *                            XBI(J),1)
                     TI = -WDOTCI(N-J+1,XR(J,J),XI(J,J),1,XBR(J),
     *                            XBI(J),1)
                     CALL WDIV(TR,TI,XR(J,J),XI(J,J),TR,TI)
                     CALL WAXPY(N-J+1,TR,TI,XR(J,J),XI(J,J),1,XBR(J),
     *                          XBI(J),1)
  250             CONTINUE
                  XR(J,J) = TEMPR
                  XI(J,J) = TEMPI
  260          CONTINUE
  270       CONTINUE
  280    CONTINUE
  290 CONTINUE
      RETURN
      END
