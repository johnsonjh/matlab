      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C 
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      DOUBLE PRECISION AR(NM,N),AI(NM,N),ORTR(IGH),ORTI(IGH)
      DOUBLE PRECISION F,G,H,FI,FR,SCALE
      DOUBLE PRECISION FLOP,PYTHAG
C 
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C 
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
C 
C     ON INPUT.
C 
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C 
C        N IS THE ORDER OF THE MATRIX.
C 
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C 
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C 
C     ON OUTPUT.
C 
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C 
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C 
C     ------------------------------------------------------------------
C 
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C 
      DO 180 M = KP1, LA
         H = 0.0D0
         ORTR(M) = 0.0D0
         ORTI(M) = 0.0D0
         SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = FLOP(SCALE + DABS(AR(I,M-1)) + DABS(AI(I,M-1)))
C 
         IF (SCALE .EQ. 0.0D0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORTR(I) = FLOP(AR(I,M-1)/SCALE)
            ORTI(I) = FLOP(AI(I,M-1)/SCALE)
            H = FLOP(H + ORTR(I)*ORTR(I) + ORTI(I)*ORTI(I))
  100    CONTINUE
C 
         G = FLOP(DSQRT(H))
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0D0) GO TO 103
         H = FLOP(H + F*G)
         G = FLOP(G/F)
         ORTR(M) = FLOP((1.0D0 + G)*ORTR(M))
         ORTI(M) = FLOP((1.0D0 + G)*ORTI(M))
         GO TO 105
C 
  103    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H)*A ..........
  105    DO 130 J = M, N
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               FR = FLOP(FR + ORTR(I)*AR(I,J) + ORTI(I)*AI(I,J))
               FI = FLOP(FI + ORTR(I)*AI(I,J) - ORTI(I)*AR(I,J))
  110       CONTINUE
C 
            FR = FLOP(FR/H)
            FI = FLOP(FI/H)
C 
            DO 120 I = M, IGH
               AR(I,J) = FLOP(AR(I,J) - FR*ORTR(I) + FI*ORTI(I))
               AI(I,J) = FLOP(AI(I,J) - FR*ORTI(I) - FI*ORTR(I))
  120       CONTINUE
C 
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               FR = FLOP(FR + ORTR(J)*AR(I,J) - ORTI(J)*AI(I,J))
               FI = FLOP(FI + ORTR(J)*AI(I,J) + ORTI(J)*AR(I,J))
  140       CONTINUE
C 
            FR = FLOP(FR/H)
            FI = FLOP(FI/H)
C 
            DO 150 J = M, IGH
               AR(I,J) = FLOP(AR(I,J) - FR*ORTR(J) - FI*ORTI(J))
               AI(I,J) = FLOP(AI(I,J) + FR*ORTI(J) - FI*ORTR(J))
  150       CONTINUE
C 
  160    CONTINUE
C 
         ORTR(M) = FLOP(SCALE*ORTR(M))
         ORTI(M) = FLOP(SCALE*ORTI(M))
         AR(M,M-1) = FLOP(-G*AR(M,M-1))
         AI(M,M-1) = FLOP(-G*AI(M,M-1))
  180 CONTINUE
C 
  200 RETURN
      END
