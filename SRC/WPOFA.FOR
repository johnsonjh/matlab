      SUBROUTINE WPOFA(AR,AI,LDA,N,INFO)
      DOUBLE PRECISION AR(LDA,1),AI(LDA,1)
      DOUBLE PRECISION S,TR,TI,WDOTCR,WDOTCI
      DO 30 J = 1, N
         INFO = J
         S = 0.0D0
         JM1 = J-1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 K = 1, JM1
            TR = AR(K,J)-WDOTCR(K-1,AR(1,K),AI(1,K),1,AR(1,J),AI(1,J),1)
            TI = AI(K,J)-WDOTCI(K-1,AR(1,K),AI(1,K),1,AR(1,J),AI(1,J),1)
            CALL WDIV(TR,TI,AR(K,K),AI(K,K),TR,TI)
            AR(K,J) = TR
            AI(K,J) = TI
            S = S + TR*TR + TI*TI
   10    CONTINUE
   20    CONTINUE
         S = AR(J,J) - S
         IF (S.LE.0.0D0 .OR. AI(J,J).NE.0.0D0) GO TO 40
         AR(J,J) = DSQRT(S)
   30 CONTINUE
      INFO = 0
   40 RETURN
      END
