      DOUBLE PRECISION FUNCTION ROUNDM(X)
      DOUBLE PRECISION X,Y,Z,E,H
      DATA H/1.0D9/
      Z = DABS(X)
      Y = Z + 1.0D0
      IF (Y .EQ. Z) GO TO 40
      Y = 0.0D0
      E = H
   10 IF (E .GE. Z) GO TO 20
         E = 2.0D0*E
         GO TO 10
   20 IF (E .LE. H) GO TO 30
         IF (E .LE. Z) Y = Y + E
         IF (E .LE. Z) Z = Z - E
         E = E/2.0D0
         GO TO 20
   30 Z = IDINT(Z + 0.5D0)
      Y = Y + Z
      IF (X .LT. 0.0D0) Y = -Y
      ROUNDM = Y
      RETURN
   40 ROUNDM = X
      RETURN
      END
