      SUBROUTINE MC01SD( DP, P, S, T, MANT, E, IWORK, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To scale the coefficients of the real polynomial P(x) such that
C     the coefficients of the scaled polynomial Q(x) = sP(tx) have
C     minimal variation, where s and t are real scalars.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     DP      (input) INTEGER
C             The degree of the polynomial P(x).  DP >= 0.
C
C     P       (input/output) DOUBLE PRECISION array, dimension (DP+1)
C             On entry, this array must contain the coefficients of P(x)
C             in increasing powers of x.
C             On exit, this array contains the coefficients of the
C             scaled polynomial Q(x) in increasing powers of x.
C
C     S       (output) INTEGER
C             The exponent of the floating-point representation of the
C             scaling factor s = BASE**S, where BASE is the base of the
C             machine representation of floating-point numbers (see
C             LAPACK Library Routine DLAMCH).
C
C     T       (output) INTEGER
C             The exponent of the floating-point representation of the
C             scaling factor t = BASE**T.
C
C     MANT    (output) DOUBLE PRECISION array, dimension (DP+1)
C             This array contains the mantissas of the standard
C             floating-point representation of the coefficients of the
C             scaled polynomial Q(x) in increasing powers of x.
C
C     E       (output) INTEGER array, dimension (DP+1)
C             This array contains the exponents of the standard
C             floating-point representation of the coefficients of the
C             scaled polynomial Q(x) in increasing powers of x.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (DP+1)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if on entry, P(x) is the zero polynomial.
C
C     METHOD
C
C     Define the variation of the coefficients of the real polynomial
C
C                                         2                DP
C        P(x) = p(0) + p(1) * x + p(2) * x  + ... + p(DP) x
C
C     whose non-zero coefficients can be represented as
C                          e(i)
C        p(i) = m(i) * BASE     (where 1 <= ABS(m(i)) < BASE)
C
C     by
C
C        V = max(e(i)) - min(e(i)),
C
C     where max and min are taken over the indices i for which p(i) is
C     non-zero.
C                                        DP         i    i
C     For the scaled polynomial P(cx) = SUM p(i) * c  * x  with
C                                       i=0
C                j
C     c  = (BASE) , the variation V(j) is given by
C
C       V(j) = max(e(i) + j * i) - min(e(i) + j * i).
C
C     Using the fact that V(j) is a convex function of j, the routine
C     determines scaling factors s = (BASE)**S and t = (BASE)**T such
C     that the coefficients of the scaled polynomial Q(x) = sP(tx)
C     satisfy the following conditions:
C
C       (a) 1 <= q(0) < BASE and
C
C       (b) the variation of the coefficients of Q(x) is minimal.
C
C     Further details can be found in [1].
C
C     REFERENCES
C
C     [1] Dunaway, D.K.
C         Calculation of Zeros of a Real Polynomial through
C         Factorization using Euclid's Algorithm.
C         SIAM J. Numer. Anal., 11, pp. 1087-1104, 1974.
C
C     NUMERICAL ASPECTS
C
C     Since the scaling is performed on the exponents of the floating-
C     point representation of the coefficients of P(x), no rounding
C     errors occur during the computation of the coefficients of Q(x).
C
C     FURTHER COMMENTS
C
C     The scaling factors s and t are BASE dependent.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01GD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, polynomial operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DP, INFO, S, T
C     .. Array Arguments ..
      INTEGER           E(*), IWORK(*)
      DOUBLE PRECISION  MANT(*), P(*)
C     .. Local Scalars ..
      LOGICAL           OVFLOW
      INTEGER           BETA, DV, I, INC, J, LB, M, UB, V0, V1
C     .. External Functions ..
      INTEGER           MC01SX
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, MC01SX
C     .. External Subroutines ..
      EXTERNAL          MC01SW, MC01SY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, NINT
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      IF( DP.LT.0 ) THEN
         INFO = -1
C
C        Error return.
C
         CALL XERBLA( 'MC01SD', -INFO )
         RETURN
      END IF
C
      INFO = 0
      LB = 1
C     WHILE ( LB <= DP+1 and P(LB) = 0 ) DO
   20 IF ( LB.LE.DP+1 ) THEN
         IF ( P(LB).EQ.ZERO ) THEN
            LB = LB + 1
            GO TO 20
         END IF
      END IF
C     END WHILE 20
C
C     LB = MIN( i: P(i) non-zero).
C
      IF ( LB.EQ.DP+2 ) THEN
         INFO = 1
         RETURN
      END IF
C
      UB = DP + 1
C     WHILE ( P(UB) = 0 ) DO
   40 IF ( P(UB).EQ.ZERO ) THEN
         UB = UB - 1
         GO TO 40
      END IF
C     END WHILE 40
C
C     UB = MAX(i: P(i) non-zero).
C
      BETA = DLAMCH( 'Base' )
C
      DO 60 I = 1, DP + 1
         CALL MC01SW( P(I), BETA, MANT(I), E(I) )
   60 CONTINUE
C
C     First prescaling.
C
      M = E(LB)
      IF ( M.NE.0 ) THEN
C
         DO 80 I = LB, UB
            IF ( MANT(I).NE.ZERO ) E(I) = E(I) - M
   80    CONTINUE
C
      END IF
      S = -M
C
C     Second prescaling.
C
      IF ( UB.GT.1 ) M = NINT( DBLE( E(UB) )/DBLE( UB-1 ) )
C
      DO 100 I = LB, UB
         IF ( MANT(I).NE.ZERO ) E(I) = E(I) - M*(I-1)
  100 CONTINUE
C
      T = -M
C
      V0 = MC01SX( LB, UB, E, MANT )
      J = 1
C
      DO 120 I = LB, UB
         IF ( MANT(I).NE.ZERO ) IWORK(I) = E(I) + (I-1)
  120 CONTINUE
C
      V1 = MC01SX( LB, UB, IWORK, MANT )
      DV = V1 - V0
      IF ( DV.NE.0 ) THEN
         IF ( DV.GT.0 ) THEN
            J = 0
            INC = -1
            V1 = V0
            DV = -DV
C
            DO 130 I = LB, UB
               IWORK(I) = E(I)
  130       CONTINUE
C
         ELSE
            INC = 1
         END IF
C        WHILE ( DV < 0 ) DO
  140    IF ( DV.LT.0 ) THEN
            V0 = V1
C
            DO 150 I = LB, UB
               E(I) = IWORK(I)
  150       CONTINUE
C
            J = J + INC
C
            DO 160 I = LB, UB
               IWORK(I) = E(I) + INC*(I-1 )
  160       CONTINUE
C
            V1 = MC01SX( LB, UB, IWORK, MANT )
            DV = V1 - V0
            GO TO 140
         END IF
C        END WHILE 140
         T = T + J - INC
      END IF
C
C     Evaluation of the output parameters.
C
      DO 180 I = LB, UB
         CALL MC01SY( MANT(I), E(I), BETA, P(I), OVFLOW )
  180 CONTINUE
C
      RETURN
C *** Last line of MC01SD ***
      END
