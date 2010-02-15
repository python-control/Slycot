      SUBROUTINE MC01TD( DICO, DP, P, STABLE, NZ, DWORK, IWARN, INFO )
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
C     To determine whether or not a given polynomial P(x) with real
C     coefficients is stable, either in the continuous-time or discrete-
C     time case.
C
C     A polynomial is said to be stable in the continuous-time case
C     if all its zeros lie in the left half-plane, and stable in the
C     discrete-time case if all its zeros lie inside the unit circle.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Indicates whether the stability test to be applied to
C             P(x) is in the continuous-time or discrete-time case as
C             follows:
C             = 'C':  Continuous-time case;
C             = 'D':  Discrete-time case.
C
C     Input/Output Parameters
C
C     DP      (input/output) INTEGER
C             On entry, the degree of the polynomial P(x).  DP >= 0.
C             On exit, if P(DP+1) = 0.0 on entry, then DP contains the
C             index of the highest power of x for which P(DP+1) <> 0.0.
C
C     P       (input) DOUBLE PRECISION array, dimension (DP+1)
C             This array must contain the coefficients of P(x) in
C             increasing powers of x.
C
C     STABLE  (output) LOGICAL
C             Contains the value .TRUE. if P(x) is stable and the value
C             .FALSE. otherwise (see also NUMERICAL ASPECTS).
C
C     NZ      (output) INTEGER
C             If INFO = 0, contains the number of unstable zeros - that
C             is, the number of zeros of P(x) in the right half-plane if
C             DICO = 'C' or the number of zeros of P(x) outside the unit
C             circle if DICO = 'D' (see also NUMERICAL ASPECTS).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (2*DP+2)
C             The leading (DP+1) elements of DWORK contain the Routh
C             coefficients, if DICO = 'C', or the constant terms of
C             the Schur-Cohn transforms, if DICO = 'D'.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = k:  if the degree of the polynomial P(x) has been
C                   reduced to (DB - k) because P(DB+1-j) = 0.0 on entry
C                   for j = 0, 1,..., k-1 and P(DB+1-k) <> 0.0.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if on entry, P(x) is the zero polynomial;
C             = 2:  if the polynomial P(x) is most probably unstable,
C                   although it may be stable with one or more zeros
C                   very close to either the imaginary axis if
C                   DICO = 'C' or the unit circle if DICO = 'D'.
C                   The number of unstable zeros (NZ) is not determined.
C
C     METHOD
C
C     The stability of the real polynomial
C                                         2                DP
C        P(x) = p(0) + p(1) * x + p(2) * x  + ... + p(DP) x
C
C     is determined as follows.
C
C     In the continuous-time case (DICO = 'C') the Routh algorithm
C     (see [1]) is used. The routine computes the Routh coefficients and
C     if they are non-zero then the number of sign changes in the
C     sequence of the coefficients is equal to the number of zeros with
C     positive imaginary part.
C
C     In the discrete-time case (DICO = 'D') the Schur-Cohn
C     algorithm (see [2] and [3]) is applied to the reciprocal
C     polynomial
C                                                2               DP
C        Q(x) = p(DP) + p(DP-1) * x + p(DP-2) * x  + ... + p(0) x  .
C
C     The routine computes the constant terms of the Schur transforms
C     and if all of them are non-zero then the number of zeros of P(x)
C     with modulus greater than unity is obtained from the sequence of
C     constant terms.
C
C     REFERENCES
C
C     [1] Gantmacher, F.R.
C         Applications of the Theory of Matrices.
C         Interscience Publishers, New York, 1959.
C
C     [2] Kucera, V.
C         Discrete Linear Control. The Algorithmic Approach.
C         John Wiley & Sons, Chichester, 1979.
C
C     [3] Henrici, P.
C         Applied and Computational Complex Analysis (Vol. 1).
C         John Wiley & Sons, New York, 1974.
C
C     NUMERICAL ASPECTS
C
C     The algorithm used by the routine is numerically stable.
C
C     Note that if some of the Routh coefficients (DICO = 'C') or
C     some of the constant terms of the Schur-Cohn transforms (DICO =
C     'D') are small relative to EPS (the machine precision), then
C     the number of unstable zeros (and hence the value of STABLE) may
C     be incorrect.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01HD by F. Delebecque and
C     A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, polynomial operations,
C     stability, stability criteria, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO
      LOGICAL           STABLE
      INTEGER           DP, INFO, IWARN, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), P(*)
C     .. Local Scalars ..
      LOGICAL           DICOC
      INTEGER           I, K, K1, K2, SIGNUM
      DOUBLE PRECISION  ALPHA, P1, PK1
C     .. External Functions ..
      INTEGER           IDAMAX
      LOGICAL           LSAME
      EXTERNAL          IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DRSCL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
C
      IWARN = 0
      INFO  = 0
      DICOC = LSAME( DICO, 'C' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.DICOC .AND. .NOT.LSAME( DICO, 'D' ) ) THEN
         INFO = -1
      ELSE IF( DP.LT.0 ) THEN
         INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC01TD', -INFO )
         RETURN
      END IF
C
C     WHILE (DP >= 0 and P(DP+1) = 0 ) DO
   20 IF ( DP.GE.0 ) THEN
         IF ( P(DP+1).EQ.ZERO ) THEN
            DP = DP - 1
            IWARN = IWARN + 1
            GO TO 20
         END IF
      END IF
C     END WHILE 20
C
      IF ( DP.EQ.-1 ) THEN
         INFO = 1
         RETURN
      END IF
C
C     P(x) is not the zero polynomial and its degree is exactly DP.
C
      IF ( DICOC ) THEN
C
C        Continuous-time case.
C
C        Compute the Routh coefficients and the number of sign changes.
C
         CALL DCOPY( DP+1, P, 1, DWORK, 1 )
         NZ = 0
         K = DP
C        WHILE ( K > 0 and DWORK(K) non-zero) DO
   40    IF ( K.GT.0 ) THEN
            IF ( DWORK(K).EQ.ZERO ) THEN
               INFO = 2
            ELSE
               ALPHA = DWORK(K+1)/DWORK(K)
               IF ( ALPHA.LT.ZERO ) NZ = NZ + 1
               K = K - 1
C
               DO 60 I = K, 2, -2
                  DWORK(I) = DWORK(I) - ALPHA*DWORK(I-1)
   60          CONTINUE
C
               GO TO 40
            END IF
         END IF
C        END WHILE 40
      ELSE
C
C        Discrete-time case.
C
C        To apply [3], section 6.8, on the reciprocal of polynomial
C        P(x) the elements of the array P are copied in DWORK in
C        reverse order.
C
         CALL DCOPY( DP+1, P, 1, DWORK, -1 )
C                                                           K-1
C        DWORK(K),...,DWORK(DP+1), are the coefficients of T   P(x)
C        scaled with a factor alpha(K) in order to avoid over- or
C        underflow,
C                                                    i-1
C        DWORK(i), i = 1,...,K, contains alpha(i) * T   P(0).
C
         SIGNUM = ONE
         NZ = 0
         K  = 1
C        WHILE ( K <= DP and DWORK(K) non-zero ) DO
   80    IF ( ( K.LE.DP ) .AND. ( INFO.EQ.0 ) ) THEN
C                                        K
C           Compute the coefficients of T P(x).
C
            K1 = DP - K + 2
            K2 = DP + 2
            ALPHA = DWORK(K-1+IDAMAX( K1, DWORK(K), 1 ))
            IF ( ALPHA.EQ.ZERO ) THEN
               INFO = 2
            ELSE
               CALL DCOPY( K1, DWORK(K), 1, DWORK(K2), 1 )
               CALL DRSCL( K1, ALPHA, DWORK(K2), 1 )
               P1  = DWORK(K2)
               PK1 = DWORK(K2+K1-1)
C
               DO 100 I = 1, K1 - 1
                  DWORK(K+I) = P1*DWORK(DP+1+I) - PK1*DWORK(K2+K1-I)
  100          CONTINUE
C
C              Compute the number of unstable zeros.
C
               K = K + 1
               IF ( DWORK(K).EQ.ZERO ) THEN
                  INFO = 2
               ELSE
                  SIGNUM = SIGNUM*SIGN( ONE, DWORK(K) )
                  IF ( SIGNUM.LT.ZERO ) NZ = NZ + 1
               END IF
               GO TO 80
            END IF
C           END WHILE 80
         END IF
      END IF
C
      IF ( ( INFO.EQ.0 ) .AND. ( NZ.EQ.0 ) ) THEN
         STABLE = .TRUE.
      ELSE
         STABLE = .FALSE.
      END IF
C
      RETURN
C *** Last line of MC01TD ***
      END
