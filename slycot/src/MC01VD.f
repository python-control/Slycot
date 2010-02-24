      SUBROUTINE MC01VD( A, B, C, Z1RE, Z1IM, Z2RE, Z2IM, INFO )
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
C     To compute the roots of a quadratic equation with real
C     coefficients.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     A       (input) DOUBLE PRECISION
C             The value of the coefficient of the quadratic term.
C
C     B       (input) DOUBLE PRECISION
C             The value of the coefficient of the linear term.
C
C     C       (input) DOUBLE PRECISION
C             The value of the coefficient of the constant term.
C
C     Z1RE    (output) DOUBLE PRECISION
C     Z1IM    (output) DOUBLE PRECISION
C             The real and imaginary parts, respectively, of the largest
C             root in magnitude.
C
C     Z2RE    (output) DOUBLE PRECISION
C     Z2IM    (output) DOUBLE PRECISION
C             The real and imaginary parts, respectively, of the
C             smallest root in magnitude.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if on entry, either A = B = 0.0 or A = 0.0 and the
C                   root -C/B overflows; in this case Z1RE, Z1IM, Z2RE
C                   and Z2IM are unassigned;
C             = 2:  if on entry, A = 0.0; in this case Z1RE contains
C                   BIG and Z1IM contains zero, where BIG is a
C                   representable number near the overflow threshold
C                   of the machine (see LAPACK Library Routine DLAMCH);
C             = 3:  if on entry, either C = 0.0 and the root -B/A
C                   overflows or A, B and C are non-zero and the largest
C                   real root in magnitude cannot be computed without
C                   overflow; in this case Z1RE contains BIG and Z1IM
C                   contains zero;
C             = 4:  if the roots cannot be computed without overflow; in
C                   this case Z1RE, Z1IM, Z2RE and Z2IM are unassigned.
C
C     METHOD
C
C     The routine computes the roots (r1 and r2) of the real quadratic
C     equation
C             2
C        a * x  + b * x + c = 0
C
C     as
C             - b - SIGN(b) * SQRT(b * b - 4 * a * c)             c
C        r1 = ---------------------------------------  and r2 = ------
C                              2 * a                            a * r1
C
C     unless a = 0, in which case
C
C             -c
C        r1 = --.
C              b
C
C     Precautions are taken to avoid overflow and underflow wherever
C     possible.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01JD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Quadratic equation, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, FOUR
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, FOUR=4.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO
      DOUBLE PRECISION  A, B, C, Z1IM, Z1RE, Z2IM, Z2RE
C     .. Local Scalars ..
      LOGICAL           OVFLOW
      INTEGER           BETA, EA, EAPLEC, EB, EB2, EC, ED
      DOUBLE PRECISION  ABSA, ABSB, ABSC, BIG, M1, M2, MA, MB, MC, MD,
     $                  SFMIN, W
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          MC01SW, MC01SY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MOD, SIGN, SQRT
C     .. Executable Statements ..
C
C     Detect special cases.
C
      INFO  = 0
      BETA  = DLAMCH( 'Base' )
      SFMIN = DLAMCH( 'Safe minimum' )
      BIG   = ONE/SFMIN
      IF ( A.EQ.ZERO ) THEN
         IF ( B.EQ.ZERO ) THEN
            INFO = 1
         ELSE
            OVFLOW = .FALSE.
            Z2RE   =  ZERO
            IF ( C.NE.ZERO ) THEN
               ABSB = ABS( B )
               IF ( ABSB.GE.ONE ) THEN
                  IF ( ABS( C ).GE.ABSB*SFMIN ) Z2RE = -C/B
               ELSE
                  IF ( ABS( C ).LE.ABSB*BIG ) THEN
                     Z2RE = -C/B
                  ELSE
                     OVFLOW = .TRUE.
                     Z2RE   = BIG
                     IF ( SIGN( ONE, B )*SIGN( ONE, C ).GT.ZERO )
     $                  Z2RE = -BIG
                  END IF
               END IF
            END IF
            IF ( OVFLOW ) THEN
               INFO = 1
            ELSE
               Z1RE = BIG
               Z1IM = ZERO
               Z2IM = ZERO
               INFO = 2
            END IF
         END IF
         RETURN
      END IF
C
      IF ( C.EQ.ZERO ) THEN
         OVFLOW = .FALSE.
         Z1RE   =  ZERO
         IF ( B.NE.ZERO ) THEN
            ABSA = ABS( A )
            IF ( ABSA.GE.ONE ) THEN
               IF ( ABS( B ).GE.ABSA*SFMIN ) Z1RE = -B/A
            ELSE
               IF ( ABS( B ).LE.ABSA*BIG ) THEN
                  Z1RE = -B/A
               ELSE
                  OVFLOW = .TRUE.
                  Z1RE   = BIG
               END IF
            END IF
         END IF
         IF ( OVFLOW ) INFO = 3
         Z1IM = ZERO
         Z2RE = ZERO
         Z2IM = ZERO
         RETURN
      END IF
C
C     A and C are non-zero.
C
      IF ( B.EQ.ZERO ) THEN
         OVFLOW = .FALSE.
         ABSC   = SQRT( ABS( C ) )
         ABSA   = SQRT( ABS( A ) )
         W =  ZERO
         IF ( ABSA.GE.ONE ) THEN
            IF ( ABSC.GE.ABSA*SFMIN ) W = ABSC/ABSA
         ELSE
            IF ( ABSC.LE.ABSA*BIG ) THEN
               W = ABSC/ABSA
            ELSE
               OVFLOW = .TRUE.
               W = BIG
            END IF
         END IF
         IF ( OVFLOW ) THEN
            INFO = 4
         ELSE
            IF ( SIGN( ONE, A )*SIGN( ONE, C ).GT.ZERO ) THEN
               Z1RE = ZERO
               Z2RE = ZERO
               Z1IM =  W
               Z2IM = -W
            ELSE
               Z1RE =  W
               Z2RE = -W
               Z1IM = ZERO
               Z2IM = ZERO
            END IF
         END IF
         RETURN
      END IF
C
C     A, B and C are non-zero.
C
      CALL MC01SW( A, BETA, MA, EA )
      CALL MC01SW( B, BETA, MB, EB )
      CALL MC01SW( C, BETA, MC, EC )
C
C     Compute a 'near' floating-point representation of the discriminant
C     D = MD * BETA**ED.
C
      EAPLEC = EA + EC
      EB2 = 2*EB
      IF ( EAPLEC.GT.EB2 ) THEN
         CALL MC01SY( MB*MB, EB2-EAPLEC, BETA, W, OVFLOW )
         W = W - FOUR*MA*MC
         CALL MC01SW( W, BETA, MD, ED )
         ED = ED + EAPLEC
      ELSE
         CALL MC01SY( FOUR*MA*MC, EAPLEC-EB2, BETA, W, OVFLOW )
         W = MB*MB - W
         CALL MC01SW( W, BETA, MD, ED )
         ED = ED + EB2
      END IF
C
      IF ( MOD( ED, 2 ).NE.0 ) THEN
         ED = ED + 1
         MD = MD/BETA
      END IF
C
C     Complex roots.
C
      IF ( MD.LT.ZERO ) THEN
         CALL MC01SY( -MB/( 2*MA ), EB-EA, BETA, Z1RE, OVFLOW )
         IF ( OVFLOW ) THEN
            INFO = 4
         ELSE
            CALL MC01SY( SQRT( -MD )/( 2*MA ), ED/2-EA, BETA, Z1IM,
     $                   OVFLOW )
            IF ( OVFLOW ) THEN
               INFO = 4
            ELSE
               Z2RE =  Z1RE
               Z2IM = -Z1IM
            END IF
         END IF
         RETURN
      END IF
C
C     Real roots.
C
      MD = SQRT( MD )
      ED = ED/2
      IF ( ED.GT.EB ) THEN
         CALL MC01SY( ABS( MB ), EB-ED, BETA, W, OVFLOW )
         W = W + MD
         M1 = -SIGN( ONE, MB )*W/( 2*MA )
         CALL MC01SY( M1, ED-EA, BETA, Z1RE, OVFLOW )
         IF ( OVFLOW ) THEN
            Z1RE = BIG
            INFO = 3
         END IF
         M2 = -SIGN( ONE, MB )*2*MC/W
         CALL MC01SY( M2, EC-ED, BETA, Z2RE, OVFLOW )
      ELSE
         CALL MC01SY( MD, ED-EB, BETA, W, OVFLOW )
         W  = W + ABS( MB )
         M1 = -SIGN( ONE, MB )*W/( 2*MA )
         CALL MC01SY( M1, EB-EA, BETA, Z1RE, OVFLOW )
         IF ( OVFLOW ) THEN
            Z1RE = BIG
            INFO = 3
         END IF
         M2 = -SIGN( ONE, MB )*2*MC/W
         CALL MC01SY( M2, EC-EB, BETA, Z2RE, OVFLOW )
      END IF
      Z1IM = ZERO
      Z2IM = ZERO
C
      RETURN
C *** Last line of MC01VD ***
      END
