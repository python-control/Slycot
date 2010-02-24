      SUBROUTINE TD05AD( UNITF, OUTPUT, NP1, MP1, W, A, B, VALR, VALI,
     $                   INFO )
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
C     Given a complex valued rational function of frequency (transfer
C     function) G(jW) this routine will calculate its complex value or
C     its magnitude and phase for a specified frequency value.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UNITF   CHARACTER*1
C             Indicates the choice of frequency unit as follows:
C             = 'R':  Input frequency W in radians/second;
C             = 'H':  Input frequency W in hertz.
C
C     OUTPUT  CHARACTER*1
C             Indicates the choice of co-ordinates for output as folows:
C             = 'C':  Cartesian co-ordinates (output real and imaginary
C                     parts of G(jW));
C             = 'P':  Polar co-ordinates (output magnitude and phase
C                     of G(jW)).
C
C     Input/Output Parameters
C
C     NP1     (input) INTEGER
C             The order of the denominator + 1, i.e. N + 1.  NP1 >= 1.
C
C     MP1     (input) INTEGER
C             The order of the numerator + 1, i.e. M + 1.  MP1 >= 1.
C
C     W       (input) DOUBLE PRECISION
C             The frequency value W for which the transfer function is
C             to be evaluated.
C
C     A       (input) DOUBLE PRECISION array, dimension (NP1)
C             This array must contain the vector of denominator
C             coefficients in ascending order of powers. That is, A(i)
C             must contain the coefficient of (jW)**(i-1) for i = 1,
C             2,...,NP1.
C
C     B       (input) DOUBLE PRECISION array, dimension (MP1)
C             This array must contain the vector of numerator
C             coefficients in ascending order of powers. That is, B(i)
C             must contain the coefficient of (jW)**(i-1) for i = 1,
C             2,...,MP1.
C
C     VALR    (output) DOUBLE PRECISION
C             If OUTPUT = 'C', VALR contains the real part of G(jW).
C             If OUTPUT = 'P', VALR contains the magnitude of G(jW)
C                              in dBs.
C
C     VALI    (output) DOUBLE PRECISION
C             If OUTPUT = 'C', VALI contains the imaginary part of
C                              G(jW).
C             If OUTPUT = 'P', VALI contains the phase of G(jW) in
C                              degrees.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the frequency value W is a pole of G(jW), or all
C                   the coefficients of the A polynomial are zero.
C
C     METHOD
C
C     By substituting the values of A, B and W in the following
C     formula:
C
C            B(1)+B(2)*(jW)+B(3)*(jW)**2+...+B(MP1)*(jW)**(MP1-1)
C     G(jW) = ---------------------------------------------------.
C            A(1)+A(2)*(jW)+A(3)*(jW)**2+...+A(NP1)*(jW)**(NP1-1)
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0(N+M) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TD01AD by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, March 1981.
C
C     REVISIONS
C
C     February 1997.
C     February 22, 1998 (changed the name of TD01MD).
C
C     KEYWORDS
C
C     Elementary polynomial operations, frequency response, matrix
C     fraction, polynomial matrix, state-space representation, transfer
C     matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, EIGHT, TWENTY, NINETY, ONE80, THRE60
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, EIGHT=8.0D0,
     $                    TWENTY=20.0D0, NINETY=90.0D0, ONE80 = 180.0D0,
     $                    THRE60=360.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         OUTPUT, UNITF
      INTEGER           INFO, MP1, NP1
      DOUBLE PRECISION  VALI, VALR, W
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*)
C     .. Local Scalars ..
      LOGICAL           LOUTPU, LUNITF
      INTEGER           I, IPHASE, M, M2, N, N2, NPZERO, NZZERO
      DOUBLE PRECISION  BIMAG, BREAL, G, TIMAG, TREAL, TWOPI, W2, WC
      COMPLEX*16        ZTEMP
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAPY2
      COMPLEX*16        ZLADIV
      EXTERNAL          DLAPY2, LSAME, ZLADIV
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, ATAN, DBLE, DCMPLX, DIMAG, LOG10, MAX, MOD,
     $                  SIGN
C     .. Executable Statements ..
C
      INFO = 0
      LUNITF = LSAME( UNITF,  'H' )
      LOUTPU = LSAME( OUTPUT, 'P' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LUNITF .AND. .NOT.LSAME( UNITF,  'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOUTPU .AND. .NOT.LSAME( OUTPUT, 'C' ) ) THEN
         INFO = -2
      ELSE IF( NP1.LT.1 ) THEN
         INFO = -3
      ELSE IF( MP1.LT.1 ) THEN
         INFO = -4
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TD05AD', -INFO )
         RETURN
      END IF
C
      M = MP1 - 1
      N = NP1 - 1
      WC = W
      TWOPI = EIGHT*ATAN( ONE )
      IF ( LUNITF ) WC = WC*TWOPI
      W2 = WC**2
C
C     Determine the orders z (NZZERO) and p (NPZERO) of the factors
C     (jW)**k in the numerator and denominator polynomials, by counting
C     the zero trailing coefficients.  The value of G(jW) will then be
C     computed as (jW)**(z-p)*m(jW)/n(jW), for appropriate m and n.
C
      I = 0
C
   10 CONTINUE
      I = I + 1
      IF ( I.LE.M ) THEN
         IF ( B(I).EQ.ZERO ) GO TO 10
      END IF
C
      NZZERO = I - 1
      I = 0
C
   20 CONTINUE
      I = I + 1
      IF ( I.LE.N ) THEN
         IF ( A(I).EQ.ZERO ) GO TO 20
      END IF
C
      NPZERO = I - 1
      IPHASE = NZZERO - NPZERO
C
      M2 = MOD( M - NZZERO, 2 )
C
C     Add real parts of the numerator m(jW).
C
      TREAL = B(MP1-M2)
C
      DO 30 I = M - 1 - M2, NZZERO + 1, -2
         TREAL = B(I) - W2*TREAL
   30 CONTINUE
C
C     Add imaginary parts of the numerator m(jW).
C
      IF ( M.EQ.0 ) THEN
         TIMAG = ZERO
      ELSE
         TIMAG = B(M+M2)
C
         DO 40 I = M + M2 - 2, NZZERO + 2, -2
            TIMAG = B(I) - W2*TIMAG
   40    CONTINUE
C
         TIMAG = TIMAG*WC
      END IF
C
      N2 = MOD( N - NPZERO, 2 )
C
C     Add real parts of the denominator n(jW).
C
      BREAL = A(NP1-N2)
C
      DO 50 I = N - 1 - N2, NPZERO + 1, -2
         BREAL = A(I) - W2*BREAL
   50 CONTINUE
C
C     Add imaginary parts of the denominator n(jW).
C
      IF ( N.EQ.0 ) THEN
         BIMAG = ZERO
      ELSE
         BIMAG = A(N+N2)
C
         DO 60 I = N + N2 - 2, NPZERO + 2, -2
            BIMAG = A(I) - W2*BIMAG
   60    CONTINUE
C
         BIMAG = BIMAG*WC
      END IF
C
      IF ( ( MAX( ABS( BREAL ), ABS( BIMAG ) ).EQ.ZERO ) .OR.
     $     ( W.EQ.ZERO .AND. IPHASE.LT.0 ) ) THEN
C
C        Error return:  The specified frequency W is a pole of G(jW),
C              or all the coefficients of the A polynomial are zero.
C
         INFO = 1
      ELSE
C
C        Evaluate the complex number W**(z-p)*m(jW)/n(jW).
C
         ZTEMP =
     $      ZLADIV( DCMPLX( TREAL, TIMAG ), DCMPLX( BREAL, BIMAG ) )
         VALR = DBLE(  ZTEMP )*WC**IPHASE
         VALI = DIMAG( ZTEMP )*WC**IPHASE
C
         IF ( .NOT.LOUTPU ) THEN
C
C           Cartesian co-ordinates: Update the result for j**(z-p).
C
            I = MOD( ABS( IPHASE ), 4 )
            IF ( ( IPHASE.GT.0 .AND. I.GT.1 ) .OR.
     $           ( IPHASE.LT.0 .AND. ( I.EQ.1 .OR. I.EQ.2) ) ) THEN
               VALR = -VALR
               VALI = -VALI
            END IF
C
            IF ( MOD( I, 2 ).NE.0 ) THEN
               G    =  VALR
               VALR = -VALI
               VALI =  G
            END IF
C
         ELSE
C
C           Polar co-ordinates: Compute the magnitude and phase.
C
            G = DLAPY2( VALR, VALI )
C
            IF ( VALR.EQ.ZERO ) THEN
               VALI = SIGN( NINETY, VALI )
            ELSE
               VALI = ( ATAN( VALI/VALR )/TWOPI )*THRE60
               IF ( VALI.EQ.ZERO .AND. NZZERO.EQ.M .AND. NPZERO.EQ.N
     $                           .AND. B(NZZERO+1)*A(NPZERO+1).LT.ZERO )
     $            VALI = ONE80
            END IF
C
            VALR = TWENTY*LOG10( G )
C
            IF ( IPHASE.NE.0 )
     $         VALI = VALI + DBLE( NZZERO - NPZERO )*NINETY
         END IF
C
      END IF
C
      RETURN
C *** Last line of TD05AD ***
      END
