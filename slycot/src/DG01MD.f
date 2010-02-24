      SUBROUTINE DG01MD( INDI, N, XR, XI, INFO )
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
C     To compute the discrete Fourier transform, or inverse transform,
C     of a complex signal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     INDI    CHARACTER*1
C             Indicates whether a Fourier transform or inverse Fourier
C             transform is to be performed as follows:
C             = 'D':  (Direct) Fourier transform;
C             = 'I':  Inverse Fourier transform.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of complex samples.  N must be a power of 2.
C             N >= 2.
C
C     XR      (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the real part of either
C             the complex signal z if INDI = 'D', or f(z) if INDI = 'I'.
C             On exit, this array contains either the real part of the
C             computed Fourier transform f(z) if INDI = 'D', or the
C             inverse Fourier transform z of f(z) if INDI = 'I'.
C
C     XI      (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the imaginary part of
C             either z if INDI = 'D', or f(z) if INDI = 'I'.
C             On exit, this array contains either the imaginary part of
C             f(z) if INDI = 'D', or z if INDI = 'I'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     If INDI = 'D', then the routine performs a discrete Fourier
C     transform on the complex signal Z(i), i = 1,2,...,N. If the result
C     is denoted by FZ(k), k = 1,2,...,N, then the relationship between
C     Z and FZ is given by the formula:
C
C                     N            ((k-1)*(i-1))
C            FZ(k) = SUM ( Z(i) * V              ),
C                    i=1
C                                     2
C     where V = exp( -2*pi*j/N ) and j  = -1.
C
C     If INDI = 'I', then the routine performs an inverse discrete
C     Fourier transform on the complex signal FZ(k), k = 1,2,...,N. If
C     the result is denoted by Z(i), i = 1,2,...,N, then the
C     relationship between Z and FZ is given by the formula:
C
C                    N             ((k-1)*(i-1))
C            Z(i) = SUM ( FZ(k) * W              ),
C                   k=1
C
C     where W = exp( 2*pi*j/N ).
C
C     Note that a discrete Fourier transform, followed by an inverse
C     discrete Fourier transform, will result in a signal which is a
C     factor N larger than the original input signal.
C
C     REFERENCES
C
C     [1] Rabiner, L.R. and Rader, C.M.
C         Digital Signal Processing.
C         IEEE Press, 1972.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0( N*log(N) ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine DG01AD by R. Dekeyser, State
C     University of Gent, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Complex signals, digital signal processing, fast Fourier
C     transform.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, EIGHT
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                    TWO = 2.0D0, EIGHT = 8.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         INDI
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  XI(*), XR(*)
C     .. Local Scalars ..
      LOGICAL           LINDI
      INTEGER           I, J, K, L, M
      DOUBLE PRECISION  PI2, TI, TR, WHELP, WI, WR, WSTPI, WSTPR
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, DBLE, MOD, SIN
C     .. Executable Statements ..
C
      INFO = 0
      LINDI = LSAME( INDI, 'D' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LINDI .AND. .NOT.LSAME( INDI, 'I' ) ) THEN
         INFO = -1
      ELSE
         J = 0
         IF( N.GE.2 ) THEN
            J = N
C           WHILE ( MOD( J, 2 ).EQ.0 ) DO
   10       CONTINUE
            IF ( MOD( J, 2 ).EQ.0 ) THEN
               J = J/2
               GO TO 10
            END IF
C           END WHILE 10
         END IF
         IF ( J.NE.1 ) INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DG01MD', -INFO )
         RETURN
      END IF
C
C     Inplace shuffling of data.
C
      J = 1
C
      DO 30 I = 1, N
         IF ( J.GT.I ) THEN
            TR = XR(I)
            TI = XI(I)
            XR(I) = XR(J)
            XI(I) = XI(J)
            XR(J) = TR
            XI(J) = TI
         END IF
         K = N/2
C        REPEAT
   20    IF ( J.GT.K ) THEN
            J = J - K
            K = K/2
            IF ( K.GE.2 ) GO TO 20
         END IF
C        UNTIL ( K.LT.2 )
         J = J + K
   30 CONTINUE
C
C     Transform by decimation in time.
C
      PI2 = EIGHT*ATAN( ONE )
      IF ( LINDI ) PI2 = -PI2
C
      I = 1
C
C     WHILE ( I.LT.N ) DO
C
   40 IF ( I.LT.N ) THEN
         L = 2*I
         WHELP = PI2/DBLE( L )
         WSTPI = SIN( WHELP )
         WHELP = SIN( HALF*WHELP )
         WSTPR = -TWO*WHELP*WHELP
         WR = ONE
         WI = ZERO
C
         DO 60 J = 1, I
C
            DO 50 K = J, N, L
               M = K + I
               TR = WR*XR(M) - WI*XI(M)
               TI = WR*XI(M) + WI*XR(M)
               XR(M) = XR(K) - TR
               XI(M) = XI(K) - TI
               XR(K) = XR(K) + TR
               XI(K) = XI(K) + TI
   50       CONTINUE
C
            WHELP = WR
            WR = WR + WR*WSTPR - WI*WSTPI
            WI = WI + WHELP*WSTPI + WI*WSTPR
   60    CONTINUE
C
         I = L
         GO TO 40
C        END WHILE 40
      END IF
C
      RETURN
C *** Last line of DG01MD ***
      END
