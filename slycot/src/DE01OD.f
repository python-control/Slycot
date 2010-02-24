      SUBROUTINE DE01OD( CONV, N, A, B, INFO )
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
C     To compute the convolution or deconvolution of two real signals
C     A and B.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     CONV    CHARACTER*1
C             Indicates whether convolution or deconvolution is to be
C             performed as follows:
C             = 'C':  Convolution;
C             = 'D':  Deconvolution.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of samples.  N must be a power of 2.  N >= 2.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the first signal.
C             On exit, this array contains the convolution (if
C             CONV = 'C') or deconvolution (if CONV = 'D') of the two
C             signals.
C
C     B       (input) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the second signal.
C             NOTE that this array is overwritten.
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
C     This routine computes the convolution or deconvolution of two real
C     signals A and B using an FFT algorithm (SLICOT Library routine
C     DG01MD).
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
C     Supersedes Release 2.0 routine DE01CD by R. Dekeyser, State
C     University of Gent, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Convolution, deconvolution, digital signal processing, fast
C     Fourier transform, real signals.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         ( ZERO = 0.0D0, HALF=0.5D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         CONV
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*)
C     .. Local Scalars ..
      LOGICAL           LCONV
      INTEGER           J, KJ, ND2P1
      DOUBLE PRECISION  AC, AS, AST, BC, BS, CI, CR
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DG01MD, DLADIV, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MOD
C     .. Executable Statements ..
C
      INFO = 0
      LCONV = LSAME( CONV, 'C' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LCONV .AND. .NOT.LSAME( CONV, 'D' ) ) THEN
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
         CALL XERBLA( 'DE01OD', -INFO )
         RETURN
      END IF
C
C     Fourier transform.
C
      CALL DG01MD( 'Direct', N, A, B, INFO )
C
      IF ( LCONV ) THEN
         AST = A(1)*B(1)
      ELSE
         IF ( B(1).EQ.ZERO ) THEN
            AST = ZERO
         ELSE
            AST = A(1)/B(1)
         END IF
      END IF
C
      ND2P1 = N/2 + 1
      J = ND2P1
C
      DO 20 KJ = ND2P1, N
C
C        Components of the transform of function A.
C
         AC = HALF*( A(J) + A(KJ) )
         AS = HALF*( B(J) - B(KJ) )
C
C        Components of the transform of function B.
C
         BC = HALF*( B(KJ) + B(J) )
         BS = HALF*( A(KJ) - A(J) )
C
C        Deconvolution by complex division if CONV = 'D';
C        Convolution by complex multiplication if CONV = 'C'.
C
         IF ( LCONV ) THEN
            CR = AC*BC - AS*BS
            CI = AS*BC + AC*BS
         ELSE
            IF ( MAX( ABS( BC ), ABS( BS ) ).EQ.ZERO ) THEN
               CR = ZERO
               CI = ZERO
            ELSE
               CALL DLADIV( AC, AS, BC, BS, CR, CI )
            END IF
         END IF
C
         A(J)  =  CR
         B(J)  =  CI
         A(KJ) =  CR
         B(KJ) = -CI
         J = J - 1
   20 CONTINUE
      A(1) = AST
      B(1) = ZERO
C
C     Inverse Fourier transform.
C
      CALL DG01MD( 'Inverse', N, A, B, INFO )
C
      CALL DSCAL( N, ONE/DBLE( N ), A, 1 )
C
      RETURN
C *** Last line of DE01OD ***
      END
