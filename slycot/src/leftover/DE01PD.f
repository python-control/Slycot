      SUBROUTINE DE01PD( CONV, WGHT, N, A, B, W, INFO )
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
C     A and B using the Hartley transform.
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
C     WGHT    CHARACTER*1
C             Indicates whether the precomputed weights are available
C             or not, as follows:
C             = 'A':  available;
C             = 'N':  not available.
C             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is
C             set to 'A' on exit.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of samples.  N must be a power of 2.  N >= 0.
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
C     W       (input/output) DOUBLE PRECISION array,
C                            dimension (N - LOG2(N))
C             On entry with WGHT = 'A', this array must contain the long
C             weight vector computed by a previous call of this routine
C             or of the SLICOT Library routine DG01OD.f, with the same
C             value of N. If WGHT = 'N', the contents of this array on
C             entry is ignored.
C             On exit, this array contains the long weight vector.
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
C     This routine computes the convolution or deconvolution of two
C     real signals A and B using three scrambled Hartley transforms
C     (SLICOT Library routine DG01OD).
C
C     REFERENCES
C
C     [1] Van Loan, Charles.
C         Computational frameworks for the fast Fourier transform.
C         SIAM, 1992.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(N log(N)) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, April 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
C
C     KEYWORDS
C
C     Convolution, deconvolution, digital signal processing,
C     fast Hartley transform, real signals.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ONE, TWO
      PARAMETER         ( HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         CONV, WGHT
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*), W(*)
C     .. Local Scalars ..
      LOGICAL           LCONV, LWGHT
      INTEGER           J, L, LEN, M, P1, R1
      DOUBLE PRECISION  T1, T2, T3
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DG01OD, DLADIV, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MOD
C     .. Executable Statements ..
C
      INFO  = 0
      LCONV = LSAME( CONV, 'C' )
      LWGHT = LSAME( WGHT, 'A' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LCONV .AND. .NOT.LSAME( CONV, 'D' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LWGHT .AND. .NOT.LSAME( WGHT, 'N' ) ) THEN
         INFO = -2
      ELSE
         M = 0
         J = 0
         IF( N.GE.1 ) THEN
            J = N
C           WHILE ( MOD( J, 2 ).EQ.0 ) DO
   10       CONTINUE
            IF ( MOD( J, 2 ).EQ.0 ) THEN
               J = J/2
               M = M + 1
               GO TO 10
            END IF
C           END WHILE 10
            IF ( J.NE.1 ) INFO = -3
         ELSE IF ( N.LT.0 ) THEN
            INFO = -3
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DE01PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.LE.0 ) THEN
         RETURN
      ELSE IF ( N.EQ.1 ) THEN
         IF ( LCONV ) THEN
            A(1) = A(1)*B(1)
         ELSE
            A(1) = A(1)/B(1)
         END IF
         RETURN
      END IF
C
C     Scrambled Hartley transforms of A and B.
C
      CALL DG01OD( 'OutputScrambled', WGHT, N, A, W, INFO )
      CALL DG01OD( 'OutputScrambled', WGHT, N, B, W, INFO )
C
C     Something similar to a Hadamard product/quotient.
C
      LEN = 1
      IF( LCONV )  THEN
         A(1) = TWO*A(1)*B(1)
         A(2) = TWO*A(2)*B(2)
C
         DO 30 L = 1, M - 1
            LEN = 2*LEN
            R1  = 2*LEN
C
            DO 20 P1 = LEN + 1, LEN + LEN/2
               T1 = B(P1) + B(R1)
               T2 = B(P1) - B(R1)
               T3 = T2*A(P1)
               A(P1) = T1*A(P1) + T2*A(R1)
               A(R1) = T1*A(R1) - T3
               R1 = R1 - 1
   20       CONTINUE
C
   30    CONTINUE
C
      ELSE
C
         A(1) = HALF*A(1)/B(1)
         A(2) = HALF*A(2)/B(2)
C
         DO 50 L = 1, M - 1
            LEN = 2*LEN
            R1  = 2*LEN
C
            DO 40 P1 = LEN + 1, LEN + LEN/2
               CALL DLADIV( A(P1), A(R1), B(P1)+B(R1), B(R1)-B(P1), T1,
     $                      T2 )
               A(P1) = T1
               A(R1) = T2
               R1 = R1 - 1
   40       CONTINUE
C
   50    CONTINUE
C
      END IF
C
C     Transposed Hartley transform of A.
C
      CALL DG01OD( 'InputScrambled', WGHT, N, A, W, INFO )
      IF ( LCONV ) THEN
         CALL DSCAL( N, HALF/DBLE( N ), A, 1 )
      ELSE
         CALL DSCAL( N, TWO/DBLE( N ), A, 1 )
      END IF
C
      RETURN
C *** Last line of DE01PD ***
      END
