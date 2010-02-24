      SUBROUTINE MB02UV( N, A, LDA, IPIV, JPIV, INFO )
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
C     To compute an LU factorization, using complete pivoting, of the
C     N-by-N matrix A. The factorization has the form A = P * L * U * Q,
C     where P and Q are permutation matrices, L is lower triangular with
C     unit diagonal elements and U is upper triangular.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A to be factored.
C             On exit, the leading N-by-N part of this array contains
C             the factors L and U from the factorization A = P*L*U*Q;
C             the unit diagonal elements of L are not stored. If U(k, k)
C             appears to be less than SMIN, U(k, k) is given the value
C             of SMIN, giving a nonsingular perturbed system.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1, N).
C
C     IPIV    (output) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= i <= N, row i of the
C             matrix has been interchanged with row IPIV(i).
C
C     JPIV    (output) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= j <= N, column j of the
C             matrix has been interchanged with column JPIV(j).
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = k:  U(k, k) is likely to produce owerflow if one tries
C                   to solve for x in Ax = b. So U is perturbed to get
C                   a nonsingular system. This is a warning.
C
C     FURTHER COMMENTS
C
C     In the interests of speed, this routine does not check the input
C     for errors. It should only be used to factorize matrices A of
C     very small order.
C
C     CONTRIBUTOR
C
C     Bo Kagstrom and Peter Poromaa, Univ. of Umea, Sweden, Nov. 1993.
C
C     REVISIONS
C
C     April 1998 (T. Penzl).
C     Sep. 1998 (V. Sima).
C     March 1999 (V. Sima).
C     March 2004 (V. Sima).
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
C     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * )
C     .. Local Scalars ..
      INTEGER            I, IP, IPV, JP, JPV
      DOUBLE PRECISION   BIGNUM, EPS, SMIN, SMLNUM, XMAX
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     .. External Subroutines ..
      EXTERNAL           DGER, DLABAD, DSCAL, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     .. Executable Statements ..
C
C     Set constants to control owerflow.

      INFO = 0
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = DLAMCH( 'Safe minimum' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
C     Find max element in matrix A.
C
      IPV = 1
      JPV = 1
      XMAX = ZERO
      DO 40 JP = 1, N
         DO 20 IP = 1, N
            IF ( ABS( A(IP, JP) ) .GT. XMAX ) THEN
               XMAX = ABS( A(IP, JP) )
               IPV  = IP
               JPV  = JP
            ENDIF
   20    CONTINUE
   40 CONTINUE
      SMIN =  MAX( EPS * XMAX, SMLNUM )
C
C     Swap rows.
C
      IF ( IPV .NE. 1 ) CALL DSWAP( N, A(IPV, 1), LDA, A(1, 1), LDA )
      IPIV(1) = IPV
C
C     Swap columns.
C
      IF ( JPV .NE. 1 ) CALL DSWAP( N, A(1, JPV), 1, A(1, 1), 1 )
      JPIV(1) = JPV
C
C     Check for singularity.
C
      IF ( ABS( A(1, 1) ) .LT. SMIN ) THEN
         INFO = 1
         A(1, 1) = SMIN
      ENDIF
      IF ( N.GT.1 ) THEN
         CALL DSCAL( N - 1, ONE / A(1, 1), A(2, 1), 1 )
         CALL DGER( N - 1, N - 1, -ONE, A(2, 1), 1, A(1, 2), LDA,
     $              A(2, 2), LDA )
      ENDIF
C
C     Factorize the rest of A with complete pivoting.
C     Set pivots less than SMIN to SMIN.
C
      DO 100 I = 2, N - 1
C
C        Find max element in remaining matrix.
C
         IPV = I
         JPV = I
         XMAX = ZERO
         DO 80 JP = I, N
            DO 60 IP = I, N
               IF ( ABS( A(IP, JP) ) .GT. XMAX ) THEN
                  XMAX = ABS( A(IP, JP) )
                  IPV = IP
                  JPV = JP
               ENDIF
   60       CONTINUE
   80    CONTINUE
C
C        Swap rows.
C
         IF ( IPV .NE. I ) CALL DSWAP( N, A(IPV, 1), LDA, A(I, 1), LDA )
         IPIV(I) = IPV
C
C        Swap columns.
C
         IF ( JPV .NE. I ) CALL DSWAP( N, A(1, JPV), 1, A(1, I), 1 )
         JPIV(I) = JPV
C
C        Check for almost singularity.
C
         IF ( ABS( A(I, I) ) .LT. SMIN ) THEN
            INFO = I
            A(I, I) = SMIN
         ENDIF
         CALL DSCAL( N - I, ONE / A(I, I), A(I + 1, I), 1 )
         CALL DGER( N - I, N - I, -ONE, A(I + 1, I), 1, A(I, I + 1),
     $              LDA, A(I + 1, I + 1), LDA )
  100 CONTINUE
      IF ( ABS( A(N, N) ) .LT. SMIN ) THEN
         INFO = N
         A(N, N) = SMIN
      ENDIF
C
      RETURN
C *** Last line of MB02UV ***
      END
