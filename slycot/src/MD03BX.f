      SUBROUTINE MD03BX( M, N, FNORM, J, LDJ, E, JNORMS, GNORM, IPVT,
     $                   DWORK, LDWORK, INFO )
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
C     To compute the QR factorization with column pivoting of an
C     m-by-n matrix J (m >= n), that is, J*P = Q*R, where Q is a matrix
C     with orthogonal columns, P a permutation matrix, and R an upper
C     trapezoidal matrix with diagonal elements of nonincreasing
C     magnitude, and to apply the transformation Q' on the error
C     vector e (in-situ). The 1-norm of the scaled gradient is also
C     returned. The matrix J could be the Jacobian of a nonlinear least
C     squares problem.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the Jacobian matrix J.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the Jacobian matrix J.
C             M >= N >= 0.
C
C     FNORM   (input) DOUBLE PRECISION
C             The Euclidean norm of the vector e.  FNORM >= 0.
C
C     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, N)
C             On entry, the leading M-by-N part of this array must
C             contain the Jacobian matrix J.
C             On exit, the leading N-by-N upper triangular part of this
C             array contains the upper triangular factor R of the
C             Jacobian matrix. Note that for efficiency of the later
C             calculations, the matrix R is delivered with the leading
C             dimension MAX(1,N), possibly much smaller than the value
C             of LDJ on entry.
C
C     LDJ     (input/output) INTEGER
C             The leading dimension of array J.
C             On entry, LDJ >= MAX(1,M).
C             On exit,  LDJ >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (M)
C             On entry, this array must contain the error vector e.
C             On exit, this array contains the updated vector Q'*e.
C
C     JNORMS  (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the Euclidean norms of the columns of
C             the Jacobian matrix, considered in the initial order.
C
C     GNORM   (output) DOUBLE PRECISION
C             If FNORM > 0, the 1-norm of the scaled vector
C             J'*Q'*e/FNORM, with each element i further divided by
C             JNORMS(i) (if JNORMS(i) is nonzero).
C             If FNORM = 0, the returned value of GNORM is 0.
C
C     IPVT    (output) INTEGER array, dimension (N)
C             This array defines the permutation matrix P such that
C             J*P = Q*R. Column j of P is column IPVT(j) of the identity
C             matrix.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 1,      if N = 0 or M = 1;
C             LDWORK >= 4*N+1,  if N > 1.
C             For optimum performance LDWORK should be larger.
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
C     The algorithm uses QR factorization with column pivoting of the
C     matrix J, J*P = Q*R, and applies the orthogonal matrix Q' to the
C     vector e.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, Jacobian matrix, matrix algebra,
C     matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDJ, LDWORK, M, N
      DOUBLE PRECISION  FNORM, GNORM
C     .. Array Arguments ..
      INTEGER           IPVT(*)
      DOUBLE PRECISION  DWORK(*), E(*), J(*), JNORMS(*)
C     .. Local Scalars ..
      INTEGER           I, ITAU, JWORK, L, WRKOPT
      DOUBLE PRECISION  SUM
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DGEQP3, DLACPY, DORMQR, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      IF ( M.LT.0 ) THEN
         INFO = -1
      ELSEIF ( N.LT.0.OR. M.LT.N ) THEN
         INFO = -2
      ELSEIF ( FNORM.LT.ZERO ) THEN
         INFO = -3
      ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE
         IF ( N.EQ.0 .OR. M.EQ.1 ) THEN
            JWORK = 1
         ELSE
            JWORK = 4*N + 1
         END IF
         IF ( LDWORK.LT.JWORK )
     $      INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MD03BX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      GNORM = ZERO
      IF ( N.EQ.0 ) THEN
         LDJ = 1
         DWORK(1) = ONE
         RETURN
      ELSEIF ( M.EQ.1 ) THEN
         JNORMS(1) = ABS( J(1) )
         IF ( FNORM*J(1).NE.ZERO )
     $      GNORM = ABS( E(1)/FNORM )
         LDJ      = 1
         IPVT(1)  = 1
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Initialize the column pivoting indices.
C
      DO 10 I = 1, N
         IPVT(I) = 0
   10 CONTINUE
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      ITAU   = 1
      JWORK  = ITAU + N
      WRKOPT = 1
C
C     Compute the QR factorization with pivoting of J, and apply Q' to
C     the vector e.
C
C     Workspace: need:    4*N + 1;
C                prefer:  3*N + ( N+1 )*NB.
C
      CALL DGEQP3( M, N, J, LDJ, IPVT, DWORK(ITAU), DWORK(JWORK),
     $             LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C     Workspace: need:    N + 1;
C                prefer:  N + NB.
C
      CALL DORMQR( 'Left', 'Transpose', M, 1, N, J, LDJ, DWORK(ITAU), E,
     $             M, DWORK(JWORK), LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
      IF ( LDJ.GT.N ) THEN
C
C        Reshape the array J to have the leading dimension N.
C        This destroys the details of the orthogonal matrix Q.
C
         CALL DLACPY( 'Upper', N, N, J, LDJ, J, N )
         LDJ = N
      END IF
C
C     Compute the norm of the scaled gradient and original column norms.
C
      IF ( FNORM.NE.ZERO ) THEN
C
         DO 20 I = 1, N
            L = IPVT(I)
            JNORMS(L) = DNRM2( I, J((I-1)*LDJ+1), 1 )
            IF ( JNORMS(L).NE.ZERO ) THEN
               SUM = DDOT( I, J((I-1)*LDJ+1), 1, E, 1 )/FNORM
               GNORM = MAX( GNORM, ABS( SUM/JNORMS(L) ) )
            END IF
   20    CONTINUE
C
      ELSE
C
         DO 30 I = 1, N
            L = IPVT(I)
            JNORMS(L) = DNRM2( I, J((I-1)*LDJ+1), 1 )
   30    CONTINUE
C
      END IF
C
      DWORK(1) = WRKOPT
      RETURN
C
C *** Last line of MD03BX ***
      END
