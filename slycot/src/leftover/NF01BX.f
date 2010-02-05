      SUBROUTINE NF01BX( N, IPAR, LIPAR, DPAR, LDPAR, J, LDJ, X, INCX,
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
C     To compute (J'*J + c*I)*x, where J is an m-by-n real matrix, c is
C     a real scalar, I is the n-by-n identity matrix, and x is a real
C     n-vector.
C
C     NOTE: this routine must have the same arguments as SLICOT Library
C     routine NF01BW.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of columns of the Jacobian matrix J.  N >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters describing the structure of the
C             matrix J, as follows:
C             IPAR(1) must contain the number of rows M of the Jacobian
C                     matrix J.  M >= 0.
C             IPAR is provided for compatibility with SLICOT Library
C             routine MD03AD.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 1.
C
C     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C             The real parameters needed for solving the problem.
C             The entry DPAR(1) must contain the real scalar c.
C
C     LDPAR   (input) INTEGER
C             The length of the array DPAR.  LDPAR >= 1.
C
C     J       (input) DOUBLE PRECISION array, dimension (LDJ,N)
C             The leading M-by-N part of this array must contain the
C             Jacobian matrix J.
C
C     LDJ     INTEGER
C             The leading dimension of the array J.  LDJ >= MAX(1,M).
C
C     X       (input/output) DOUBLE PRECISION array, dimension
C             (1+(N-1)*abs(INCX))
C             On entry, this incremented array must contain the
C             vector x.
C             On exit, this incremented array contains the value of the
C             matrix-vector product (J'*J + c*I)*x.
C
C     INCX    (input) INTEGER
C             The increment for the elements of X.  INCX <> 0.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= M.
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
C     The associativity of matrix multiplications is used; the result
C     is obtained as:  x_out = J'*( J*x ) + c*x.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Mar. 2002, Oct. 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INCX, INFO, LDJ, LDPAR, LDWORK, LIPAR, N
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DPAR(*), DWORK(*), J(LDJ,*), X(*)
C     .. Local Scalars ..
      INTEGER           M
      DOUBLE PRECISION  C
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSEIF ( LIPAR.LT.1 ) THEN
         INFO = -3
      ELSEIF ( LDPAR.LT.1 ) THEN
         INFO = -5
      ELSEIF ( INCX.EQ.0 ) THEN
         INFO = -9
      ELSE
         M = IPAR(1)
         IF ( M.LT.0 ) THEN
            INFO = -2
         ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
            INFO = -7
         ELSEIF ( LDWORK.LT.M ) THEN
            INFO = -11
         ENDIF
      ENDIF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'NF01BX', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      C = DPAR(1)
      IF ( M.EQ.0 ) THEN
C
C        Special case, void J: x <-- c*x.
C
         CALL DSCAL( N, C, X, INCX )
         RETURN
      END IF
C
      CALL DGEMV( 'NoTranspose', M, N, ONE, J, LDJ, X, INCX, ZERO,
     $            DWORK, 1 )
      CALL DGEMV( 'Transpose', M, N, ONE, J, LDJ, DWORK, 1, C, X, INCX )
      RETURN
C
C *** Last line of NF01BX ***
      END
