      SUBROUTINE MD03BB( COND, N, IPAR, LIPAR, R, LDR, IPVT, DIAG, QTB,
     $                   DELTA, PAR, RANKS, X, RX, TOL, DWORK, LDWORK,
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
C     To determine a value for the parameter PAR such that if x solves
C     the system
C
C           A*x = b ,     sqrt(PAR)*D*x = 0 ,
C
C     in the least squares sense, where A is an m-by-n matrix, D is an
C     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if
C     DELTA is a positive number, DXNORM is the Euclidean norm of D*x,
C     then either PAR is zero and
C
C           ( DXNORM - DELTA ) .LE. 0.1*DELTA ,
C
C     or PAR is positive and
C
C           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA .
C
C     It is assumed that a QR factorization, with column pivoting, of A
C     is available, that is, A*P = Q*R, where P is a permutation matrix,
C     Q has orthogonal columns, and R is an upper triangular matrix
C     with diagonal elements of nonincreasing magnitude.
C     The routine needs the full upper triangle of R, the permutation
C     matrix P, and the first n components of Q'*b (' denotes the
C     transpose). On output, MD03BB also provides an upper triangular
C     matrix S such that
C
C           P'*(A'*A + PAR*D*D)*P = S'*S .
C
C     Matrix S is used in the solution process.
C
C     This routine is an interface to SLICOT Library routine MD03BY,
C     for solving standard nonlinear least squares problems using SLICOT
C     routine MD03BD.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COND    CHARACTER*1
C             Specifies whether the condition of the matrices R and S
C             should be estimated, as follows:
C             = 'E' :  use incremental condition estimation for R and S;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of R and S for zero values;
C             = 'U' :  use the rank already stored in RANKS (for R).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix R.  N >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters describing the structure of the
C             matrix R. IPAR and LIPAR are not used by this routine,
C             but are provided for compatibility with SLICOT Library
C             routine MD03BD.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 0.
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR, N)
C             On entry, the leading N-by-N upper triangular part of this
C             array must contain the upper triangular matrix R.
C             On exit, the full upper triangle is unaltered, and the
C             strict lower triangle contains the strict upper triangle
C             (transposed) of the upper triangular matrix S.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     IPVT    (input) INTEGER array, dimension (N)
C             This array must define the permutation matrix P such that
C             A*P = Q*R. Column j of P is column IPVT(j) of the identity
C             matrix.
C
C     DIAG    (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the diagonal elements of the
C             matrix D.  DIAG(I) <> 0, I = 1,...,N.
C
C     QTB     (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the first n elements of the
C             vector Q'*b.
C
C     DELTA   (input) DOUBLE PRECISION
C             An upper bound on the Euclidean norm of D*x.  DELTA > 0.
C
C     PAR     (input/output) DOUBLE PRECISION
C             On entry, PAR must contain an initial estimate of the
C             Levenberg-Marquardt parameter.  PAR >= 0.
C             On exit, it contains the final estimate of this parameter.
C
C     RANKS   (input or output) INTEGER array, dimension (1)
C             On entry, if COND = 'U' and N > 0, this array must contain
C             the numerical rank of the matrix R.
C             On exit, this array contains the numerical rank of the
C             matrix S.
C             RANKS is defined as an array for compatibility with SLICOT
C             Library routine MD03BD.
C
C     X       (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the least squares solution of the
C             system A*x = b, sqrt(PAR)*D*x = 0.
C
C     RX      (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the matrix-vector product -R*P'*x.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             rank of the matrices R and S. If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             the reciprocal condition number;  a (sub)matrix whose
C             estimated condition number is less than 1/TOL is
C             considered to be of full rank.  If the user sets TOL <= 0,
C             then an implicitly computed, default tolerance, defined by
C             TOLDEF = N*EPS,  is used instead, where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH).
C             This parameter is not relevant if COND = 'U' or 'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, the first N elements of this array contain the
C             diagonal elements of the upper triangular matrix S.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 4*N, if COND =  'E';
C             LDWORK >= 2*N, if COND <> 'E'.
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
C     This routine calls SLICOT Library routine MD03BY to perform the
C     calculations.
C
C     FURTHER COMMENTS
C
C     For efficiency, the arguments are not checked. This is done in
C     the routine MD03BY (except for LIPAR).
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
C     Linear system of equations, matrix operations, plane rotations.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER         COND
      INTEGER           INFO, LDR, LDWORK, LIPAR, N
      DOUBLE PRECISION  DELTA, PAR, TOL
C     .. Array Arguments ..
      INTEGER           IPAR(*), IPVT(*), RANKS(*)
      DOUBLE PRECISION  DIAG(*), DWORK(*), QTB(*), R(LDR,*), RX(*), X(*)
C     .. External Subroutines ..
      EXTERNAL          MD03BY
C     ..
C     .. Executable Statements ..
C
      CALL MD03BY( COND, N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR,
     $             RANKS(1), X, RX, TOL, DWORK, LDWORK, INFO )
      RETURN
C
C *** Last line of MD03BB ***
      END
