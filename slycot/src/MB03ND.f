      INTEGER FUNCTION MB03ND( N, THETA, Q2, E2, PIVMIN, INFO )
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
C     To find the number of singular values of the bidiagonal matrix
C
C              |q(1) e(1)  .    ...    0   |
C              | 0   q(2) e(2)         .   |
C          J = | .                     .   |
C              | .                   e(N-1)|
C              | 0   ...     ...   0  q(N) |
C
C     which are less than or equal to a given bound THETA.
C
C     This routine is intended to be called only by other SLICOT
C     routines.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the bidiagonal matrix J.  N >= 0.
C
C     THETA   (input) DOUBLE PRECISION
C             Given bound.
C             Note: If THETA < 0.0 on entry, then MB03ND is set to 0
C                   as the singular values of J are non-negative.
C
C     Q2      (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the squares of the diagonal
C             elements q(1),q(2),...,q(N) of the bidiagonal matrix J.
C             That is, Q2(i) = J(i,i)**2 for i = 1,2,...,N.
C
C     E2      (input) DOUBLE PRECISION array, dimension (N-1)
C             This array must contain the squares of the superdiagonal
C             elements e(1),e(2),...,e(N-1) of the bidiagonal matrix J.
C             That is, E2(k) = J(k,k+1)**2 for k = 1,2,...,N-1.
C
C     PIVMIN  (input) DOUBLE PRECISION
C             The minimum absolute value of a "pivot" in the Sturm
C             sequence loop.
C             PIVMIN >= max( max( |q(i)|, |e(k)| )**2*sf_min, sf_min ),
C             where i = 1,2,...,N, k = 1,2,...,N-1, and sf_min is at
C             least the smallest number that can divide one without
C             overflow (see LAPACK Library routine DLAMCH).
C             Note that this condition is not checked by the routine.
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
C     The computation of the number of singular values s(i) of J which
C     are less than or equal to THETA is based on applying Sylvester's
C     Law of Inertia, or equivalently, Sturm sequences [1,p.52] to the
C     unreduced symmetric tridiagonal matrices associated with J as
C     follows. Let T be the following 2N-by-2N symmetric matrix
C     associated with J:
C
C               | 0   J'|
C          T =  |       |.
C               | J   0 |
C
C     (The eigenvalues of T are given by s(1),s(2),...,s(N),-s(1),-s(2),
C     ...,-s(N)). Then, by permuting the rows and columns of T into the
C     order 1, N+1, 2, N+2, ..., N, 2N it follows that T is orthogonally
C     similar to the tridiagonal matrix T" with zeros on its diagonal
C     and q(1), e(1), q(2), e(2), ..., e(N-1), q(N) on its offdiagonals
C     [3,4]. If q(1),q(2),...,q(N) and e(1),e(2),...,e(N-1) are nonzero,
C     Sylvester's Law of Inertia may be applied directly to T".
C     Otherwise, T" is block diagonal and each diagonal block (which is
C     then unreduced) must be analysed separately by applying
C     Sylvester's Law of Inertia.
C
C     REFERENCES
C
C     [1] Parlett, B.N.
C         The Symmetric Eigenvalue Problem.
C         Prentice Hall, Englewood Cliffs, New Jersey, 1980.
C
C     [2] Demmel, J. and Kahan, W.
C         Computing Small Singular Values of Bidiagonal Matrices with
C         Guaranteed High Relative Accuracy.
C         Technical Report, Courant Inst., New York, March 1988.
C
C     [3] Van Huffel, S. and Vandewalle, J.
C         The Partial Total Least-Squares Algorithm.
C         J. Comput. and Appl. Math., 21, pp. 333-341, 1988.
C
C     [4] Golub, G.H. and Kahan, W.
C         Calculating the Singular Values and Pseudo-inverse of a
C         Matrix.
C         SIAM J. Numer. Anal., Ser. B, 2, pp. 205-224, 1965.
C
C     [5] Demmel, J.W., Dhillon, I. and Ren, H.
C         On the Correctness of Parallel Bisection in Floating Point.
C         Computer Science Division Technical Report UCB//CSD-94-805,
C         University of California, Berkeley, CA 94720, March 1994.
C
C     NUMERICAL ASPECTS
C
C     The singular values s(i) could also be obtained with the use of
C     the symmetric tridiagonal matrix T = J'J, whose eigenvalues are
C     the squared singular values of J [4,p.213]. However, the method
C     actually used by the routine is more accurate and equally
C     efficient (see [2]).
C
C     To avoid overflow, matrix J should be scaled so that its largest
C     element is no greater than  overflow**(1/2) * underflow**(1/4)
C     in absolute value (and not much smaller than that, for maximal
C     accuracy).
C
C     With respect to accuracy the following condition holds (see [2]):
C
C     If the established value is denoted by p, then at least p
C     singular values of J are less than or equal to
C     THETA/(1 - (3 x N - 1.5) x EPS) and no more than p singular values
C     are less than or equal to
C     THETA x (1 - (6 x N-2) x EPS)/(1 - (3 x N - 1.5) x EPS).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB03BD by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     July 10, 1997.
C
C     KEYWORDS
C
C     Bidiagonal matrix, singular values.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, N
      DOUBLE PRECISION  PIVMIN, THETA
C     .. Array Arguments ..
      DOUBLE PRECISION  E2(*), Q2(*)
C     .. Local Scalars ..
      INTEGER           J, NUMEIG
      DOUBLE PRECISION  R, T
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
C     Test the input scalar arguments.  PIVMIN is not checked.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
C
C        Error return.
C
         MB03ND = ZERO
         CALL XERBLA( 'MB03ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. THETA.LT.ZERO ) THEN
         MB03ND = 0
         RETURN
      END IF
C
      NUMEIG = N
      T = -THETA
      R = T
      IF ( ABS( R ).LT.PIVMIN ) R = -PIVMIN
C
      DO 20 J = 1, N - 1
         R = T - Q2(J)/R
         IF ( ABS( R ).LT.PIVMIN ) R = -PIVMIN
         IF ( R.GT.ZERO ) NUMEIG = NUMEIG - 1
         R = T - E2(J)/R
         IF ( ABS( R ).LT.PIVMIN ) R = -PIVMIN
         IF ( R.GT.ZERO ) NUMEIG = NUMEIG - 1
   20 CONTINUE
C
      R = T - Q2(N)/R
      IF ( ABS( R ).LT.PIVMIN ) R = -PIVMIN
      IF ( R.GT.ZERO ) NUMEIG = NUMEIG - 1
      MB03ND = NUMEIG
C
      RETURN
C *** Last line of MB03ND ***
      END
