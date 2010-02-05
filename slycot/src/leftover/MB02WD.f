      SUBROUTINE MB02WD( FORM, F, N, IPAR, LIPAR, DPAR, LDPAR, ITMAX,
     $                   A, LDA, B, INCB, X, INCX, TOL, DWORK, LDWORK,
     $                   IWARN, INFO )
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
C     To solve the system of linear equations Ax = b, with A symmetric,
C     positive definite, or, in the implicit form, f(A, x) = b, where
C     y = f(A, x) is a symmetric positive definite linear mapping
C     from x to y, using the conjugate gradient (CG) algorithm without
C     preconditioning.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     FORM     CHARACTER*1
C              Specifies the form of the system of equations, as
C              follows:
C              = 'U' :  Ax = b, the upper triagular part of A is used;
C              = 'L' :  Ax = b, the lower triagular part of A is used;
C              = 'F' :  the implicit, function form, f(A, x) = b.
C
C     Function Parameters
C
C     F       EXTERNAL
C             If FORM = 'F', then F is a subroutine which calculates the
C             value of f(A, x), for given A and x.
C             If FORM <> 'F', then F is not called.
C
C             F must have the following interface:
C
C             SUBROUTINE F( N, IPAR, LIPAR, DPAR, LDPAR, A, LDA, X,
C            $              INCX, DWORK, LDWORK, INFO )
C
C             where
C
C             N       (input) INTEGER
C                     The dimension of the vector x.  N >= 0.
C
C             IPAR    (input) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the matrix A.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 0.
C
C             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C                     The real parameters needed for solving the
C                     problem.
C
C             LDPAR   (input) INTEGER
C                     The length of the array DPAR.  LDPAR >= 0.
C
C             A       (input) DOUBLE PRECISION array, dimension
C                     (LDA, NC), where NC is the number of columns.
C                     The leading NR-by-NC part of this array must
C                     contain the (compressed) representation of the
C                     matrix A, where NR is the number of rows of A
C                     (function of IPAR entries).
C
C             LDA     (input) INTEGER
C                     The leading dimension of the array A.
C                     LDA >= MAX(1,NR).
C
C             X       (input/output) DOUBLE PRECISION array, dimension
C                     (1+(N-1)*INCX)
C                     On entry, this incremented array must contain the
C                     vector x.
C                     On exit, this incremented array contains the value
C                     of the function f, y = f(A, x).
C
C             INCX    (input) INTEGER
C                     The increment for the elements of X.  INCX > 0.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine F.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine F).
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input scalar argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine F. The LAPACK Library routine XERBLA
C                     should be used in conjunction with negative INFO.
C                     INFO must be zero if the subroutine finished
C                     successfully.
C
C             Parameters marked with "(input)" must not be changed.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the vector x.  N >= 0.
C             If FORM = 'U' or FORM = 'L', N is also the number of rows
C             and columns of the matrix A.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             If FORM = 'F', the integer parameters describing the
C             structure of the matrix A.
C             This parameter is ignored if FORM = 'U' or FORM = 'L'.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 0.
C
C     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C             If FORM = 'F', the real parameters needed for solving
C             the problem.
C             This parameter is ignored if FORM = 'U' or FORM = 'L'.
C
C     LDPAR   (input) INTEGER
C             The length of the array DPAR.  LDPAR >= 0.
C
C     ITMAX   (input) INTEGER
C             The maximal number of iterations to do.  ITMAX >= 0.
C
C     A       (input) DOUBLE PRECISION array,
C                     dimension (LDA, NC), if FORM = 'F',
C                     dimension (LDA, N),  otherwise.
C             If FORM = 'F', the leading NR-by-NC part of this array
C             must contain the (compressed) representation of the
C             matrix A, where NR and NC are the number of rows and
C             columns, respectively, of the matrix A. The array A is
C             not referenced by this routine itself, except in the
C             calls to the routine F.
C             If FORM <> 'F', the leading N-by-N part of this array
C             must contain the matrix A, assumed to be symmetric;
C             only the triangular part specified by FORM is referenced.
C
C     LDA     (input) INTEGER
C             The leading dimension of array A.
C             LDA >= MAX(1,NR), if FORM = 'F';
C             LDA >= MAX(1,N),  if FORM = 'U' or FORM = 'L'.
C
C     B       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCB)
C             The incremented vector b.
C
C     INCB    (input) INTEGER
C             The increment for the elements of B.  INCB > 0.
C
C     X       (input/output) DOUBLE PRECISION array, dimension
C             (1+(N-1)*INCX)
C             On entry, this incremented array must contain an initial
C             approximation of the solution. If an approximation is not
C             known, setting all elements of x to zero is recommended.
C             On exit, this incremented array contains the computed
C             solution x of the system of linear equations.
C
C     INCX    (input) INTEGER
C             The increment for the elements of X.  INCX > 0.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If TOL > 0, absolute tolerance for the iterative process.
C             The algorithm will stop if || Ax - b ||_2 <= TOL. Since
C             it is advisable to use a relative tolerance, say TOLER,
C             TOL should be chosen as TOLER*|| b ||_2.
C             If TOL <= 0, a default relative tolerance,
C             TOLDEF = N*EPS*|| b ||_2,  is used, where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the number of
C             iterations performed and DWORK(2) returns the remaining
C             residual, || Ax - b ||_2.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(2,3*N + DWORK(F)),  if FORM = 'F',
C                       where DWORK(F) is the workspace needed by F;
C             LDWORK >= MAX(2,3*N),       if FORM = 'U' or FORM = 'L'.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  the algorithm finished after ITMAX > 0 iterations,
C                   without achieving the desired precision TOL;
C             = 2:  ITMAX is zero; in this case, DWORK(2) is not set.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, then F returned with INFO = i.
C
C     METHOD
C
C     The following CG iteration is used for solving Ax = b:
C
C     Start: q(0) = r(0) = Ax - b
C
C                   < q(k),  r(k) >
C     ALPHA(k) = - ----------------
C                   < q(k), Aq(k) >
C     x(k+1)   = x(k) - ALPHA(k) * q(k)
C     r(k+1)   = r(k) - ALPHA(k) * Aq(k)
C                 < r(k+1), r(k+1) >
C     BETA(k)  = --------------------
C                 < r(k)  , r(k)   >
C     q(k+1)   = r(k+1) + BETA(k) * q(k)
C
C     where <.,.> denotes the scalar product.
C
C     REFERENCES
C
C     [1] Golub, G.H. and van Loan, C.F.
C         Matrix Computations. Third Edition.
C         M. D. Johns Hopkins University Press, Baltimore, pp. 520-528,
C         1996.
C
C     [2] Luenberger, G.
C         Introduction to Linear and Nonlinear Programming.
C         Addison-Wesley, Reading, MA, p.187, York, 1973.
C
C     NUMERICAL ASPECTS
C
C     Since the residuals are orthogonal in the scalar product
C     <x, y> = y'Ax, the algorithm is theoretically finite. But rounding
C     errors cause a loss of orthogonality, so a finite termination
C     cannot be guaranteed. However, one can prove [2] that
C
C        || x-x_k ||_A := sqrt( (x-x_k)' * A * (x-x_k) )
C
C                                             sqrt( kappa_2(A) ) - 1
C                      <=  2 || x-x_0 ||_A * ------------------------ ,
C                                             sqrt( kappa_2(A) ) + 1
C
C     where kappa_2 is the condition number.
C
C     The approximate number of floating point operations is
C        (k*(N**2 + 15*N) + N**2 + 3*N)/2, if FORM <> 'F',
C        k*(f + 7*N) + f,                  if FORM =  'F',
C     where k is the number of CG iterations performed, and f is the
C     number of floating point operations required by the subroutine F.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     March, 2002.
C
C     KEYWORDS
C
C     Conjugate gradients, convergence, linear system of equations,
C     matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         FORM
      INTEGER           INCB, INCX, INFO, ITMAX, IWARN, LDA, LDPAR,
     $                  LDWORK, LIPAR, N
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(*), DPAR(*), DWORK(*), X(*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BETA, RES, RESOLD, TOLDEF
      INTEGER           AQ, DWLEFT, K, R
      LOGICAL           MAT
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DDOT, DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, DSYMV, F, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      MAT = LSAME( FORM, 'U' ) .OR. LSAME( FORM, 'L' )
C
C     Check the scalar input parameters.
C
      IWARN = 0
      INFO  = 0
      IF( .NOT.( MAT .OR. LSAME( FORM, 'F' ) ) ) THEN
         INFO = -1
      ELSEIF ( N.LT.0 ) THEN
         INFO = -3
      ELSEIF ( .NOT. MAT .AND. LIPAR.LT.0 ) THEN
         INFO = -5
      ELSEIF ( .NOT. MAT .AND. LDPAR.LT.0 ) THEN
         INFO = -7
      ELSEIF ( ITMAX.LT.0 ) THEN
         INFO = -8
      ELSEIF ( LDA.LT.1 .OR. ( MAT .AND. LDA.LT.N ) ) THEN
         INFO = -10
      ELSEIF ( INCB.LE.0 ) THEN
         INFO = -12
      ELSEIF ( INCX.LE.0 ) THEN
         INFO = -14
      ELSEIF ( LDWORK.LT.MAX( 2, 3*N ) ) THEN
         INFO = -17
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02WD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ZERO
         DWORK(2) = ZERO
         RETURN
      ENDIF
C
      IF ( ITMAX.EQ.0 ) THEN
         DWORK(1) = ZERO
         IWARN = 2
         RETURN
      ENDIF
C
C     Set default tolerance, if needed.
C
      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO )
     $   TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )*DNRM2( N, B, INCB )
C
C     Initialize local variables.
C
      K = 0
C
C     Vector q is stored in DWORK(1), A*q or f(A, q) in DWORK(AQ),
C     and r in DWORK(R). The workspace for F starts in DWORK(DWLEFT).
C
      AQ     = N + 1
      R      = N + AQ
      DWLEFT = N + R
C
C     Prepare the first iteration, initialize r and q.
C
      IF ( MAT ) THEN
         CALL DCOPY( N, B, INCB, DWORK(R), 1 )
         CALL DSYMV( FORM, N, ONE, A, LDA, X, INCX, -ONE, DWORK(R), 1 )
      ELSE
         CALL DCOPY( N, X, INCX, DWORK(R), 1 )
         CALL F( N, IPAR, LIPAR, DPAR, LDPAR, A, LDA, DWORK(R), 1,
     $           DWORK(DWLEFT), LDWORK-DWLEFT+1, INFO )
         IF ( INFO.NE.0 )
     $      RETURN
         CALL DAXPY( N, -ONE, B, INCB, DWORK(R), 1 )
      ENDIF
      CALL DCOPY( N, DWORK(R), 1, DWORK, 1 )
C
      RES = DNRM2( N, DWORK(R), 1 )
C
C     Do nothing if x is already the solution.
C
      IF ( RES.LE.TOLDEF ) GOTO 20
C
C     Begin of the iteration loop.
C
C     WHILE ( RES.GT.TOLDEF .AND. K.LE.ITMAX ) DO
   10 CONTINUE
C
C        Calculate A*q or f(A, q).
C
         IF ( MAT ) THEN
            CALL DSYMV( FORM, N, ONE, A, LDA, DWORK, 1, ZERO, DWORK(AQ),
     $                  1 )
         ELSE
            CALL DCOPY( N, DWORK, 1, DWORK(AQ), 1 )
            CALL F( N, IPAR, LIPAR, DPAR, LDPAR, A, LDA, DWORK(AQ), 1,
     $              DWORK(DWLEFT), LDWORK-DWLEFT+1, INFO )
            IF ( INFO.NE.0 )
     $         RETURN
         ENDIF
C
C        Calculate ALPHA(k).
C
         ALPHA = DDOT( N, DWORK, 1, DWORK(R),  1 ) /
     $           DDOT( N, DWORK, 1, DWORK(AQ), 1 )
C
C        x(k+1) = x(k) - ALPHA(k)*q(k).
C
         CALL DAXPY( N, -ALPHA, DWORK, 1, X, INCX )
C
C        r(k+1) = r(k) - ALPHA(k)*(A*q(k)).
C
         CALL DAXPY( N, -ALPHA, DWORK(AQ), 1, DWORK(R), 1 )
C
C        Save RES and calculate a new RES.
C
         RESOLD = RES
         RES = DNRM2( N, DWORK(R), 1 )
C
C        Exit if tolerance is reached.
C
         IF ( RES.LE.TOLDEF ) GOTO 20
C
C        Calculate BETA(k).
C
         BETA = ( RES/RESOLD )**2
C
C        q(k+1) = r(k+1) + BETA(k)*q(k).
C
         CALL DSCAL( N, BETA, DWORK, 1 )
         CALL DAXPY( N, ONE,  DWORK(R), 1, DWORK, 1 )
C
C        End of the iteration loop.
C
         K = K + 1
         IF ( K.LT.ITMAX ) GOTO 10
C     END WHILE 10
C
C     Tolerance was not reached!
C
      IWARN = 1
C
   20 CONTINUE
C
      DWORK(1) = K
      DWORK(2) = RES
C
C *** Last line of MB02WD ***
      END
