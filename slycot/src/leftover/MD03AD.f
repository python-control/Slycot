      SUBROUTINE MD03AD( XINIT, ALG, STOR, UPLO, FCN, JPJ, M, N, ITMAX,
     $                   NPRINT, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                   LDPAR2, X, NFEV, NJEV, TOL, CGTOL, DWORK,
     $                   LDWORK, IWARN, INFO )
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
C     To minimize the sum of the squares of m nonlinear functions, e, in
C     n variables, x, by a modification of the Levenberg-Marquardt
C     algorithm, using either a Cholesky-based or a conjugate gradients
C     solver. The user must provide a subroutine FCN which calculates
C     the functions and the Jacobian J (possibly by finite differences),
C     and another subroutine JPJ, which computes either J'*J + par*I
C     (if ALG = 'D'), or (J'*J + par*I)*x (if ALG = 'I'), where par is
C     the Levenberg factor, exploiting the possible structure of the
C     Jacobian matrix. Template implementations of these routines are
C     included in the SLICOT Library.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     XINIT   CHARACTER*1
C             Specifies how the variables x are initialized, as follows:
C             = 'R' :  the array X is initialized to random values; the
C                      entries DWORK(1:4) are used to initialize the
C                      random number generator: the first three values
C                      are converted to integers between 0 and 4095, and
C                      the last one is converted to an odd integer
C                      between 1 and 4095;
C             = 'G' :  the given entries of X are used as initial values
C                      of variables.
C
C     ALG     CHARACTER*1
C             Specifies the algorithm used for solving the linear
C             systems involving a Jacobian matrix J, as follows:
C             = 'D' :  a direct algorithm, which computes the Cholesky
C                      factor of the matrix J'*J + par*I is used;
C             = 'I' :  an iterative Conjugate Gradients algorithm, which
C                      only needs the matrix J, is used.
C             In both cases, matrix J is stored in a compressed form.
C
C     STOR    CHARACTER*1
C             If ALG = 'D', specifies the storage scheme for the
C             symmetric matrix J'*J, as follows:
C             = 'F' :  full storage is used;
C             = 'P' :  packed storage is used.
C             The option STOR = 'F' usually ensures a faster execution.
C             This parameter is not relevant if ALG = 'I'.
C
C     UPLO    CHARACTER*1
C             If ALG = 'D', specifies which part of the matrix J'*J
C             is stored, as follows:
C             = 'U' :  the upper triagular part is stored;
C             = 'L' :  the lower triagular part is stored.
C             The option UPLO = 'U' usually ensures a faster execution.
C             This parameter is not relevant if ALG = 'I'.
C
C     Function Parameters
C
C     FCN     EXTERNAL
C             Subroutine which evaluates the functions and the Jacobian.
C             FCN must be declared in an external statement in the user
C             calling program, and must have the following interface:
C
C             SUBROUTINE FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1,
C            $                DPAR2, LDPAR2, X, NFEVL, E, J, LDJ, JTE,
C            $                DWORK, LDWORK, INFO )
C
C             where
C
C             IFLAG   (input/output) INTEGER
C                     On entry, this parameter must contain a value
C                     defining the computations to be performed:
C                     = 0 :  Optionally, print the current iterate X,
C                            function values E, and Jacobian matrix J,
C                            or other results defined in terms of these
C                            values. See the argument NPRINT of MD03AD.
C                            Do not alter E and J.
C                     = 1 :  Calculate the functions at X and return
C                            this vector in E. Do not alter J.
C                     = 2 :  Calculate the Jacobian at X and return
C                            this matrix in J. Also return J'*e in JTE
C                            and NFEVL (see below). Do not alter E.
C                     = 3 :  Do not compute neither the functions nor
C                            the Jacobian, but return in LDJ and
C                            IPAR/DPAR1,DPAR2 (some of) the integer/real
C                            parameters needed.
C                     On exit, the value of this parameter should not be
C                     changed by FCN unless the user wants to terminate
C                     execution of MD03AD, in which case IFLAG must be
C                     set to a negative integer.
C
C             M       (input) INTEGER
C                     The number of functions.  M >= 0.
C
C             N       (input) INTEGER
C                     The number of variables.  M >= N >= 0.
C
C             IPAR    (input/output) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the Jacobian matrix or needed for problem solving.
C                     IPAR is an input parameter, except for IFLAG = 3
C                     on entry, when it is also an output parameter.
C                     On exit, if IFLAG = 3, IPAR(1) contains the length
C                     of the array J, for storing the Jacobian matrix,
C                     and the entries IPAR(2:5) contain the workspace
C                     required by FCN for IFLAG = 1, FCN for IFLAG = 2,
C                     JPJ for ALG = 'D', and JPJ for ALG = 'I',
C                     respectively.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 5.
C
C             DPAR1   (input/output) DOUBLE PRECISION array, dimension
C                     (LDPAR1,*) or (LDPAR1)
C                     A first set of real parameters needed for
C                     describing or solving the problem.
C                     DPAR1 can also be used as an additional array for
C                     intermediate results when computing the functions
C                     or the Jacobian. For control problems, DPAR1 could
C                     store the input trajectory of a system.
C
C             LDPAR1  (input) INTEGER
C                     The leading dimension or the length of the array
C                     DPAR1, as convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1,
C                     if leading dimension.)
C
C             DPAR2   (input/output) DOUBLE PRECISION array, dimension
C                     (LDPAR2,*) or (LDPAR2)
C                     A second set of real parameters needed for
C                     describing or solving the problem.
C                     DPAR2 can also be used as an additional array for
C                     intermediate results when computing the functions
C                     or the Jacobian. For control problems, DPAR2 could
C                     store the output trajectory of a system.
C
C             LDPAR2  (input) INTEGER
C                     The leading dimension or the length of the array
C                     DPAR2, as convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1,
C                     if leading dimension.)
C
C             X       (input) DOUBLE PRECISION array, dimension (N)
C                     This array must contain the value of the
C                     variables x where the functions or the Jacobian
C                     must be evaluated.
C
C             NFEVL   (input/output) INTEGER
C                     The number of function evaluations needed to
C                     compute the Jacobian by a finite difference
C                     approximation.
C                     NFEVL is an input parameter if IFLAG = 0, or an
C                     output parameter if IFLAG = 2. If the Jacobian is
C                     computed analytically, NFEVL should be set to a
C                     non-positive value.
C
C             E       (input/output) DOUBLE PRECISION array,
C                     dimension (M)
C                     This array contains the value of the (error)
C                     functions e evaluated at X.
C                     E is an input parameter if IFLAG = 0 or 2, or an
C                     output parameter if IFLAG = 1.
C
C             J       (input/output) DOUBLE PRECISION array, dimension
C                     (LDJ,NC), where NC is the number of columns
C                     needed.
C                     This array contains a possibly compressed
C                     representation of the Jacobian matrix evaluated
C                     at X. If full Jacobian is stored, then NC = N.
C                     J is an input parameter if IFLAG = 0, or an output
C                     parameter if IFLAG = 2.
C
C             LDJ     (input/output) INTEGER
C                     The leading dimension of array J.  LDJ >= 1.
C                     LDJ is essentially used inside the routines FCN
C                     and JPJ.
C                     LDJ is an input parameter, except for IFLAG = 3
C                     on entry, when it is an output parameter.
C                     It is assumed in MD03AD that LDJ is not larger
C                     than needed.
C
C             JTE     (output) DOUBLE PRECISION array, dimension (N)
C                     If IFLAG = 2, the matrix-vector product J'*e.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine FCN.
C                     On exit, if INFO = 0, DWORK(1) returns the optimal
C                     value of LDWORK.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine FCN).  LDWORK >= 1.
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input (scalar) argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine FCN. The LAPACK Library routine XERBLA
C                     should be used in conjunction with negative INFO.
C                     INFO must be zero if the subroutine finished
C                     successfully.
C
C             Parameters marked with "(input)" must not be changed.
C
C     JPJ     EXTERNAL
C             Subroutine which computes J'*J + par*I, if ALG = 'D', and
C             J'*J*x + par*x, if ALG = 'I', where J is the Jacobian as
C             described above.
C
C             JPJ must have the following interface:
C
C             SUBROUTINE JPJ( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR,
C            $                J, LDJ, JTJ, LDJTJ, DWORK, LDWORK, INFO )
C
C             if ALG = 'D', and
C
C             SUBROUTINE JPJ( N, IPAR, LIPAR, DPAR, LDPAR, J, LDJ, X,
C            $                INCX, DWORK, LDWORK, INFO )
C
C             if ALG = 'I', where
C
C             STOR    (input) CHARACTER*1
C                     Specifies the storage scheme for the symmetric
C                     matrix J'*J, as follows:
C                     = 'F' :  full storage is used;
C                     = 'P' :  packed storage is used.
C
C             UPLO    (input) CHARACTER*1
C                     Specifies which part of the matrix J'*J is stored,
C                     as follows:
C                     = 'U' :  the upper triagular part is stored;
C                     = 'L' :  the lower triagular part is stored.
C
C             N       (input) INTEGER
C                     The number of columns of the matrix J.  N >= 0.
C
C             IPAR    (input) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the Jacobian matrix.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 0.
C
C             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C                     DPAR(1) must contain an initial estimate of the
C                     Levenberg-Marquardt parameter, par.  DPAR(1) >= 0.
C
C             LDPAR   (input) INTEGER
C                     The length of the array DPAR.  LDPAR >= 1.
C
C             J       (input) DOUBLE PRECISION array, dimension
C                     (LDJ, NC), where NC is the number of columns.
C                     The leading NR-by-NC part of this array must
C                     contain the (compressed) representation of the
C                     Jacobian matrix J, where NR is the number of rows
C                     of J (function of IPAR entries).
C
C             LDJ     (input) INTEGER
C                     The leading dimension of array J.
C                     LDJ >= MAX(1,NR).
C
C             JTJ     (output) DOUBLE PRECISION array,
C                              dimension (LDJTJ,N),    if STOR = 'F',
C                              dimension (N*(N+1)/2),  if STOR = 'P'.
C                     The leading N-by-N (if STOR = 'F'), or N*(N+1)/2
C                     (if STOR = 'P') part of this array contains the
C                     upper or lower triangle of the matrix J'*J+par*I,
C                     depending on UPLO = 'U', or UPLO = 'L',
C                     respectively, stored either as a two-dimensional,
C                     or one-dimensional array, depending on STOR.
C
C             LDJTJ   (input) INTEGER
C                     The leading dimension of the array JTJ.
C                     LDJTJ >= MAX(1,N), if STOR = 'F'.
C                     LDJTJ >= 1,        if STOR = 'P'.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine JPJ.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine JPJ).
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input (scalar) argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine JPJ. The LAPACK Library routine XERBLA
C                     should be used in conjunction with negative INFO
C                     values. INFO must be zero if the subroutine
C                     finished successfully.
C
C             If ALG = 'I', the parameters in common with those for
C             ALG = 'D', have the same meaning, and the additional
C             parameters are:
C
C             X       (input/output) DOUBLE PRECISION array, dimension
C                     (1+(N-1)*INCX)
C                     On entry, this incremented array must contain the
C                     vector x.
C                     On exit, this incremented array contains the value
C                     of the matrix-vector product (J'*J + par)*x.
C
C             INCX    (input) INTEGER
C                     The increment for the elements of X.  INCX > 0.
C
C             Parameters marked with "(input)" must not be changed.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of functions.  M >= 0.
C
C     N       (input) INTEGER
C             The number of variables.  M >= N >= 0.
C
C     ITMAX   (input) INTEGER
C             The maximum number of iterations.  ITMAX >= 0.
C
C     NPRINT  (input) INTEGER
C             This parameter enables controlled printing of iterates if
C             it is positive. In this case, FCN is called with IFLAG = 0
C             at the beginning of the first iteration and every NPRINT
C             iterations thereafter and immediately prior to return,
C             with X, E, and J available for printing. If NPRINT is not
C             positive, no special calls of FCN with IFLAG = 0 are made.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters needed, for instance, for
C             describing the structure of the Jacobian matrix, which
C             are handed over to the routines FCN and JPJ.
C             The first five entries of this array are modified
C             internally by a call to FCN (with IFLAG = 3), but are
C             restored on exit.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 5.
C
C     DPAR1   (input/output) DOUBLE PRECISION array, dimension
C             (LDPAR1,*) or (LDPAR1)
C             A first set of real parameters needed for describing or
C             solving the problem. This argument is not used by MD03AD
C             routine, but it is passed to the routine FCN.
C
C     LDPAR1  (input) INTEGER
C             The leading dimension or the length of the array DPAR1, as
C             convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, if leading
C             dimension.)
C
C     DPAR2   (input/output) DOUBLE PRECISION array, dimension
C             (LDPAR2,*) or (LDPAR2)
C             A second set of real parameters needed for describing or
C             solving the problem. This argument is not used by MD03AD
C             routine, but it is passed to the routine FCN.
C
C     LDPAR2  (input) INTEGER
C             The leading dimension or the length of the array DPAR2, as
C             convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1, if leading
C             dimension.)
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, if XINIT = 'G', this array must contain the
C             vector of initial variables x to be optimized.
C             If XINIT = 'R', this array need not be set before entry,
C             and random values will be used to initialize x.
C             On exit, if INFO = 0, this array contains the vector of
C             values that (approximately) minimize the sum of squares of
C             error functions. The values returned in IWARN and
C             DWORK(1:5) give details on the iterative process.
C
C     NFEV    (output) INTEGER
C             The number of calls to FCN with IFLAG = 1. If FCN is
C             properly implemented, this includes the function
C             evaluations needed for finite difference approximation
C             of the Jacobian.
C
C     NJEV    (output) INTEGER
C             The number of calls to FCN with IFLAG = 2.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If TOL >= 0, the tolerance which measures the relative
C             error desired in the sum of squares. Termination occurs
C             when the actual relative reduction in the sum of squares
C             is at most TOL. If the user sets  TOL < 0, then  SQRT(EPS)
C             is used instead TOL, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH).
C
C     CGTOL   DOUBLE PRECISION
C             If ALG = 'I' and CGTOL > 0, the tolerance which measures
C             the relative residual of the solutions computed by the
C             conjugate gradients (CG) algorithm. Termination of a
C             CG process occurs when the relative residual is at
C             most CGTOL. If the user sets  CGTOL <= 0, then  SQRT(EPS)
C             is used instead CGTOL.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, DWORK(2) returns the residual error norm (the
C             sum of squares), DWORK(3) returns the number of iterations
C             performed, DWORK(4) returns the total number of conjugate
C             gradients iterations performed (zero, if ALG = 'D'), and
C             DWORK(5) returns the final Levenberg factor.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= max( 5, M + 2*N + size(J) +
C                            max( DW( FCN|IFLAG = 1 ) + N,
C                                 DW( FCN|IFLAG = 2 ),
C                                 DW( sol ) ) ),
C             where size(J) is the size of the Jacobian (provided by FCN
C             in IPAR(1), for IFLAG = 3), DW( f ) is the workspace
C             needed by the routine f, where f is FCN or JPJ (provided
C             by FCN in IPAR(2:5), for IFLAG = 3), and DW( sol ) is the
C             workspace needed for solving linear systems,
C             DW( sol ) = N*N + DW( JPJ ),  if ALG = 'D', STOR = 'F';
C             DW( sol ) = N*(N+1)/2 + DW( JPJ ),
C                                           if ALG = 'D', STOR = 'P';
C             DW( sol ) = 3*N + DW( JPJ ),  if ALG = 'I'.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             < 0:  the user set IFLAG = IWARN in the subroutine FCN;
C             = 0:  no warning;
C             = 1:  if the iterative process did not converge in ITMAX
C                   iterations with tolerance TOL;
C             = 2:  if ALG = 'I', and in one or more iterations of the
C                   Levenberg-Marquardt algorithm, the conjugate
C                   gradient algorithm did not finish after 3*N
C                   iterations, with the accuracy required in the
C                   call;
C             = 3:  the cosine of the angle between e and any column of
C                   the Jacobian is at most FACTOR*EPS in absolute
C                   value, where FACTOR = 100 is defined in a PARAMETER
C                   statement;
C             = 4:  TOL is too small: no further reduction in the sum
C                   of squares is possible.
C                   In all these cases, DWORK(1:5) are set as described
C                   above.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  user-defined routine FCN returned with INFO <> 0
C                   for IFLAG = 1;
C             = 2:  user-defined routine FCN returned with INFO <> 0
C                   for IFLAG = 2;
C             = 3:  SLICOT Library routine MB02XD, if ALG = 'D', or
C                   SLICOT Library routine MB02WD, if ALG = 'I' (or
C                   user-defined routine JPJ), returned with INFO <> 0.
C
C     METHOD
C
C     If XINIT = 'R', the initial value for X is set to a vector of
C     pseudo-random values uniformly distributed in [-1,1].
C
C     The Levenberg-Marquardt algorithm (described in [1]) is used for
C     optimizing the parameters. This algorithm needs the Jacobian
C     matrix J, which is provided by the subroutine FCN. The algorithm
C     tries to update x by the formula
C
C         x = x - p,
C
C     using the solution of the system of linear equations
C
C         (J'*J + PAR*I)*p = J'*e,
C
C     where I is the identity matrix, and e the error function vector.
C     The Levenberg factor PAR is decreased after each successfull step
C     and increased in the other case.
C
C     If ALG = 'D', a direct method, which evaluates the matrix product
C     J'*J + par*I and then factors it using Cholesky algorithm,
C     implemented in the SLICOT Libray routine MB02XD, is used for
C     solving the linear system above.
C
C     If ALG = 'I', the Conjugate Gradients method, described in [2],
C     and implemented in the SLICOT Libray routine MB02WD, is used for
C     solving the linear system above. The main advantage of this method
C     is that in most cases the solution of the system can be computed
C     in less time than the time needed to compute the matrix J'*J
C     This is, however, problem dependent.
C
C     REFERENCES
C
C     [1] Kelley, C.T.
C         Iterative Methods for Optimization.
C         Society for Industrial and Applied Mathematics (SIAM),
C         Philadelphia (Pa.), 1999.
C
C     [2] Golub, G.H. and van Loan, C.F.
C         Matrix Computations. Third Edition.
C         M. D. Johns Hopkins University Press, Baltimore, pp. 520-528,
C         1996.
C
C     [3] More, J.J.
C         The Levenberg-Marquardt algorithm: implementation and theory.
C         In Watson, G.A. (Ed.), Numerical Analysis, Lecture Notes in
C         Mathematics, vol. 630, Springer-Verlag, Berlin, Heidelberg
C         and New York, pp. 105-116, 1978.
C
C     NUMERICAL ASPECTS
C
C     The Levenberg-Marquardt algorithm described in [3] is scaling
C     invariant and globally convergent to (maybe local) minima.
C     According to [1], the convergence rate near a local minimum is
C     quadratic, if the Jacobian is computed analytically, and linear,
C     if the Jacobian is computed numerically.
C
C     Whether or not the direct algorithm is faster than the iterative
C     Conjugate Gradients algorithm for solving the linear systems
C     involved depends on several factors, including the conditioning
C     of the Jacobian matrix, and the ratio between its dimensions.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Mar. 2002.
C
C     KEYWORDS
C
C     Conjugate gradients, least-squares approximation,
C     Levenberg-Marquardt algorithm, matrix operations, optimization.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, FOUR, FIVE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0,
     $                    FIVE = 5.0D0 )
      DOUBLE PRECISION  FACTOR, MARQF, MINIMP, PARMAX
      PARAMETER         ( FACTOR = 10.0D0**2,    MARQF  = 2.0D0**2,
     $                    MINIMP =  2.0D0**(-3), PARMAX = 1.0D20 )
C     .. Scalar Arguments ..
      CHARACTER         ALG, STOR, UPLO, XINIT
      INTEGER           INFO, ITMAX, IWARN, LDPAR1, LDPAR2, LDWORK,
     $                  LIPAR, M, N, NFEV, NJEV, NPRINT
      DOUBLE PRECISION  CGTOL, TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  DPAR1(LDPAR1,*), DPAR2(LDPAR2,*), DWORK(*), X(*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      LOGICAL           CHOL, FULL, INIT, UPPER
      INTEGER           DWJTJ, E, I, IFLAG, INFOL, ITER, ITERCG, IW1,
     $                  IW2, IWARNL, JAC, JTE, JW1, JW2, JWORK, LDJ,
     $                  LDW, LFCN1, LFCN2, LJTJ, LJTJD, LJTJI, NFEVL,
     $                  SIZEJ, WRKOPT
      DOUBLE PRECISION  ACTRED, BIGNUM, CGTDEF, EPSMCH, FNORM, FNORM1,
     $                  GNORM, GSMIN, PAR, SMLNUM, SQREPS, TOLDEF
C     .. Local Arrays ..
      INTEGER           SEED(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DDOT, DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLABAD, DLARNV, FCN, JPJ, MB02WD, MB02XD,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN, MOD, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INIT  = LSAME( XINIT, 'R' )
      CHOL  = LSAME( ALG,   'D' )
      FULL  = LSAME( STOR,  'F' )
      UPPER = LSAME( UPLO,  'U' )
C
C     Check the scalar input parameters.
C
      IWARN = 0
      INFO  = 0
      IF( .NOT.( INIT .OR. LSAME( XINIT, 'G' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( CHOL .OR. LSAME( ALG, 'I' ) ) ) THEN
         INFO = -2
      ELSEIF ( CHOL .AND. .NOT.( FULL  .OR. LSAME( STOR, 'P' ) ) ) THEN
         INFO = -3
      ELSEIF ( CHOL .AND. .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -4
      ELSEIF ( M.LT.0 ) THEN
         INFO = -7
      ELSEIF ( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -8
      ELSEIF ( ITMAX.LT.0 ) THEN
         INFO = -9
      ELSEIF ( LIPAR.LT.5 ) THEN
         INFO = -12
      ELSEIF( LDPAR1.LT.0 ) THEN
         INFO = -14
      ELSEIF( LDPAR2.LT.0 ) THEN
         INFO = -16
      ELSEIF ( LDWORK.LT.5 ) THEN
         INFO = -23
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MD03AD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      NFEV = 0
      NJEV = 0
      IF ( MIN( N, ITMAX ).EQ.0 ) THEN
         DWORK(1) = FIVE
         DWORK(2) = ZERO
         DWORK(3) = ZERO
         DWORK(4) = ZERO
         DWORK(5) = ZERO
         RETURN
      ENDIF
C
C     Call FCN to get the size of the array J, for storing the Jacobian
C     matrix, the leading dimension LDJ and the workspace required
C     by FCN for IFLAG = 1 and IFLAG = 2, and JPJ. The entries
C     DWORK(1:4) should not be modified by the special call of FCN
C     below, if XINIT = 'R' and the values in DWORK(1:4) are explicitly
C     desired for initialization of the random number generator.
C
      IFLAG = 3
      IW1   = IPAR(1)
      IW2   = IPAR(2)
      JW1   = IPAR(3)
      JW2   = IPAR(4)
      LJTJ  = IPAR(5)
C
      CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2, LDPAR2,
     $          X, NFEVL, DWORK, DWORK, LDJ, DWORK, DWORK, LDWORK,
     $          INFOL )
C
      SIZEJ = IPAR(1)
      LFCN1 = IPAR(2)
      LFCN2 = IPAR(3)
      LJTJD = IPAR(4)
      LJTJI = IPAR(5)
C
      IPAR(1) = IW1
      IPAR(2) = IW2
      IPAR(3) = JW1
      IPAR(4) = JW2
      IPAR(5) = LJTJ
C
C     Define pointers to the array variables stored in DWORK.
C
      JAC = 1
      E   = JAC + SIZEJ
      JTE = E   + M
      IW1 = JTE + N
      IW2 = IW1 + N
      JW1 = IW2
      JW2 = IW2 + N
C
C     Check the workspace length.
C
      JWORK = JW1
      IF ( CHOL ) THEN
         IF ( FULL ) THEN
            LDW = N*N
         ELSE
            LDW = ( N*( N + 1 ) ) / 2
         ENDIF
         DWJTJ = JWORK
         JWORK = DWJTJ + LDW
         LJTJ  = LJTJD
      ELSE
         LDW  = 3*N
         LJTJ = LJTJI
      ENDIF
      IF ( LDWORK.LT.MAX( 5, SIZEJ + M + 2*N +
     $                       MAX( LFCN1 + N, LFCN2, LDW + LJTJ ) ) )
     $      THEN
         INFO = -23
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MD03AD', -INFO )
         RETURN
      ENDIF
C
C     Set default tolerances. SQREPS is the square root of the machine
C     precision, and GSMIN is used in the tests of the gradient norm.
C
      EPSMCH = DLAMCH( 'Epsilon' )
      SQREPS = SQRT( EPSMCH )
      TOLDEF = TOL
      IF ( TOLDEF.LT.ZERO )
     $   TOLDEF = SQREPS
      CGTDEF = CGTOL
      IF ( CGTDEF.LE.ZERO )
     $   CGTDEF = SQREPS
      GSMIN  = FACTOR*EPSMCH
      WRKOPT = 5
C
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
C     Initialization.
C
      IF ( INIT ) THEN
C
C        SEED is the initial state of the random number generator.
C        SEED(4) must be odd.
C
         SEED(1) = MOD(   INT( DWORK(1) ), 4096 )
         SEED(2) = MOD(   INT( DWORK(2) ), 4096 )
         SEED(3) = MOD(   INT( DWORK(3) ), 4096 )
         SEED(4) = MOD( 2*INT( DWORK(4) ) + 1, 4096 )
         CALL DLARNV( 2, SEED, N, X )
      ENDIF
C
C     Evaluate the function at the starting point and calculate
C     its norm.
C     Workspace: need:    SIZEJ + M + 2*N + LFCN1;
C                prefer:  larger.
C
      IFLAG = 1
      CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2, LDPAR2,
     $          X, NFEVL, DWORK(E), DWORK(JAC), LDJ, DWORK(JTE),
     $          DWORK(JW1), LDWORK-JW1+1, INFOL )
C
      IF ( INFOL.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(JW1) ) + JW1 - 1 )
      NFEV   = 1
      FNORM  = DNRM2( M, DWORK(E), 1 )
      ACTRED = ZERO
      ITERCG = 0
      ITER   = 0
      IWARNL = 0
      PAR    = ZERO
      IF ( IFLAG.LT.0 .OR. FNORM.EQ.ZERO )
     $   GO TO 40
C
C     Set the initial vector for the conjugate gradients algorithm.
C
      DWORK(IW1) = ZERO
      CALL DCOPY( N, DWORK(IW1), 0, DWORK(IW1), 1 )
C
C     WHILE ( nonconvergence and ITER < ITMAX ) DO
C
C     Beginning of the outer loop.
C
   10 CONTINUE
C
C        Calculate the Jacobian matrix.
C        Workspace: need:    SIZEJ + M + 2*N + LFCN2;
C                   prefer:  larger.
C
         ITER  = ITER + 1
         IFLAG = 2
         CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $             LDPAR2, X, NFEVL, DWORK(E), DWORK(JAC), LDJ,
     $             DWORK(JTE), DWORK(JW1), LDWORK-JW1+1, INFOL )
C
         IF ( INFOL.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
C
C        Compute the gradient norm.
C
         GNORM = DNRM2( N, DWORK(JTE), 1 )
         IF ( NFEVL.GT.0 )
     $      NFEV = NFEV + NFEVL
         NJEV = NJEV + 1
         IF ( GNORM.LE.GSMIN )
     $      IWARN = 3
         IF ( IWARN.NE.0 )
     $      GO TO 40
         IF ( ITER.EQ.1 ) THEN
            WRKOPT = MAX( WRKOPT, INT( DWORK(JW1) ) + JW1 - 1 )
            PAR    = MIN( GNORM, SQRT( PARMAX ) )
         END IF
         IF ( IFLAG.LT.0 )
     $      GO TO 40
C
C        If requested, call FCN to enable printing of iterates.
C
         IF ( NPRINT.GT.0 ) THEN
            IFLAG = 0
            IF ( MOD( ITER-1, NPRINT ).EQ.0 ) THEN
               CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                   LDPAR2, X, NFEV, DWORK(E), DWORK(JAC), LDJ,
     $                   DWORK(JTE), DWORK(JW1), LDWORK-JW1+1, INFOL )
C
               IF ( IFLAG.LT.0 )
     $            GO TO 40
            END IF
         END IF
C
C        Beginning of the inner loop.
C
   20    CONTINUE
C
C           Store the Levenberg factor in DWORK(E) (which is no longer
C           needed), to pass it to JPJ routine.
C
            DWORK(E) = PAR
C
C           Solve (J'*J + PAR*I)*x = J'*e, and store x in DWORK(IW1).
C           Additional workspace:
C                      N*N + DW(JPJ),          if ALG = 'D', STOR = 'F';
C                      N*( N + 1)/2 + DW(JPJ), if ALG = 'D', STOR = 'P';
C                      3*N + DW(JPJ),          if ALG = 'I'.
C
            IF ( CHOL ) THEN
               CALL DCOPY( N, DWORK(JTE), 1, DWORK(IW1), 1 )
               CALL MB02XD( 'Function', STOR, UPLO, JPJ, M, N, 1, IPAR,
     $                      LIPAR, DWORK(E), 1, DWORK(JAC), LDJ,
     $                      DWORK(IW1), N, DWORK(DWJTJ), N,
     $                      DWORK(JWORK), LDWORK-JWORK+1, INFOL )
            ELSE
               CALL MB02WD( 'Function', JPJ, N, IPAR, LIPAR, DWORK(E),
     $                      1, 3*N, DWORK(JAC), LDJ, DWORK(JTE), 1,
     $                      DWORK(IW1), 1, CGTOL*GNORM, DWORK(JWORK),
     $                      LDWORK-JWORK+1, IWARN, INFOL )
               ITERCG = ITERCG + INT( DWORK(JWORK) )
               IWARNL = MAX( 2*IWARN, IWARNL )
            ENDIF
C
            IF ( INFOL.NE.0 ) THEN
               INFO = 3
               RETURN
            ENDIF
C
C           Compute updated X.
C
            DO 30 I = 0, N - 1
               DWORK(IW2+I) = X(I+1) - DWORK(IW1+I)
   30       CONTINUE
C
C           Evaluate the function at x - p and calculate its norm.
C           Workspace:  need:    SIZEJ + M + 3*N + LFCN1;
C                       prefer:  larger.
C
            IFLAG = 1
            CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                LDPAR2, DWORK(IW2), NFEVL, DWORK(E), DWORK(JAC),
     $                LDJ, DWORK(JTE), DWORK(JW2), LDWORK-JW2+1, INFOL )
C
            IF ( INFOL.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
C
            NFEV = NFEV + 1
            IF ( IFLAG.LT.0 )
     $         GO TO 40
            FNORM1 = DNRM2( M, DWORK(E), 1 )
C
C           Now, check whether this step was successful and update the
C           Levenberg factor.
C
            IF ( FNORM.LT.FNORM1 ) THEN
C
C              Unsuccessful step: increase PAR.
C
               ACTRED = ONE
               IF ( PAR.GT.PARMAX ) THEN
                  IF ( PAR/MARQF.LE.BIGNUM )
     $               PAR = PAR*MARQF
               ELSE
                  PAR = PAR*MARQF
               END IF
C
            ELSE
C
C              Successful step: update PAR, X, and FNORM.
C
               ACTRED = ONE - ( FNORM1/FNORM )**2
               IF ( ( FNORM - FNORM1 )*( FNORM + FNORM1 ) .LT.
     $              MINIMP*DDOT( N, DWORK(IW1), 1,
     $                           DWORK(JTE), 1 ) ) THEN
                  IF ( PAR.GT.PARMAX ) THEN
                     IF ( PAR/MARQF.LE.BIGNUM )
     $                  PAR = PAR*MARQF
                  ELSE
                     PAR = PAR*MARQF
                  END IF
               ELSE
                  PAR = MAX( PAR/MARQF, SMLNUM )
               ENDIF
               CALL DCOPY( N, DWORK(IW2), 1, X, 1 )
               FNORM = FNORM1
            ENDIF
C
            IF ( ( ACTRED.LE.TOLDEF ) .OR. ( ITER.GT.ITMAX ) .OR.
     $           ( PAR.GT.PARMAX ) )
     $         GO TO 40
            IF ( ACTRED.LE.EPSMCH ) THEN
               IWARN = 4
               GO TO 40
            ENDIF
C
C           End of the inner loop. Repeat if unsuccessful iteration.
C
            IF ( FNORM.LT.FNORM1 )
     $         GO TO 20
C
C        End of the outer loop.
C
         GO TO 10
C
C     END WHILE 10
C
   40 CONTINUE
C
C     Termination, either normal or user imposed.
C
      IF ( ACTRED.GT.TOLDEF )
     $   IWARN = 1
      IF ( IWARNL.NE.0 )
     $   IWARN = 2
C
      IF ( IFLAG.LT.0 )
     $   IWARN = IFLAG
      IF ( NPRINT.GT.0 ) THEN
         IFLAG = 0
         CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $             LDPAR2, X, NFEV, DWORK(E), DWORK(JAC), LDJ,
     $             DWORK(JTE), DWORK(JW1), LDWORK-JW1+1, INFOL )
         IF ( IFLAG.LT.0 )
     $      IWARN = IFLAG
      END IF
C
      DWORK(1) = WRKOPT
      DWORK(2) = FNORM
      DWORK(3) = ITER
      DWORK(4) = ITERCG
      DWORK(5) = PAR
C
      RETURN
C *** Last line of MD03AD ***
      END
