      SUBROUTINE MD03BD( XINIT, SCALE, COND, FCN, QRFACT, LMPARM, M, N,
     $                   ITMAX, FACTOR, NPRINT, IPAR, LIPAR, DPAR1,
     $                   LDPAR1, DPAR2, LDPAR2, X, DIAG, NFEV, NJEV,
     $                   FTOL, XTOL, GTOL, TOL, IWORK, DWORK, LDWORK,
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
C     To minimize the sum of the squares of m nonlinear functions, e, in
C     n variables, x, by a modification of the Levenberg-Marquardt
C     algorithm. The user must provide a subroutine FCN which calculates
C     the functions and the Jacobian (possibly by finite differences).
C     In addition, specialized subroutines QRFACT, for QR factorization
C     with pivoting of the Jacobian, and LMPARM, for the computation of
C     Levenberg-Marquardt parameter, exploiting the possible structure
C     of the Jacobian matrix, should be provided. Template
C     implementations of these routines are included in SLICOT Library.
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
C     SCALE   CHARACTER*1
C             Specifies how the variables will be scaled, as follows:
C             = 'I' :  use internal scaling;
C             = 'S' :  use specified scaling factors, given in DIAG.
C
C     COND    CHARACTER*1
C             Specifies whether the condition of the linear systems
C             involved should be estimated, as follows:
C             = 'E' :  use incremental condition estimation to find the
C                      numerical rank;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of matrices for zero values.
C
C     Function Parameters
C
C     FCN     EXTERNAL
C             Subroutine which evaluates the functions and the Jacobian.
C             FCN must be declared in an external statement in the user
C             calling program, and must have the following interface:
C
C             SUBROUTINE FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1,
C            $                DPAR2, LDPAR2, X, NFEVL, E, J, LDJ, DWORK,
C            $                LDWORK, INFO )
C
C             where
C
C             IFLAG   (input/output) INTEGER
C                     On entry, this parameter must contain a value
C                     defining the computations to be performed:
C                     = 0 :  Optionally, print the current iterate X,
C                            function values E, and Jacobian matrix J,
C                            or other results defined in terms of these
C                            values. See the argument NPRINT of MD03BD.
C                            Do not alter E and J.
C                     = 1 :  Calculate the functions at X and return
C                            this vector in E. Do not alter J.
C                     = 2 :  Calculate the Jacobian at X and return
C                            this matrix in J. Also return NFEVL
C                            (see below). Do not alter E.
C                     = 3 :  Do not compute neither the functions nor
C                            the Jacobian, but return in LDJ and
C                            IPAR/DPAR1,DPAR2 (some of) the integer/real
C                            parameters needed.
C                     On exit, the value of this parameter should not be
C                     changed by FCN unless the user wants to terminate
C                     execution of MD03BD, in which case IFLAG must be
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
C                     QRFACT, and LMPARM, respectively.
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
C                     LDJ is essentially used inside the routines FCN,
C                     QRFACT and LMPARM.
C                     LDJ is an input parameter, except for IFLAG = 3
C                     on entry, when it is an output parameter.
C                     It is assumed in MD03BD that LDJ is not larger
C                     than needed.
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
C     QRFACT  EXTERNAL
C             Subroutine which computes the QR factorization with
C             (block) column pivoting of the Jacobian matrix, J*P = Q*R.
C             QRFACT must be declared in an external statement in the
C             calling program, and must have the following interface:
C
C             SUBROUTINE QRFACT( N, IPAR, LIPAR, FNORM, J, LDJ, E,
C            $                   JNORMS, GNORM, IPVT, DWORK, LDWORK,
C            $                   INFO )
C
C             where
C
C             N       (input) INTEGER
C                     The number of columns of the Jacobian matrix J.
C                     N >= 0.
C
C             IPAR    (input) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the Jacobian matrix.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 0.
C
C             FNORM   (input) DOUBLE PRECISION
C                     The Euclidean norm of the vector e.  FNORM >= 0.
C
C             J       (input/output) DOUBLE PRECISION array, dimension
C                     (LDJ, NC), where NC is the number of columns.
C                     On entry, the leading NR-by-NC part of this array
C                     must contain the (compressed) representation
C                     of the Jacobian matrix J, where NR is the number
C                     of rows of J (function of IPAR entries).
C                     On exit, the leading N-by-NC part of this array
C                     contains a (compressed) representation of the
C                     upper triangular factor R of the Jacobian matrix.
C                     For efficiency of the later calculations, the
C                     matrix R is delivered with the leading dimension
C                     MAX(1,N), possibly much smaller than the value
C                     of LDJ on entry.
C
C             LDJ     (input/output) INTEGER
C                     The leading dimension of array J.
C                     On entry, LDJ >= MAX(1,NR).
C                     On exit,  LDJ >= MAX(1,N).
C
C             E       (input/output) DOUBLE PRECISION array, dimension
C                     (NR)
C                     On entry, this array contains the error vector e.
C                     On exit, this array contains the updated vector
C                     Z*Q'*e, where Z is a block row permutation matrix
C                     (possibly identity) used in the QR factorization
C                     of J. (See, for example, the SLICOT Library
C                     routine NF01BS, Section METHOD.)
C
C             JNORMS  (output) DOUBLE PRECISION array, dimension (N)
C                     This array contains the Euclidean norms of the
C                     columns of the Jacobian matrix (in the original
C                     order).
C
C             GNORM   (output) DOUBLE PRECISION
C                     If FNORM > 0, the 1-norm of the scaled vector
C                     J'*e/FNORM, with each element i further divided
C                     by JNORMS(i) (if JNORMS(i) is nonzero).
C                     If FNORM = 0, the returned value of GNORM is 0.
C
C             IPVT    (output) INTEGER array, dimension (N)
C                     This array defines the permutation matrix P such
C                     that J*P = Q*R. Column j of P is column IPVT(j) of
C                     the identity matrix.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine QRFACT.
C                     On exit, if INFO = 0, DWORK(1) returns the optimal
C                     value of LDWORK.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine QRFACT).  LDWORK >= 1.
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input (scalar) argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine QRFACT. The LAPACK Library routine
C                     XERBLA should be used in conjunction with negative
C                     INFO. INFO must be zero if the subroutine finished
C                     successfully.
C
C             Parameters marked with "(input)" must not be changed.
C
C     LMPARM  EXTERNAL
C             Subroutine which determines a value for the Levenberg-
C             Marquardt parameter PAR such that if x solves the system
C
C                   J*x = b ,     sqrt(PAR)*D*x = 0 ,
C
C             in the least squares sense, where J is an m-by-n matrix,
C             D is an n-by-n nonsingular diagonal matrix, and b is an
C             m-vector, and if DELTA is a positive number, DXNORM is
C             the Euclidean norm of D*x, then either PAR is zero and
C
C                   ( DXNORM - DELTA ) .LE. 0.1*DELTA ,
C
C             or PAR is positive and
C
C                   ABS( DXNORM - DELTA ) .LE. 0.1*DELTA .
C
C             It is assumed that a block QR factorization, with column
C             pivoting, of J is available, that is, J*P = Q*R, where P
C             is a permutation matrix, Q has orthogonal columns, and
C             R is an upper triangular matrix (possibly stored in a
C             compressed form), with diagonal elements of nonincreasing
C             magnitude for each block. On output, LMPARM also provides
C             a (compressed) representation of an upper triangular
C             matrix S, such that
C
C                   P'*(J'*J + PAR*D*D)*P = S'*S .
C
C             LMPARM must be declared in an external statement in the
C             calling program, and must have the following interface:
C
C             SUBROUTINE LMPARM( COND, N, IPAR, LIPAR, R, LDR, IPVT,
C            $                   DIAG, QTB, DELTA, PAR, RANKS, X, RX,
C            $                   TOL, DWORK, LDWORK, INFO )
C
C             where
C
C             COND    CHARACTER*1
C                     Specifies whether the condition of the linear
C                     systems involved should be estimated, as follows:
C                     = 'E' :  use incremental condition estimation
C                              to find the numerical rank;
C                     = 'N' :  do not use condition estimation, but
C                              check the diagonal entries for zero
C                              values;
C                     = 'U' :  use the ranks already stored in RANKS
C                              (for R).
C
C             N       (input) INTEGER
C                     The order of the matrix R.  N >= 0.
C
C             IPAR    (input) INTEGER array, dimension (LIPAR)
C                     The integer parameters describing the structure of
C                     the Jacobian matrix.
C
C             LIPAR   (input) INTEGER
C                     The length of the array IPAR.  LIPAR >= 0.
C
C             R       (input/output) DOUBLE PRECISION array, dimension
C                     (LDR, NC), where NC is the number of columns.
C                     On entry, the leading N-by-NC part of this array
C                     must contain the (compressed) representation (Rc)
C                     of the upper triangular matrix R.
C                     On exit, the full upper triangular part of R
C                     (in representation Rc), is unaltered, and the
C                     remaining part contains (part of) the (compressed)
C                     representation of the transpose of the upper
C                     triangular matrix S.
C
C             LDR     (input) INTEGER
C                     The leading dimension of array R.
C                     LDR >= MAX(1,N).
C
C             IPVT    (input) INTEGER array, dimension (N)
C                     This array must define the permutation matrix P
C                     such that J*P = Q*R. Column j of P is column
C                     IPVT(j) of the identity matrix.
C
C             DIAG    (input) DOUBLE PRECISION array, dimension (N)
C                     This array must contain the diagonal elements of
C                     the matrix D.  DIAG(I) <> 0, I = 1,...,N.
C
C             QTB     (input) DOUBLE PRECISION array, dimension (N)
C                     This array must contain the first n elements of
C                     the vector Q'*b.
C
C             DELTA   (input) DOUBLE PRECISION
C                     An upper bound on the Euclidean norm of D*x.
C                     DELTA > 0.
C
C             PAR     (input/output) DOUBLE PRECISION
C                     On entry, PAR must contain an initial estimate of
C                     the Levenberg-Marquardt parameter.  PAR >= 0.
C                     On exit, it contains the final estimate of this
C                     parameter.
C
C             RANKS   (input or output) INTEGER array, dimension (r),
C                     where r is the number of diagonal blocks R_k in R,
C                     corresponding to the block column structure of J.
C                     On entry, if COND = 'U' and N > 0, this array must
C                     contain the numerical ranks of the submatrices
C                     R_k, k = 1:r. The number r is defined in terms of
C                     the entries of IPAR.
C                     On exit, if N > 0, this array contains the
C                     numerical ranks of the submatrices S_k, k = 1:r.
C
C             X       (output) DOUBLE PRECISION array, dimension (N)
C                     This array contains the least squares solution of
C                     the system J*x = b, sqrt(PAR)*D*x = 0.
C
C             RX      (output) DOUBLE PRECISION array, dimension (N)
C                     This array contains the matrix-vector product
C                     -R*P'*x.
C
C             TOL     (input) DOUBLE PRECISION
C                     If COND = 'E', the tolerance to be used for
C                     finding the ranks of the submatrices R_k and S_k.
C                     If the user sets TOL > 0, then the given value of
C                     TOL is used as a lower bound for the reciprocal
C                     condition number;  a (sub)matrix whose estimated
C                     condition number is less than 1/TOL is considered
C                     to be of full rank.  If the user sets TOL <= 0,
C                     then an implicitly computed, default tolerance,
C                     defined by TOLDEF = N*EPS,  is used instead,
C                     where EPS is the machine precision (see LAPACK
C                     Library routine DLAMCH).
C                     This parameter is not relevant if COND = 'U'
C                     or 'N'.
C
C             DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C                     The workspace array for subroutine LMPARM.
C                     On exit, if INFO = 0, DWORK(1) returns the optimal
C                     value of LDWORK.
C
C             LDWORK  (input) INTEGER
C                     The size of the array DWORK (as large as needed
C                     in the subroutine LMPARM).  LDWORK >= 1.
C
C             INFO    INTEGER
C                     Error indicator, set to a negative value if an
C                     input (scalar) argument is erroneous, and to
C                     positive values for other possible errors in the
C                     subroutine LMPARM. The LAPACK Library routine
C                     XERBLA should be used in conjunction with negative
C                     INFO. INFO must be zero if the subroutine finished
C                     successfully.
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
C     FACTOR  (input) DOUBLE PRECISION
C             The value used in determining the initial step bound. This
C             bound is set to the product of FACTOR and the Euclidean
C             norm of DIAG*X if nonzero, or else to FACTOR itself.
C             In most cases FACTOR should lie in the interval (.1,100).
C             A generally recommended value is 100.  FACTOR > 0.
C
C     NPRINT  (input) INTEGER
C             This parameter enables controlled printing of iterates if
C             it is positive. In this case, FCN is called with IFLAG = 0
C             at the beginning of the first iteration and every NPRINT
C             iterations thereafter and immediately prior to return,
C             with X, E, and J available for printing. Note that when
C             called immediately prior to return, J normally contains
C             the result returned by QRFACT and LMPARM (the compressed
C             R and S factors). If NPRINT is not positive, no special
C             calls of FCN with IFLAG = 0 are made.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters needed, for instance, for
C             describing the structure of the Jacobian matrix, which
C             are handed over to the routines FCN, QRFACT and LMPARM.
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
C             solving the problem. This argument is not used by MD03BD
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
C             solving the problem. This argument is not used by MD03BD
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
C             DWORK(1:4) give details on the iterative process.
C
C     DIAG    (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, if SCALE = 'S', this array must contain some
C             positive entries that serve as multiplicative scale
C             factors for the variables x.  DIAG(I) > 0, I = 1,...,N.
C             If SCALE = 'I', DIAG is internally set.
C             On exit, this array contains the scale factors used
C             (or finally used, if SCALE = 'I').
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
C     FTOL    DOUBLE PRECISION
C             If FTOL >= 0, the tolerance which measures the relative
C             error desired in the sum of squares. Termination occurs
C             when both the actual and predicted relative reductions in
C             the sum of squares are at most FTOL. If the user sets
C             FTOL < 0,  then  SQRT(EPS)  is used instead FTOL, where
C             EPS is the machine precision (see LAPACK Library routine
C             DLAMCH).
C
C     XTOL    DOUBLE PRECISION
C             If XTOL >= 0, the tolerance which measures the relative
C             error desired in the approximate solution. Termination
C             occurs when the relative error between two consecutive
C             iterates is at most XTOL. If the user sets  XTOL < 0,
C             then  SQRT(EPS)  is used instead XTOL.
C
C     GTOL    DOUBLE PRECISION
C             If GTOL >= 0, the tolerance which measures the
C             orthogonality desired between the function vector e and
C             the columns of the Jacobian J. Termination occurs when
C             the cosine of the angle between e and any column of the
C             Jacobian J is at most GTOL in absolute value. If the user
C             sets  GTOL < 0,  then  EPS  is used instead GTOL.
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             ranks of the matrices of linear systems to be solved. If
C             the user sets TOL > 0, then the given value of TOL is used
C             as a lower bound for the reciprocal condition number;  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*EPS,  is used instead.
C             This parameter is not relevant if COND = 'N'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N+r), where r is the number
C             of diagonal blocks R_k in R (see description of LMPARM).
C             On output, if INFO = 0, the first N entries of this array
C             define a permutation matrix P such that J*P = Q*R, where
C             J is the final calculated Jacobian, Q is an orthogonal
C             matrix (not stored), and R is upper triangular with
C             diagonal elements of nonincreasing magnitude (possibly
C             for each block column of J). Column j of P is column
C             IWORK(j) of the identity matrix. If INFO = 0, the entries
C             N+1:N+r of this array contain the ranks of the final
C             submatrices S_k (see description of LMPARM).
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, DWORK(2) returns the residual error norm (the
C             sum of squares), DWORK(3) returns the number of iterations
C             performed, and DWORK(4) returns the final Levenberg
C             factor. If INFO = 0, N > 0, and IWARN >= 0, the elements
C             DWORK(5) to DWORK(4+M) contain the final matrix-vector
C             product Z*Q'*e, and the elements DWORK(5+M) to
C             DWORK(4+M+N*NC) contain the (compressed) representation of
C             final upper triangular matrices R and S (if IWARN <> 4).
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= max( 4, M + max( size(J) +
C                                        max( DW( FCN|IFLAG = 1 ),
C                                             DW( FCN|IFLAG = 2 ),
C                                             DW( QRFACT ) + N ),
C                                        N*NC + N +
C                                        max( M + DW( FCN|IFLAG = 1 ),
C                                             N + DW( LMPARM ) ) ) ),
C             where size(J) is the size of the Jacobian (provided by FCN
C             in IPAR(1), for IFLAG = 3), and DW( f ) is the workspace
C             needed by the routine f, where f is FCN, QRFACT, or LMPARM
C             (provided by FCN in IPAR(2:5), for IFLAG = 3).
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             < 0:  the user set IFLAG = IWARN in the subroutine FCN;
C             = 1:  both actual and predicted relative reductions in
C                   the sum of squares are at most FTOL;
C             = 2:  relative error between two consecutive iterates is
C                   at most XTOL;
C             = 3:  conditions for IWARN = 1 and IWARN = 2 both hold;
C             = 4:  the cosine of the angle between e and any column of
C                   the Jacobian is at most GTOL in absolute value;
C             = 5:  the number of iterations has reached ITMAX without
C                   satisfying any convergence condition;
C             = 6:  FTOL is too small: no further reduction in the sum
C                   of squares is possible;
C             = 7:  XTOL is too small: no further improvement in the
C                   approximate solution x is possible;
C             = 8:  GTOL is too small: e is orthogonal to the columns of
C                   the Jacobian to machine precision.
C             In all these cases, DWORK(1:4) are set as described above.
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
C             = 3:  user-defined routine QRFACT returned with INFO <> 0;
C             = 4:  user-defined routine LMPARM returned with INFO <> 0.
C
C     METHOD
C
C     If XINIT = 'R', the initial value for x is set to a vector of
C     pseudo-random values uniformly distributed in (-1,1).
C
C     The Levenberg-Marquardt algorithm (described in [1,3]) is used for
C     optimizing the variables x. This algorithm needs the Jacobian
C     matrix J, which is provided by the subroutine FCN. A trust region
C     method is used. The algorithm tries to update x by the formula
C
C         x = x - p,
C
C     using an approximate solution of the system of linear equations
C
C         (J'*J + PAR*D*D)*p = J'*e,
C
C     with e the error function vector, and D a diagonal nonsingular
C     matrix, where either PAR = 0 and
C
C         ( norm( D*x ) - DELTA ) <= 0.1*DELTA ,
C
C     or PAR > 0 and
C
C         ABS( norm( D*x ) - DELTA ) <= 0.1*DELTA .
C
C     DELTA is the radius of the trust region. If the Gauss-Newton
C     direction is not acceptable, then an iterative algorithm obtains
C     improved lower and upper bounds for the Levenberg-Marquardt
C     parameter PAR. Only a few iterations are generally needed for
C     convergence of the algorithm. The trust region radius DELTA
C     and the Levenberg factor PAR are updated based on the ratio
C     between the actual and predicted reduction in the sum of squares.
C
C     REFERENCES
C
C     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E.
C         User's Guide for MINPACK-1.
C         Applied Math. Division, Argonne National Laboratory, Argonne,
C         Illinois, Report ANL-80-74, 1980.
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
C     The convergence rate near a local minimum is quadratic, if the
C     Jacobian is computed analytically, and linear, if the Jacobian
C     is computed numerically.
C
C     FURTHER COMMENTS
C
C     This routine is a more general version of the subroutines LMDER
C     and LMDER1 from the MINPACK package [1], which enables to exploit
C     the structure of the problem, and optionally use condition
C     estimation. Unstructured problems could be solved as well.
C
C     Template SLICOT Library implementations for FCN, QRFACT and
C     LMPARM routines are:
C     MD03BF, MD03BA, and MD03BB, respectively, for standard problems;
C     NF01BF, NF01BS, and NF01BP, respectively, for optimizing the
C     parameters of Wiener systems (structured problems).
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, Feb. 15, 2004.
C
C     KEYWORDS
C
C     Least-squares approximation,  Levenberg-Marquardt algorithm,
C     matrix operations, optimization.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, FOUR, P1, P5, P25, P75, P0001
      PARAMETER         ( ZERO = 0.0D0,  ONE   = 1.0D0,  FOUR = 4.0D0,
     $                    P1   = 1.0D-1, P5    = 5.0D-1, P25 = 2.5D-1,
     $                    P75  = 7.5D-1, P0001 = 1.0D-4 )
C     .. Scalar Arguments ..
      CHARACTER         COND, SCALE, XINIT
      INTEGER           INFO, ITMAX, IWARN, LDPAR1, LDPAR2, LDWORK,
     $                  LIPAR, M, N, NFEV, NJEV, NPRINT
      DOUBLE PRECISION  FACTOR, FTOL, GTOL, TOL, XTOL
C     .. Array Arguments ..
      INTEGER           IPAR(*), IWORK(*)
      DOUBLE PRECISION  DIAG(*), DPAR1(*), DPAR2(*), DWORK(*), X(*)
C     .. Local Scalars ..
      LOGICAL           BADSCL, INIT, ISCAL, SSCAL
      INTEGER           E, IFLAG, INFOL, ITER, IW1, IW2, IW3, J, JAC,
     $                  JW1, JW2, JWORK, L, LDJ, LDJSAV, LFCN1, LFCN2,
     $                  LLMP, LQRF, NC, NFEVL, SIZEJ, WRKOPT
      DOUBLE PRECISION  ACTRED, DELTA, DIRDER, EPSMCH, FNORM, FNORM1,
     $                  FTDEF, GNORM, GTDEF, PAR, PNORM, PRERED, RATIO,
     $                  TEMP, TEMP1, TEMP2, TOLDEF, XNORM, XTDEF
C     .. Local Arrays ..
      INTEGER           SEED(4)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLARNV, FCN, LMPARM, QRFACT, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN, MOD, SQRT
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INIT  = LSAME( XINIT, 'R' )
      ISCAL = LSAME( SCALE, 'I' )
      SSCAL = LSAME( SCALE, 'S' )
      INFO  = 0
      IWARN = 0
      IF( .NOT.( INIT .OR. LSAME( XINIT, 'G' ) ) ) THEN
         INFO = -1
      ELSEIF( .NOT.( ISCAL .OR. SSCAL ) ) THEN
         INFO = -2
      ELSEIF( .NOT.( LSAME( COND, 'E' ) .OR. LSAME( COND, 'N' ) ) ) THEN
         INFO = -3
      ELSEIF( M.LT.0 ) THEN
         INFO = -7
      ELSEIF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -8
      ELSEIF( ITMAX.LT.0 ) THEN
         INFO = -9
      ELSEIF( FACTOR.LE.ZERO ) THEN
         INFO = -10
      ELSEIF( LIPAR.LT.5 ) THEN
         INFO = -13
      ELSEIF( LDPAR1.LT.0 ) THEN
         INFO = -15
      ELSEIF( LDPAR2.LT.0 ) THEN
         INFO = -17
      ELSEIF ( LDWORK.LT.4 ) THEN
         INFO = -28
      ELSEIF ( SSCAL ) THEN
         BADSCL = .FALSE.
C
         DO 10 J = 1, N
            BADSCL = BADSCL .OR. DIAG(J).LE.ZERO
   10    CONTINUE
C
         IF ( BADSCL )
     $      INFO = -19
      END IF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MD03BD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      NFEV = 0
      NJEV = 0
      IF ( N.EQ.0 ) THEN
         DWORK(1) = FOUR
         DWORK(2) = ZERO
         DWORK(3) = ZERO
         DWORK(4) = ZERO
         RETURN
      END IF
C
C     Call FCN to get the size of the array J, for storing the Jacobian
C     matrix, the leading dimension LDJ and the workspace required
C     by FCN for IFLAG = 1 and IFLAG = 2, QRFACT and LMPARM. The
C     entries DWORK(1:4) should not be modified by the special call of
C     FCN below, if XINIT = 'R' and the values in DWORK(1:4) are
C     explicitly desired for initialization of the random number
C     generator.
C
      IFLAG = 3
      IW1   = IPAR(1)
      IW2   = IPAR(2)
      IW3   = IPAR(3)
      JW1   = IPAR(4)
      JW2   = IPAR(5)
C
      CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2, LDPAR2,
     $          X, NFEVL, DWORK, DWORK, LDJSAV, DWORK, LDWORK, INFOL )
      SIZEJ = IPAR(1)
      LFCN1 = IPAR(2)
      LFCN2 = IPAR(3)
      LQRF  = IPAR(4)
      LLMP  = IPAR(5)
      IF ( LDJSAV.GT.0 ) THEN
         NC = SIZEJ/LDJSAV
      ELSE
         NC = SIZEJ
      END IF
C
      IPAR(1) = IW1
      IPAR(2) = IW2
      IPAR(3) = IW3
      IPAR(4) = JW1
      IPAR(5) = JW2
C
C     Check the workspace length.
C
      E     = 1
      JAC   = E   + M
      JW1   = JAC + SIZEJ
      JW2   = JW1 + N
      IW1   = JAC + N*NC
      IW2   = IW1 + N
      IW3   = IW2 + N
      JWORK = IW2 + M
C
      L = MAX( 4, M + MAX( SIZEJ + MAX( LFCN1, LFCN2, N + LQRF ),
     $                     N*NC + N + MAX( M + LFCN1, N + LLMP ) ) )
      IF ( LDWORK.LT.L ) THEN
         INFO = -28
         CALL XERBLA( 'MD03BD', -INFO )
         RETURN
      ENDIF
C
C     Set default tolerances. EPSMCH is the machine precision.
C
      EPSMCH = DLAMCH( 'Epsilon' )
      FTDEF  = FTOL
      XTDEF  = XTOL
      GTDEF  = GTOL
      TOLDEF =  TOL
      IF ( MIN( FTDEF, XTDEF, GTDEF, TOLDEF ).LE.ZERO ) THEN
         IF ( FTDEF.LT.ZERO )
     $      FTDEF = SQRT( EPSMCH )
         IF ( XTDEF.LT.ZERO )
     $      XTDEF = SQRT( EPSMCH )
         IF ( GTDEF.LT.ZERO )
     $      GTDEF = EPSMCH
         IF ( TOLDEF.LE.ZERO )
     $      TOLDEF = DBLE( N )*EPSMCH
      ENDIF
      WRKOPT = 1
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
C     Initialize Levenberg-Marquardt parameter and iteration counter.
C
      PAR  = ZERO
      ITER = 1
C
C     Evaluate the function at the starting point
C     and calculate its norm.
C     Workspace: need:    M + SIZEJ + LFCN1;
C                prefer:  larger.
C
      IFLAG = 1
      CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2, LDPAR2,
     $          X, NFEVL, DWORK(E), DWORK(JAC), LDJ, DWORK(JW1),
     $          LDWORK-JW1+1, INFOL )
C
      IF ( INFOL.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(JW1) ) + JW1 - 1 )
      NFEV   = 1
      FNORM  = DNRM2( M, DWORK(E), 1 )
      IF ( IFLAG.LT.0 .OR. FNORM.EQ.ZERO )
     $   GO TO 90
C
C     Beginning of the outer loop.
C
   20 CONTINUE
C
C        Calculate the Jacobian matrix.
C        Workspace: need:    M + SIZEJ + LFCN2;
C                   prefer:  larger.
C
         LDJ   = LDJSAV
         IFLAG = 2
         CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $             LDPAR2, X, NFEVL, DWORK(E), DWORK(JAC), LDJ,
     $             DWORK(JW1), LDWORK-JW1+1, INFOL )
C
         IF ( INFOL.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
         IF ( ITER.EQ.1 )
     $      WRKOPT = MAX( WRKOPT, INT( DWORK(JW1) ) + JW1 - 1 )
         IF ( NFEVL.GT.0 )
     $      NFEV = NFEV + NFEVL
         NJEV = NJEV + 1
         IF ( IFLAG.LT.0 )
     $      GO TO 90
C
C        If requested, call FCN to enable printing of iterates.
C
         IF ( NPRINT.GT.0 ) THEN
            IFLAG = 0
            IF ( MOD( ITER-1, NPRINT ).EQ.0 ) THEN
               CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                   LDPAR2, X, NFEV, DWORK(E), DWORK(JAC), LDJ,
     $                   DWORK(JW1), LDWORK-JW1+1, INFOL )
C
               IF ( IFLAG.LT.0 )
     $            GO TO 90
            END IF
         END IF
C
C        Compute the QR factorization of the Jacobian.
C        Workspace: need:    M + SIZEJ + N + LQRF;
C                   prefer:  larger.
C
         CALL QRFACT( N, IPAR, LIPAR, FNORM, DWORK(JAC), LDJ, DWORK(E),
     $                DWORK(JW1), GNORM, IWORK, DWORK(JW2),
     $                LDWORK-JW2+1, INFOL )
         IF ( INFOL.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
C
C        On the first iteration and if SCALE = 'I', scale according
C        to the norms of the columns of the initial Jacobian.
C
         IF ( ITER.EQ.1 ) THEN
            WRKOPT = MAX( WRKOPT, INT( DWORK(JW2) ) + JW2 - 1 )
            IF ( ISCAL ) THEN
C
               DO 30 J = 1, N
                  DIAG(J) = DWORK(JW1+J-1)
                  IF ( DIAG(J).EQ.ZERO )
     $               DIAG(J) = ONE
   30          CONTINUE
C
            END IF
C
C           On the first iteration, calculate the norm of the scaled
C           x and initialize the step bound DELTA.
C
            DO 40 J = 1, N
               DWORK(IW1+J-1) = DIAG(J)*X(J)
   40       CONTINUE
C
            XNORM = DNRM2( N, DWORK(IW1), 1 )
            DELTA = FACTOR*XNORM
            IF ( DELTA.EQ.ZERO )
     $         DELTA = FACTOR
         ELSE
C
C           Rescale if necessary.
C
            IF ( ISCAL ) THEN
C
               DO 50 J = 1, N
                  DIAG(J) = MAX( DIAG(J), DWORK(JW1+J-1) )
   50          CONTINUE
C
            END IF
         END IF
C
C        Test for convergence of the gradient norm.
C
         IF ( GNORM.LE.GTDEF )
     $      IWARN = 4
         IF ( IWARN.NE.0 )
     $      GO TO 90
C
C        Beginning of the inner loop.
C
   60    CONTINUE
C
C           Determine the Levenberg-Marquardt parameter and the
C           direction p, and compute -R*P'*p.
C           Workspace:  need:    M + N*NC + 2*N + LLMP;
C                       prefer:  larger.
C
            CALL LMPARM( COND, N, IPAR, LIPAR, DWORK(JAC), LDJ,
     $                   IWORK, DIAG, DWORK(E), DELTA, PAR, IWORK(N+1),
     $                   DWORK(IW1), DWORK(IW2), TOLDEF, DWORK(IW3),
     $                   LDWORK-IW3+1, INFOL )
            IF ( INFOL.NE.0 ) THEN
               INFO = 4
               RETURN
            END IF
            IF ( ITER.EQ.1 )
     $         WRKOPT = MAX( WRKOPT, INT( DWORK(IW3) ) + IW3 - 1 )
C
            TEMP1 = DNRM2( N, DWORK(IW2), 1 )/FNORM
C
C           Store the direction p and x - p.
C
            DO 70 J = 0, N - 1
               DWORK(IW2+J) = DIAG(J+1)*DWORK(IW1+J)
               DWORK(IW1+J) = X(J+1)  - DWORK(IW1+J)
   70       CONTINUE
C
C           Compute the norm of scaled p and the scaled predicted
C           reduction and the scaled directional derivative.
C
            PNORM = DNRM2( N, DWORK(IW2), 1 )
            TEMP2 = ( SQRT( PAR )*PNORM )/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -( TEMP1**2 + TEMP2**2 )
C
C           On the first iteration, adjust the initial step bound.
C
            IF ( ITER.EQ.1 )
     $         DELTA = MIN( DELTA, PNORM )
C
C           Evaluate the function at x - p and calculate its norm.
C           Workspace:  need:    2*M + N*NC + N + LFCN1;
C                       prefer:  larger.
C
            IFLAG = 1
            CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                LDPAR2, DWORK(IW1), NFEVL, DWORK(IW2), DWORK(JAC),
     $                LDJ, DWORK(JWORK), LDWORK-JWORK+1, INFOL )
            IF ( INFOL.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
C
            NFEV = NFEV + 1
            IF ( IFLAG.LT.0 )
     $         GO TO 90
            FNORM1 = DNRM2( M, DWORK(IW2), 1 )
C
C           Compute the scaled actual reduction.
C
            ACTRED = -ONE
            IF ( P1*FNORM1.LT.FNORM )
     $         ACTRED = ONE - ( FNORM1/FNORM )**2
C
C           Compute the ratio of the actual to the predicted reduction.
C
            RATIO = ZERO
            IF ( PRERED.NE.ZERO )
     $         RATIO = ACTRED/PRERED
C
C           Update the step bound.
C
            IF ( RATIO.LE.P25 ) THEN
               IF ( ACTRED.GE.ZERO ) THEN
                  TEMP = P5
               ELSE
                  TEMP = P5*DIRDER/( DIRDER + P5*ACTRED )
               END IF
               IF ( P1*FNORM1.GE.FNORM .OR. TEMP.LT.P1 )
     $            TEMP = P1
               DELTA = TEMP*MIN( DELTA, PNORM/P1 )
               PAR = PAR/TEMP
            ELSE
               IF ( PAR.EQ.ZERO .OR. RATIO.GE.P75 ) THEN
                  DELTA = PNORM/P5
                  PAR   = P5*PAR
               END IF
            END IF
C
C           Test for successful iteration.
C
            IF ( RATIO.GE.P0001 ) THEN
C
C              Successful iteration. Update x, e, and their norms.
C
               DO 80 J = 1, N
                  X(J) = DWORK(IW1+J-1)
                  DWORK(IW1+J-1) = DIAG(J)*X(J)
   80          CONTINUE
C
               CALL DCOPY( M, DWORK(IW2), 1, DWORK(E), 1 )
               XNORM = DNRM2( N, DWORK(IW1), 1 )
               FNORM = FNORM1
               ITER = ITER + 1
            END IF
C
C           Tests for convergence.
C
            IF ( ABS( ACTRED ).LE.FTDEF .AND. PRERED.LE.FTDEF .AND.
     $          P5*RATIO.LE.ONE )
     $         IWARN = 1
            IF ( DELTA.LE.XTDEF*XNORM )
     $         IWARN = 2
            IF ( ABS( ACTRED ).LE.FTDEF .AND. PRERED.LE.FTDEF .AND.
     $          P5*RATIO.LE.ONE .AND. IWARN.EQ.2 )
     $         IWARN = 3
            IF ( IWARN.NE.0 )
     $         GO TO 90
C
C           Tests for termination and stringent tolerances.
C
            IF ( ITER.GE.ITMAX )
     $         IWARN = 5
            IF ( ABS( ACTRED ).LE.EPSMCH .AND. PRERED.LE.EPSMCH .AND.
     $          P5*RATIO.LE.ONE )
     $         IWARN = 6
            IF ( DELTA.LE.EPSMCH*XNORM )
     $         IWARN = 7
            IF ( GNORM.LE.EPSMCH )
     $         IWARN = 8
            IF ( IWARN.NE.0 )
     $         GO TO 90
C
C           End of the inner loop. Repeat if unsuccessful iteration.
C
            IF ( RATIO.LT.P0001 ) GO TO 60
C
C        End of the outer loop.
C
         GO TO 20
C
   90 CONTINUE
C
C     Termination, either normal or user imposed.
C     Note that DWORK(JAC) normally contains the results returned by
C     QRFACT and LMPARM (the compressed R and S factors).
C
      IF ( IFLAG.LT.0 )
     $   IWARN = IFLAG
      IF ( NPRINT.GT.0 ) THEN
         IFLAG = 0
         CALL FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $             LDPAR2, X, NFEV, DWORK(E), DWORK(JAC), LDJ,
     $             DWORK(JWORK), LDWORK-JWORK+1, INFOL )
         IF ( IFLAG.LT.0 )
     $      IWARN = IFLAG
      END IF
C
      IF ( IWARN.GE.0 ) THEN
         DO 100 J = M + N*NC, 1, -1
            DWORK(4+J) = DWORK(J)
  100    CONTINUE
      END IF
      DWORK(1) = WRKOPT
      DWORK(2) = FNORM
      DWORK(3) = ITER
      DWORK(4) = PAR
C
      RETURN
C *** Last line of MD03BD ***
      END
