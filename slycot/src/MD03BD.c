/* MD03BD.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int md03bd_(char *xinit, char *scale, char *cond, S_fp fcn, 
	S_fp qrfact, S_fp lmparm, integer *m, integer *n, integer *itmax, 
	doublereal *factor, integer *nprint, integer *ipar, integer *lipar, 
	doublereal *dpar1, integer *ldpar1, doublereal *dpar2, integer *
	ldpar2, doublereal *x, doublereal *diag, integer *nfev, integer *njev,
	 doublereal *ftol, doublereal *xtol, doublereal *gtol, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen xinit_len, ftnlen scale_len, ftnlen 
	cond_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer e, j, l, nc, iw1, iw2, iw3, jw1, jw2, jac, ldj;
    static doublereal par;
    static integer seed[4];
    static logical init;
    static integer iter, llmp, lqrf;
    static doublereal temp;
    static integer lfcn1, lfcn2;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1, temp2;
    static integer iflag;
    static doublereal ftdef, delta, gtdef;
    static logical iscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical sscal;
    static integer infol, nfevl;
    static doublereal xtdef, ratio;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fnorm, gnorm;
    static integer sizej, jwork;
    static doublereal pnorm, xnorm, fnorm1;
    static logical badscl;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal actred, dirder, toldef, epsmch, prered;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ldjsav;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
    static integer wrkopt;


/*     SLICOT RELEASE 5.0. */

/*     Copyright (c) 2002-2009 NICONET e.V. */

/*     This program is free software: you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation, either version 2 of */
/*     the License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see */
/*     <http://www.gnu.org/licenses/>. */

/*     PURPOSE */

/*     To minimize the sum of the squares of m nonlinear functions, e, in */
/*     n variables, x, by a modification of the Levenberg-Marquardt */
/*     algorithm. The user must provide a subroutine FCN which calculates */
/*     the functions and the Jacobian (possibly by finite differences). */
/*     In addition, specialized subroutines QRFACT, for QR factorization */
/*     with pivoting of the Jacobian, and LMPARM, for the computation of */
/*     Levenberg-Marquardt parameter, exploiting the possible structure */
/*     of the Jacobian matrix, should be provided. Template */
/*     implementations of these routines are included in SLICOT Library. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     XINIT   CHARACTER*1 */
/*             Specifies how the variables x are initialized, as follows: */
/*             = 'R' :  the array X is initialized to random values; the */
/*                      entries DWORK(1:4) are used to initialize the */
/*                      random number generator: the first three values */
/*                      are converted to integers between 0 and 4095, and */
/*                      the last one is converted to an odd integer */
/*                      between 1 and 4095; */
/*             = 'G' :  the given entries of X are used as initial values */
/*                      of variables. */

/*     SCALE   CHARACTER*1 */
/*             Specifies how the variables will be scaled, as follows: */
/*             = 'I' :  use internal scaling; */
/*             = 'S' :  use specified scaling factors, given in DIAG. */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the linear systems */
/*             involved should be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation to find the */
/*                      numerical rank; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of matrices for zero values. */

/*     Function Parameters */

/*     FCN     EXTERNAL */
/*             Subroutine which evaluates the functions and the Jacobian. */
/*             FCN must be declared in an external statement in the user */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, */
/*            $                DPAR2, LDPAR2, X, NFEVL, E, J, LDJ, DWORK, */
/*            $                LDWORK, INFO ) */

/*             where */

/*             IFLAG   (input/output) INTEGER */
/*                     On entry, this parameter must contain a value */
/*                     defining the computations to be performed: */
/*                     = 0 :  Optionally, print the current iterate X, */
/*                            function values E, and Jacobian matrix J, */
/*                            or other results defined in terms of these */
/*                            values. See the argument NPRINT of MD03BD. */
/*                            Do not alter E and J. */
/*                     = 1 :  Calculate the functions at X and return */
/*                            this vector in E. Do not alter J. */
/*                     = 2 :  Calculate the Jacobian at X and return */
/*                            this matrix in J. Also return NFEVL */
/*                            (see below). Do not alter E. */
/*                     = 3 :  Do not compute neither the functions nor */
/*                            the Jacobian, but return in LDJ and */
/*                            IPAR/DPAR1,DPAR2 (some of) the integer/real */
/*                            parameters needed. */
/*                     On exit, the value of this parameter should not be */
/*                     changed by FCN unless the user wants to terminate */
/*                     execution of MD03BD, in which case IFLAG must be */
/*                     set to a negative integer. */

/*             M       (input) INTEGER */
/*                     The number of functions.  M >= 0. */

/*             N       (input) INTEGER */
/*                     The number of variables.  M >= N >= 0. */

/*             IPAR    (input/output) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix or needed for problem solving. */
/*                     IPAR is an input parameter, except for IFLAG = 3 */
/*                     on entry, when it is also an output parameter. */
/*                     On exit, if IFLAG = 3, IPAR(1) contains the length */
/*                     of the array J, for storing the Jacobian matrix, */
/*                     and the entries IPAR(2:5) contain the workspace */
/*                     required by FCN for IFLAG = 1, FCN for IFLAG = 2, */
/*                     QRFACT, and LMPARM, respectively. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 5. */

/*             DPAR1   (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDPAR1,*) or (LDPAR1) */
/*                     A first set of real parameters needed for */
/*                     describing or solving the problem. */
/*                     DPAR1 can also be used as an additional array for */
/*                     intermediate results when computing the functions */
/*                     or the Jacobian. For control problems, DPAR1 could */
/*                     store the input trajectory of a system. */

/*             LDPAR1  (input) INTEGER */
/*                     The leading dimension or the length of the array */
/*                     DPAR1, as convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, */
/*                     if leading dimension.) */

/*             DPAR2   (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDPAR2,*) or (LDPAR2) */
/*                     A second set of real parameters needed for */
/*                     describing or solving the problem. */
/*                     DPAR2 can also be used as an additional array for */
/*                     intermediate results when computing the functions */
/*                     or the Jacobian. For control problems, DPAR2 could */
/*                     store the output trajectory of a system. */

/*             LDPAR2  (input) INTEGER */
/*                     The leading dimension or the length of the array */
/*                     DPAR2, as convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1, */
/*                     if leading dimension.) */

/*             X       (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the value of the */
/*                     variables x where the functions or the Jacobian */
/*                     must be evaluated. */

/*             NFEVL   (input/output) INTEGER */
/*                     The number of function evaluations needed to */
/*                     compute the Jacobian by a finite difference */
/*                     approximation. */
/*                     NFEVL is an input parameter if IFLAG = 0, or an */
/*                     output parameter if IFLAG = 2. If the Jacobian is */
/*                     computed analytically, NFEVL should be set to a */
/*                     non-positive value. */

/*             E       (input/output) DOUBLE PRECISION array, */
/*                     dimension (M) */
/*                     This array contains the value of the (error) */
/*                     functions e evaluated at X. */
/*                     E is an input parameter if IFLAG = 0 or 2, or an */
/*                     output parameter if IFLAG = 1. */

/*             J       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDJ,NC), where NC is the number of columns */
/*                     needed. */
/*                     This array contains a possibly compressed */
/*                     representation of the Jacobian matrix evaluated */
/*                     at X. If full Jacobian is stored, then NC = N. */
/*                     J is an input parameter if IFLAG = 0, or an output */
/*                     parameter if IFLAG = 2. */

/*             LDJ     (input/output) INTEGER */
/*                     The leading dimension of array J.  LDJ >= 1. */
/*                     LDJ is essentially used inside the routines FCN, */
/*                     QRFACT and LMPARM. */
/*                     LDJ is an input parameter, except for IFLAG = 3 */
/*                     on entry, when it is an output parameter. */
/*                     It is assumed in MD03BD that LDJ is not larger */
/*                     than needed. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine FCN. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine FCN).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine FCN. The LAPACK Library routine XERBLA */
/*                     should be used in conjunction with negative INFO. */
/*                     INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     QRFACT  EXTERNAL */
/*             Subroutine which computes the QR factorization with */
/*             (block) column pivoting of the Jacobian matrix, J*P = Q*R. */
/*             QRFACT must be declared in an external statement in the */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE QRFACT( N, IPAR, LIPAR, FNORM, J, LDJ, E, */
/*            $                   JNORMS, GNORM, IPVT, DWORK, LDWORK, */
/*            $                   INFO ) */

/*             where */

/*             N       (input) INTEGER */
/*                     The number of columns of the Jacobian matrix J. */
/*                     N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             FNORM   (input) DOUBLE PRECISION */
/*                     The Euclidean norm of the vector e.  FNORM >= 0. */

/*             J       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDJ, NC), where NC is the number of columns. */
/*                     On entry, the leading NR-by-NC part of this array */
/*                     must contain the (compressed) representation */
/*                     of the Jacobian matrix J, where NR is the number */
/*                     of rows of J (function of IPAR entries). */
/*                     On exit, the leading N-by-NC part of this array */
/*                     contains a (compressed) representation of the */
/*                     upper triangular factor R of the Jacobian matrix. */
/*                     For efficiency of the later calculations, the */
/*                     matrix R is delivered with the leading dimension */
/*                     MAX(1,N), possibly much smaller than the value */
/*                     of LDJ on entry. */

/*             LDJ     (input/output) INTEGER */
/*                     The leading dimension of array J. */
/*                     On entry, LDJ >= MAX(1,NR). */
/*                     On exit,  LDJ >= MAX(1,N). */

/*             E       (input/output) DOUBLE PRECISION array, dimension */
/*                     (NR) */
/*                     On entry, this array contains the error vector e. */
/*                     On exit, this array contains the updated vector */
/*                     Z*Q'*e, where Z is a block row permutation matrix */
/*                     (possibly identity) used in the QR factorization */
/*                     of J. (See, for example, the SLICOT Library */
/*                     routine NF01BS, Section METHOD.) */

/*             JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the Euclidean norms of the */
/*                     columns of the Jacobian matrix (in the original */
/*                     order). */

/*             GNORM   (output) DOUBLE PRECISION */
/*                     If FNORM > 0, the 1-norm of the scaled vector */
/*                     J'*e/FNORM, with each element i further divided */
/*                     by JNORMS(i) (if JNORMS(i) is nonzero). */
/*                     If FNORM = 0, the returned value of GNORM is 0. */

/*             IPVT    (output) INTEGER array, dimension (N) */
/*                     This array defines the permutation matrix P such */
/*                     that J*P = Q*R. Column j of P is column IPVT(j) of */
/*                     the identity matrix. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine QRFACT. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine QRFACT).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine QRFACT. The LAPACK Library routine */
/*                     XERBLA should be used in conjunction with negative */
/*                     INFO. INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     LMPARM  EXTERNAL */
/*             Subroutine which determines a value for the Levenberg- */
/*             Marquardt parameter PAR such that if x solves the system */

/*                   J*x = b ,     sqrt(PAR)*D*x = 0 , */

/*             in the least squares sense, where J is an m-by-n matrix, */
/*             D is an n-by-n nonsingular diagonal matrix, and b is an */
/*             m-vector, and if DELTA is a positive number, DXNORM is */
/*             the Euclidean norm of D*x, then either PAR is zero and */

/*                   ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*             or PAR is positive and */

/*                   ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*             It is assumed that a block QR factorization, with column */
/*             pivoting, of J is available, that is, J*P = Q*R, where P */
/*             is a permutation matrix, Q has orthogonal columns, and */
/*             R is an upper triangular matrix (possibly stored in a */
/*             compressed form), with diagonal elements of nonincreasing */
/*             magnitude for each block. On output, LMPARM also provides */
/*             a (compressed) representation of an upper triangular */
/*             matrix S, such that */

/*                   P'*(J'*J + PAR*D*D)*P = S'*S . */

/*             LMPARM must be declared in an external statement in the */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE LMPARM( COND, N, IPAR, LIPAR, R, LDR, IPVT, */
/*            $                   DIAG, QTB, DELTA, PAR, RANKS, X, RX, */
/*            $                   TOL, DWORK, LDWORK, INFO ) */

/*             where */

/*             COND    CHARACTER*1 */
/*                     Specifies whether the condition of the linear */
/*                     systems involved should be estimated, as follows: */
/*                     = 'E' :  use incremental condition estimation */
/*                              to find the numerical rank; */
/*                     = 'N' :  do not use condition estimation, but */
/*                              check the diagonal entries for zero */
/*                              values; */
/*                     = 'U' :  use the ranks already stored in RANKS */
/*                              (for R). */

/*             N       (input) INTEGER */
/*                     The order of the matrix R.  N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             R       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDR, NC), where NC is the number of columns. */
/*                     On entry, the leading N-by-NC part of this array */
/*                     must contain the (compressed) representation (Rc) */
/*                     of the upper triangular matrix R. */
/*                     On exit, the full upper triangular part of R */
/*                     (in representation Rc), is unaltered, and the */
/*                     remaining part contains (part of) the (compressed) */
/*                     representation of the transpose of the upper */
/*                     triangular matrix S. */

/*             LDR     (input) INTEGER */
/*                     The leading dimension of array R. */
/*                     LDR >= MAX(1,N). */

/*             IPVT    (input) INTEGER array, dimension (N) */
/*                     This array must define the permutation matrix P */
/*                     such that J*P = Q*R. Column j of P is column */
/*                     IPVT(j) of the identity matrix. */

/*             DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the diagonal elements of */
/*                     the matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*             QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the first n elements of */
/*                     the vector Q'*b. */

/*             DELTA   (input) DOUBLE PRECISION */
/*                     An upper bound on the Euclidean norm of D*x. */
/*                     DELTA > 0. */

/*             PAR     (input/output) DOUBLE PRECISION */
/*                     On entry, PAR must contain an initial estimate of */
/*                     the Levenberg-Marquardt parameter.  PAR >= 0. */
/*                     On exit, it contains the final estimate of this */
/*                     parameter. */

/*             RANKS   (input or output) INTEGER array, dimension (r), */
/*                     where r is the number of diagonal blocks R_k in R, */
/*                     corresponding to the block column structure of J. */
/*                     On entry, if COND = 'U' and N > 0, this array must */
/*                     contain the numerical ranks of the submatrices */
/*                     R_k, k = 1:r. The number r is defined in terms of */
/*                     the entries of IPAR. */
/*                     On exit, if N > 0, this array contains the */
/*                     numerical ranks of the submatrices S_k, k = 1:r. */

/*             X       (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the least squares solution of */
/*                     the system J*x = b, sqrt(PAR)*D*x = 0. */

/*             RX      (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the matrix-vector product */
/*                     -R*P'*x. */

/*             TOL     (input) DOUBLE PRECISION */
/*                     If COND = 'E', the tolerance to be used for */
/*                     finding the ranks of the submatrices R_k and S_k. */
/*                     If the user sets TOL > 0, then the given value of */
/*                     TOL is used as a lower bound for the reciprocal */
/*                     condition number;  a (sub)matrix whose estimated */
/*                     condition number is less than 1/TOL is considered */
/*                     to be of full rank.  If the user sets TOL <= 0, */
/*                     then an implicitly computed, default tolerance, */
/*                     defined by TOLDEF = N*EPS,  is used instead, */
/*                     where EPS is the machine precision (see LAPACK */
/*                     Library routine DLAMCH). */
/*                     This parameter is not relevant if COND = 'U' */
/*                     or 'N'. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine LMPARM. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine LMPARM).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine LMPARM. The LAPACK Library routine */
/*                     XERBLA should be used in conjunction with negative */
/*                     INFO. INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of functions.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of variables.  M >= N >= 0. */

/*     ITMAX   (input) INTEGER */
/*             The maximum number of iterations.  ITMAX >= 0. */

/*     FACTOR  (input) DOUBLE PRECISION */
/*             The value used in determining the initial step bound. This */
/*             bound is set to the product of FACTOR and the Euclidean */
/*             norm of DIAG*X if nonzero, or else to FACTOR itself. */
/*             In most cases FACTOR should lie in the interval (.1,100). */
/*             A generally recommended value is 100.  FACTOR > 0. */

/*     NPRINT  (input) INTEGER */
/*             This parameter enables controlled printing of iterates if */
/*             it is positive. In this case, FCN is called with IFLAG = 0 */
/*             at the beginning of the first iteration and every NPRINT */
/*             iterations thereafter and immediately prior to return, */
/*             with X, E, and J available for printing. Note that when */
/*             called immediately prior to return, J normally contains */
/*             the result returned by QRFACT and LMPARM (the compressed */
/*             R and S factors). If NPRINT is not positive, no special */
/*             calls of FCN with IFLAG = 0 are made. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed, for instance, for */
/*             describing the structure of the Jacobian matrix, which */
/*             are handed over to the routines FCN, QRFACT and LMPARM. */
/*             The first five entries of this array are modified */
/*             internally by a call to FCN (with IFLAG = 3), but are */
/*             restored on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 5. */

/*     DPAR1   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR1,*) or (LDPAR1) */
/*             A first set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03BD */
/*             routine, but it is passed to the routine FCN. */

/*     LDPAR1  (input) INTEGER */
/*             The leading dimension or the length of the array DPAR1, as */
/*             convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, if leading */
/*             dimension.) */

/*     DPAR2   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR2,*) or (LDPAR2) */
/*             A second set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03BD */
/*             routine, but it is passed to the routine FCN. */

/*     LDPAR2  (input) INTEGER */
/*             The leading dimension or the length of the array DPAR2, as */
/*             convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1, if leading */
/*             dimension.) */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if XINIT = 'G', this array must contain the */
/*             vector of initial variables x to be optimized. */
/*             If XINIT = 'R', this array need not be set before entry, */
/*             and random values will be used to initialize x. */
/*             On exit, if INFO = 0, this array contains the vector of */
/*             values that (approximately) minimize the sum of squares of */
/*             error functions. The values returned in IWARN and */
/*             DWORK(1:4) give details on the iterative process. */

/*     DIAG    (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if SCALE = 'S', this array must contain some */
/*             positive entries that serve as multiplicative scale */
/*             factors for the variables x.  DIAG(I) > 0, I = 1,...,N. */
/*             If SCALE = 'I', DIAG is internally set. */
/*             On exit, this array contains the scale factors used */
/*             (or finally used, if SCALE = 'I'). */

/*     NFEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 1. If FCN is */
/*             properly implemented, this includes the function */
/*             evaluations needed for finite difference approximation */
/*             of the Jacobian. */

/*     NJEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 2. */

/*     Tolerances */

/*     FTOL    DOUBLE PRECISION */
/*             If FTOL >= 0, the tolerance which measures the relative */
/*             error desired in the sum of squares. Termination occurs */
/*             when both the actual and predicted relative reductions in */
/*             the sum of squares are at most FTOL. If the user sets */
/*             FTOL < 0,  then  SQRT(EPS)  is used instead FTOL, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     XTOL    DOUBLE PRECISION */
/*             If XTOL >= 0, the tolerance which measures the relative */
/*             error desired in the approximate solution. Termination */
/*             occurs when the relative error between two consecutive */
/*             iterates is at most XTOL. If the user sets  XTOL < 0, */
/*             then  SQRT(EPS)  is used instead XTOL. */

/*     GTOL    DOUBLE PRECISION */
/*             If GTOL >= 0, the tolerance which measures the */
/*             orthogonality desired between the function vector e and */
/*             the columns of the Jacobian J. Termination occurs when */
/*             the cosine of the angle between e and any column of the */
/*             Jacobian J is at most GTOL in absolute value. If the user */
/*             sets  GTOL < 0,  then  EPS  is used instead GTOL. */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the matrices of linear systems to be solved. If */
/*             the user sets TOL > 0, then the given value of TOL is used */
/*             as a lower bound for the reciprocal condition number;  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*EPS,  is used instead. */
/*             This parameter is not relevant if COND = 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+r), where r is the number */
/*             of diagonal blocks R_k in R (see description of LMPARM). */
/*             On output, if INFO = 0, the first N entries of this array */
/*             define a permutation matrix P such that J*P = Q*R, where */
/*             J is the final calculated Jacobian, Q is an orthogonal */
/*             matrix (not stored), and R is upper triangular with */
/*             diagonal elements of nonincreasing magnitude (possibly */
/*             for each block column of J). Column j of P is column */
/*             IWORK(j) of the identity matrix. If INFO = 0, the entries */
/*             N+1:N+r of this array contain the ranks of the final */
/*             submatrices S_k (see description of LMPARM). */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the residual error norm (the */
/*             sum of squares), DWORK(3) returns the number of iterations */
/*             performed, and DWORK(4) returns the final Levenberg */
/*             factor. If INFO = 0, N > 0, and IWARN >= 0, the elements */
/*             DWORK(5) to DWORK(4+M) contain the final matrix-vector */
/*             product Z*Q'*e, and the elements DWORK(5+M) to */
/*             DWORK(4+M+N*NC) contain the (compressed) representation of */
/*             final upper triangular matrices R and S (if IWARN <> 4). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 4, M + max( size(J) + */
/*                                        max( DW( FCN|IFLAG = 1 ), */
/*                                             DW( FCN|IFLAG = 2 ), */
/*                                             DW( QRFACT ) + N ), */
/*                                        N*NC + N + */
/*                                        max( M + DW( FCN|IFLAG = 1 ), */
/*                                             N + DW( LMPARM ) ) ) ), */
/*             where size(J) is the size of the Jacobian (provided by FCN */
/*             in IPAR(1), for IFLAG = 3), and DW( f ) is the workspace */
/*             needed by the routine f, where f is FCN, QRFACT, or LMPARM */
/*             (provided by FCN in IPAR(2:5), for IFLAG = 3). */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             < 0:  the user set IFLAG = IWARN in the subroutine FCN; */
/*             = 1:  both actual and predicted relative reductions in */
/*                   the sum of squares are at most FTOL; */
/*             = 2:  relative error between two consecutive iterates is */
/*                   at most XTOL; */
/*             = 3:  conditions for IWARN = 1 and IWARN = 2 both hold; */
/*             = 4:  the cosine of the angle between e and any column of */
/*                   the Jacobian is at most GTOL in absolute value; */
/*             = 5:  the number of iterations has reached ITMAX without */
/*                   satisfying any convergence condition; */
/*             = 6:  FTOL is too small: no further reduction in the sum */
/*                   of squares is possible; */
/*             = 7:  XTOL is too small: no further improvement in the */
/*                   approximate solution x is possible; */
/*             = 8:  GTOL is too small: e is orthogonal to the columns of */
/*                   the Jacobian to machine precision. */
/*             In all these cases, DWORK(1:4) are set as described above. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 1; */
/*             = 2:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 2; */
/*             = 3:  user-defined routine QRFACT returned with INFO <> 0; */
/*             = 4:  user-defined routine LMPARM returned with INFO <> 0. */

/*     METHOD */

/*     If XINIT = 'R', the initial value for x is set to a vector of */
/*     pseudo-random values uniformly distributed in (-1,1). */

/*     The Levenberg-Marquardt algorithm (described in [1,3]) is used for */
/*     optimizing the variables x. This algorithm needs the Jacobian */
/*     matrix J, which is provided by the subroutine FCN. A trust region */
/*     method is used. The algorithm tries to update x by the formula */

/*         x = x - p, */

/*     using an approximate solution of the system of linear equations */

/*         (J'*J + PAR*D*D)*p = J'*e, */

/*     with e the error function vector, and D a diagonal nonsingular */
/*     matrix, where either PAR = 0 and */

/*         ( norm( D*x ) - DELTA ) <= 0.1*DELTA , */

/*     or PAR > 0 and */

/*         ABS( norm( D*x ) - DELTA ) <= 0.1*DELTA . */

/*     DELTA is the radius of the trust region. If the Gauss-Newton */
/*     direction is not acceptable, then an iterative algorithm obtains */
/*     improved lower and upper bounds for the Levenberg-Marquardt */
/*     parameter PAR. Only a few iterations are generally needed for */
/*     convergence of the algorithm. The trust region radius DELTA */
/*     and the Levenberg factor PAR are updated based on the ratio */
/*     between the actual and predicted reduction in the sum of squares. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     [2] Golub, G.H. and van Loan, C.F. */
/*         Matrix Computations. Third Edition. */
/*         M. D. Johns Hopkins University Press, Baltimore, pp. 520-528, */
/*         1996. */

/*     [3] More, J.J. */
/*         The Levenberg-Marquardt algorithm: implementation and theory. */
/*         In Watson, G.A. (Ed.), Numerical Analysis, Lecture Notes in */
/*         Mathematics, vol. 630, Springer-Verlag, Berlin, Heidelberg */
/*         and New York, pp. 105-116, 1978. */

/*     NUMERICAL ASPECTS */

/*     The Levenberg-Marquardt algorithm described in [3] is scaling */
/*     invariant and globally convergent to (maybe local) minima. */
/*     The convergence rate near a local minimum is quadratic, if the */
/*     Jacobian is computed analytically, and linear, if the Jacobian */
/*     is computed numerically. */

/*     FURTHER COMMENTS */

/*     This routine is a more general version of the subroutines LMDER */
/*     and LMDER1 from the MINPACK package [1], which enables to exploit */
/*     the structure of the problem, and optionally use condition */
/*     estimation. Unstructured problems could be solved as well. */

/*     Template SLICOT Library implementations for FCN, QRFACT and */
/*     LMPARM routines are: */
/*     MD03BF, MD03BA, and MD03BB, respectively, for standard problems; */
/*     NF01BF, NF01BS, and NF01BP, respectively, for optimizing the */
/*     parameters of Wiener systems (structured problems). */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Least-squares approximation,  Levenberg-Marquardt algorithm, */
/*     matrix operations, optimization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --dwork;
    --iwork;
    --diag;
    --x;
    --dpar2;
    --dpar1;
    --ipar;

    /* Function Body */
    init = lsame_(xinit, "R", (ftnlen)1, (ftnlen)1);
    iscal = lsame_(scale, "I", (ftnlen)1, (ftnlen)1);
    sscal = lsame_(scale, "S", (ftnlen)1, (ftnlen)1);
    *info = 0;
    *iwarn = 0;
    if (! (init || lsame_(xinit, "G", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (iscal || sscal)) {
	*info = -2;
    } else if (! (lsame_(cond, "E", (ftnlen)1, (ftnlen)1) || lsame_(cond, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -7;
    } else if (*n < 0 || *n > *m) {
	*info = -8;
    } else if (*itmax < 0) {
	*info = -9;
    } else if (*factor <= 0.) {
	*info = -10;
    } else if (*lipar < 5) {
	*info = -13;
    } else if (*ldpar1 < 0) {
	*info = -15;
    } else if (*ldpar2 < 0) {
	*info = -17;
    } else if (*ldwork < 4) {
	*info = -28;
    } else if (sscal) {
	badscl = FALSE_;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    badscl = badscl || diag[j] <= 0.;
/* L10: */
	}

	if (badscl) {
	    *info = -19;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MD03BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *nfev = 0;
    *njev = 0;
    if (*n == 0) {
	dwork[1] = 4.;
	dwork[2] = 0.;
	dwork[3] = 0.;
	dwork[4] = 0.;
	return 0;
    }

/*     Call FCN to get the size of the array J, for storing the Jacobian */
/*     matrix, the leading dimension LDJ and the workspace required */
/*     by FCN for IFLAG = 1 and IFLAG = 2, QRFACT and LMPARM. The */
/*     entries DWORK(1:4) should not be modified by the special call of */
/*     FCN below, if XINIT = 'R' and the values in DWORK(1:4) are */
/*     explicitly desired for initialization of the random number */
/*     generator. */

    iflag = 3;
    iw1 = ipar[1];
    iw2 = ipar[2];
    iw3 = ipar[3];
    jw1 = ipar[4];
    jw2 = ipar[5];

    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[1], &dwork[1], &ldjsav, &dwork[1], 
	    ldwork, &infol);
    sizej = ipar[1];
    lfcn1 = ipar[2];
    lfcn2 = ipar[3];
    lqrf = ipar[4];
    llmp = ipar[5];
    if (ldjsav > 0) {
	nc = sizej / ldjsav;
    } else {
	nc = sizej;
    }

    ipar[1] = iw1;
    ipar[2] = iw2;
    ipar[3] = iw3;
    ipar[4] = jw1;
    ipar[5] = jw2;

/*     Check the workspace length. */

    e = 1;
    jac = e + *m;
    jw1 = jac + sizej;
    jw2 = jw1 + *n;
    iw1 = jac + *n * nc;
    iw2 = iw1 + *n;
    iw3 = iw2 + *n;
    jwork = iw2 + *m;

/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
    i__5 = max(lfcn1,lfcn2), i__6 = *n + lqrf;
/* Computing MAX */
    i__7 = *m + lfcn1, i__8 = *n + llmp;
    i__3 = sizej + max(i__5,i__6), i__4 = *n * nc + *n + max(i__7,i__8);
    i__1 = 4, i__2 = *m + max(i__3,i__4);
    l = max(i__1,i__2);
    if (*ldwork < l) {
	*info = -28;
	i__1 = -(*info);
	xerbla_("MD03BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Set default tolerances. EPSMCH is the machine precision. */

    epsmch = dlamch_("Epsilon", (ftnlen)7);
    ftdef = *ftol;
    xtdef = *xtol;
    gtdef = *gtol;
    toldef = *tol;
/* Computing MIN */
    d__1 = min(ftdef,xtdef), d__1 = min(d__1,gtdef);
    if (min(d__1,toldef) <= 0.) {
	if (ftdef < 0.) {
	    ftdef = sqrt(epsmch);
	}
	if (xtdef < 0.) {
	    xtdef = sqrt(epsmch);
	}
	if (gtdef < 0.) {
	    gtdef = epsmch;
	}
	if (toldef <= 0.) {
	    toldef = (doublereal) (*n) * epsmch;
	}
    }
    wrkopt = 1;

/*     Initialization. */

    if (init) {

/*        SEED is the initial state of the random number generator. */
/*        SEED(4) must be odd. */

	seed[0] = (integer) dwork[1] % 4096;
	seed[1] = (integer) dwork[2] % 4096;
	seed[2] = (integer) dwork[3] % 4096;
	seed[3] = (((integer) dwork[4] << 1) + 1) % 4096;
	dlarnv_(&c__2, seed, n, &x[1]);
    }

/*     Initialize Levenberg-Marquardt parameter and iteration counter. */

    par = 0.;
    iter = 1;

/*     Evaluate the function at the starting point */
/*     and calculate its norm. */
/*     Workspace: need:    M + SIZEJ + LFCN1; */
/*                prefer:  larger. */

    iflag = 1;
    i__1 = *ldwork - jw1 + 1;
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac], &ldj, &dwork[jw1], 
	    &i__1, &infol);

    if (infol != 0) {
	*info = 1;
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
    wrkopt = max(i__1,i__2);
    *nfev = 1;
    fnorm = dnrm2_(m, &dwork[e], &c__1);
    if (iflag < 0 || fnorm == 0.) {
	goto L90;
    }

/*     Beginning of the outer loop. */

L20:

/*        Calculate the Jacobian matrix. */
/*        Workspace: need:    M + SIZEJ + LFCN2; */
/*                   prefer:  larger. */

    ldj = ldjsav;
    iflag = 2;
    i__1 = *ldwork - jw1 + 1;
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac], &ldj, &dwork[jw1], 
	    &i__1, &infol);

    if (infol != 0) {
	*info = 2;
	return 0;
    }
    if (iter == 1) {
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
	wrkopt = max(i__1,i__2);
    }
    if (nfevl > 0) {
	*nfev += nfevl;
    }
    ++(*njev);
    if (iflag < 0) {
	goto L90;
    }

/*        If requested, call FCN to enable printing of iterates. */

    if (*nprint > 0) {
	iflag = 0;
	if ((iter - 1) % *nprint == 0) {
	    i__1 = *ldwork - jw1 + 1;
	    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1]
		    , ldpar2, &x[1], nfev, &dwork[e], &dwork[jac], &ldj, &
		    dwork[jw1], &i__1, &infol);

	    if (iflag < 0) {
		goto L90;
	    }
	}
    }

/*        Compute the QR factorization of the Jacobian. */
/*        Workspace: need:    M + SIZEJ + N + LQRF; */
/*                   prefer:  larger. */

    i__1 = *ldwork - jw2 + 1;
    (*qrfact)(n, &ipar[1], lipar, &fnorm, &dwork[jac], &ldj, &dwork[e], &
	    dwork[jw1], &gnorm, &iwork[1], &dwork[jw2], &i__1, &infol);
    if (infol != 0) {
	*info = 3;
	return 0;
    }

/*        On the first iteration and if SCALE = 'I', scale according */
/*        to the norms of the columns of the initial Jacobian. */

    if (iter == 1) {
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jw2] + jw2 - 1;
	wrkopt = max(i__1,i__2);
	if (iscal) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		diag[j] = dwork[jw1 + j - 1];
		if (diag[j] == 0.) {
		    diag[j] = 1.;
		}
/* L30: */
	    }

	}

/*           On the first iteration, calculate the norm of the scaled */
/*           x and initialize the step bound DELTA. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dwork[iw1 + j - 1] = diag[j] * x[j];
/* L40: */
	}

	xnorm = dnrm2_(n, &dwork[iw1], &c__1);
	delta = *factor * xnorm;
	if (delta == 0.) {
	    delta = *factor;
	}
    } else {

/*           Rescale if necessary. */

	if (iscal) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		d__1 = diag[j], d__2 = dwork[jw1 + j - 1];
		diag[j] = max(d__1,d__2);
/* L50: */
	    }

	}
    }

/*        Test for convergence of the gradient norm. */

    if (gnorm <= gtdef) {
	*iwarn = 4;
    }
    if (*iwarn != 0) {
	goto L90;
    }

/*        Beginning of the inner loop. */

L60:

/*           Determine the Levenberg-Marquardt parameter and the */
/*           direction p, and compute -R*P'*p. */
/*           Workspace:  need:    M + N*NC + 2*N + LLMP; */
/*                       prefer:  larger. */

    i__1 = *ldwork - iw3 + 1;
    (*lmparm)(cond, n, &ipar[1], lipar, &dwork[jac], &ldj, &iwork[1], &diag[1]
	    , &dwork[e], &delta, &par, &iwork[*n + 1], &dwork[iw1], &dwork[
	    iw2], &toldef, &dwork[iw3], &i__1, &infol, (ftnlen)1);
    if (infol != 0) {
	*info = 4;
	return 0;
    }
    if (iter == 1) {
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[iw3] + iw3 - 1;
	wrkopt = max(i__1,i__2);
    }

    temp1 = dnrm2_(n, &dwork[iw2], &c__1) / fnorm;

/*           Store the direction p and x - p. */

    i__1 = *n - 1;
    for (j = 0; j <= i__1; ++j) {
	dwork[iw2 + j] = diag[j + 1] * dwork[iw1 + j];
	dwork[iw1 + j] = x[j + 1] - dwork[iw1 + j];
/* L70: */
    }

/*           Compute the norm of scaled p and the scaled predicted */
/*           reduction and the scaled directional derivative. */

    pnorm = dnrm2_(n, &dwork[iw2], &c__1);
    temp2 = sqrt(par) * pnorm / fnorm;
/* Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    prered = d__1 * d__1 + d__2 * d__2 / .5;
/* Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    dirder = -(d__1 * d__1 + d__2 * d__2);

/*           On the first iteration, adjust the initial step bound. */

    if (iter == 1) {
	delta = min(delta,pnorm);
    }

/*           Evaluate the function at x - p and calculate its norm. */
/*           Workspace:  need:    2*M + N*NC + N + LFCN1; */
/*                       prefer:  larger. */

    iflag = 1;
    i__1 = *ldwork - jwork + 1;
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &dwork[iw1], &nfevl, &dwork[iw2], &dwork[jac], &ldj, &
	    dwork[jwork], &i__1, &infol);
    if (infol != 0) {
	*info = 1;
	return 0;
    }

    ++(*nfev);
    if (iflag < 0) {
	goto L90;
    }
    fnorm1 = dnrm2_(m, &dwork[iw2], &c__1);

/*           Compute the scaled actual reduction. */

    actred = -1.;
    if (fnorm1 * .1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = 1. - d__1 * d__1;
    }

/*           Compute the ratio of the actual to the predicted reduction. */

    ratio = 0.;
    if (prered != 0.) {
	ratio = actred / prered;
    }

/*           Update the step bound. */

    if (ratio <= .25) {
	if (actred >= 0.) {
	    temp = .5;
	} else {
	    temp = dirder * .5 / (dirder + actred * .5);
	}
	if (fnorm1 * .1 >= fnorm || temp < .1) {
	    temp = .1;
	}
/* Computing MIN */
	d__1 = delta, d__2 = pnorm / .1;
	delta = temp * min(d__1,d__2);
	par /= temp;
    } else {
	if (par == 0. || ratio >= .75) {
	    delta = pnorm / .5;
	    par *= .5;
	}
    }

/*           Test for successful iteration. */

    if (ratio >= 1e-4) {

/*              Successful iteration. Update x, e, and their norms. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    x[j] = dwork[iw1 + j - 1];
	    dwork[iw1 + j - 1] = diag[j] * x[j];
/* L80: */
	}

	dcopy_(m, &dwork[iw2], &c__1, &dwork[e], &c__1);
	xnorm = dnrm2_(n, &dwork[iw1], &c__1);
	fnorm = fnorm1;
	++iter;
    }

/*           Tests for convergence. */

    if (abs(actred) <= ftdef && prered <= ftdef && ratio * .5 <= 1.) {
	*iwarn = 1;
    }
    if (delta <= xtdef * xnorm) {
	*iwarn = 2;
    }
    if (abs(actred) <= ftdef && prered <= ftdef && ratio * .5 <= 1. && *iwarn 
	    == 2) {
	*iwarn = 3;
    }
    if (*iwarn != 0) {
	goto L90;
    }

/*           Tests for termination and stringent tolerances. */

    if (iter >= *itmax) {
	*iwarn = 5;
    }
    if (abs(actred) <= epsmch && prered <= epsmch && ratio * .5 <= 1.) {
	*iwarn = 6;
    }
    if (delta <= epsmch * xnorm) {
	*iwarn = 7;
    }
    if (gnorm <= epsmch) {
	*iwarn = 8;
    }
    if (*iwarn != 0) {
	goto L90;
    }

/*           End of the inner loop. Repeat if unsuccessful iteration. */

    if (ratio < 1e-4) {
	goto L60;
    }

/*        End of the outer loop. */

    goto L20;

L90:

/*     Termination, either normal or user imposed. */
/*     Note that DWORK(JAC) normally contains the results returned by */
/*     QRFACT and LMPARM (the compressed R and S factors). */

    if (iflag < 0) {
	*iwarn = iflag;
    }
    if (*nprint > 0) {
	iflag = 0;
	i__1 = *ldwork - jwork + 1;
	(*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
		ldpar2, &x[1], nfev, &dwork[e], &dwork[jac], &ldj, &dwork[
		jwork], &i__1, &infol);
	if (iflag < 0) {
	    *iwarn = iflag;
	}
    }

    if (*iwarn >= 0) {
	for (j = *m + *n * nc; j >= 1; --j) {
	    dwork[j + 4] = dwork[j];
/* L100: */
	}
    }
    dwork[1] = (doublereal) wrkopt;
    dwork[2] = fnorm;
    dwork[3] = (doublereal) iter;
    dwork[4] = par;

    return 0;
/* *** Last line of MD03BD *** */
} /* md03bd_ */

