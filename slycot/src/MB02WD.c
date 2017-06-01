/* MB02WD.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b8 = 1.;
static doublereal c_b9 = -1.;
static doublereal c_b22 = 0.;

/* Subroutine */ int mb02wd_(char *form, S_fp f, integer *n, integer *ipar, 
	integer *lipar, doublereal *dpar, integer *ldpar, integer *itmax, 
	doublereal *a, integer *lda, doublereal *b, integer *incb, doublereal 
	*x, integer *incx, doublereal *tol, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen form_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer k, r__, aq;
    static logical mat;
    static doublereal res, beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), dsymv_(char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer dwleft;
    static doublereal resold;


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

/*     To solve the system of linear equations Ax = b, with A symmetric, */
/*     positive definite, or, in the implicit form, f(A, x) = b, where */
/*     y = f(A, x) is a symmetric positive definite linear mapping */
/*     from x to y, using the conjugate gradient (CG) algorithm without */
/*     preconditioning. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FORM     CHARACTER*1 */
/*              Specifies the form of the system of equations, as */
/*              follows: */
/*              = 'U' :  Ax = b, the upper triagular part of A is used; */
/*              = 'L' :  Ax = b, the lower triagular part of A is used; */
/*              = 'F' :  the implicit, function form, f(A, x) = b. */

/*     Function Parameters */

/*     F       EXTERNAL */
/*             If FORM = 'F', then F is a subroutine which calculates the */
/*             value of f(A, x), for given A and x. */
/*             If FORM <> 'F', then F is not called. */

/*             F must have the following interface: */

/*             SUBROUTINE F( N, IPAR, LIPAR, DPAR, LDPAR, A, LDA, X, */
/*            $              INCX, DWORK, LDWORK, INFO ) */

/*             where */

/*             N       (input) INTEGER */
/*                     The dimension of the vector x.  N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the matrix A. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*                     The real parameters needed for solving the */
/*                     problem. */

/*             LDPAR   (input) INTEGER */
/*                     The length of the array DPAR.  LDPAR >= 0. */

/*             A       (input) DOUBLE PRECISION array, dimension */
/*                     (LDA, NC), where NC is the number of columns. */
/*                     The leading NR-by-NC part of this array must */
/*                     contain the (compressed) representation of the */
/*                     matrix A, where NR is the number of rows of A */
/*                     (function of IPAR entries). */

/*             LDA     (input) INTEGER */
/*                     The leading dimension of the array A. */
/*                     LDA >= MAX(1,NR). */

/*             X       (input/output) DOUBLE PRECISION array, dimension */
/*                     (1+(N-1)*INCX) */
/*                     On entry, this incremented array must contain the */
/*                     vector x. */
/*                     On exit, this incremented array contains the value */
/*                     of the function f, y = f(A, x). */

/*             INCX    (input) INTEGER */
/*                     The increment for the elements of X.  INCX > 0. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine F. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine F). */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input scalar argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine F. The LAPACK Library routine XERBLA */
/*                     should be used in conjunction with negative INFO. */
/*                     INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the vector x.  N >= 0. */
/*             If FORM = 'U' or FORM = 'L', N is also the number of rows */
/*             and columns of the matrix A. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             If FORM = 'F', the integer parameters describing the */
/*             structure of the matrix A. */
/*             This parameter is ignored if FORM = 'U' or FORM = 'L'. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 0. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             If FORM = 'F', the real parameters needed for solving */
/*             the problem. */
/*             This parameter is ignored if FORM = 'U' or FORM = 'L'. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 0. */

/*     ITMAX   (input) INTEGER */
/*             The maximal number of iterations to do.  ITMAX >= 0. */

/*     A       (input) DOUBLE PRECISION array, */
/*                     dimension (LDA, NC), if FORM = 'F', */
/*                     dimension (LDA, N),  otherwise. */
/*             If FORM = 'F', the leading NR-by-NC part of this array */
/*             must contain the (compressed) representation of the */
/*             matrix A, where NR and NC are the number of rows and */
/*             columns, respectively, of the matrix A. The array A is */
/*             not referenced by this routine itself, except in the */
/*             calls to the routine F. */
/*             If FORM <> 'F', the leading N-by-N part of this array */
/*             must contain the matrix A, assumed to be symmetric; */
/*             only the triangular part specified by FORM is referenced. */

/*     LDA     (input) INTEGER */
/*             The leading dimension of array A. */
/*             LDA >= MAX(1,NR), if FORM = 'F'; */
/*             LDA >= MAX(1,N),  if FORM = 'U' or FORM = 'L'. */

/*     B       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCB) */
/*             The incremented vector b. */

/*     INCB    (input) INTEGER */
/*             The increment for the elements of B.  INCB > 0. */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*INCX) */
/*             On entry, this incremented array must contain an initial */
/*             approximation of the solution. If an approximation is not */
/*             known, setting all elements of x to zero is recommended. */
/*             On exit, this incremented array contains the computed */
/*             solution x of the system of linear equations. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X.  INCX > 0. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If TOL > 0, absolute tolerance for the iterative process. */
/*             The algorithm will stop if || Ax - b ||_2 <= TOL. Since */
/*             it is advisable to use a relative tolerance, say TOLER, */
/*             TOL should be chosen as TOLER*|| b ||_2. */
/*             If TOL <= 0, a default relative tolerance, */
/*             TOLDEF = N*EPS*|| b ||_2,  is used, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the number of */
/*             iterations performed and DWORK(2) returns the remaining */
/*             residual, || Ax - b ||_2. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(2,3*N + DWORK(F)),  if FORM = 'F', */
/*                       where DWORK(F) is the workspace needed by F; */
/*             LDWORK >= MAX(2,3*N),       if FORM = 'U' or FORM = 'L'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the algorithm finished after ITMAX > 0 iterations, */
/*                   without achieving the desired precision TOL; */
/*             = 2:  ITMAX is zero; in this case, DWORK(2) is not set. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then F returned with INFO = i. */

/*     METHOD */

/*     The following CG iteration is used for solving Ax = b: */

/*     Start: q(0) = r(0) = Ax - b */

/*                   < q(k),  r(k) > */
/*     ALPHA(k) = - ---------------- */
/*                   < q(k), Aq(k) > */
/*     x(k+1)   = x(k) - ALPHA(k) * q(k) */
/*     r(k+1)   = r(k) - ALPHA(k) * Aq(k) */
/*                 < r(k+1), r(k+1) > */
/*     BETA(k)  = -------------------- */
/*                 < r(k)  , r(k)   > */
/*     q(k+1)   = r(k+1) + BETA(k) * q(k) */

/*     where <.,.> denotes the scalar product. */

/*     REFERENCES */

/*     [1] Golub, G.H. and van Loan, C.F. */
/*         Matrix Computations. Third Edition. */
/*         M. D. Johns Hopkins University Press, Baltimore, pp. 520-528, */
/*         1996. */

/*     [2] Luenberger, G. */
/*         Introduction to Linear and Nonlinear Programming. */
/*         Addison-Wesley, Reading, MA, p.187, York, 1973. */

/*     NUMERICAL ASPECTS */

/*     Since the residuals are orthogonal in the scalar product */
/*     <x, y> = y'Ax, the algorithm is theoretically finite. But rounding */
/*     errors cause a loss of orthogonality, so a finite termination */
/*     cannot be guaranteed. However, one can prove [2] that */

/*        || x-x_k ||_A := sqrt( (x-x_k)' * A * (x-x_k) ) */

/*                                             sqrt( kappa_2(A) ) - 1 */
/*                      <=  2 || x-x_0 ||_A * ------------------------ , */
/*                                             sqrt( kappa_2(A) ) + 1 */

/*     where kappa_2 is the condition number. */

/*     The approximate number of floating point operations is */
/*        (k*(N**2 + 15*N) + N**2 + 3*N)/2, if FORM <> 'F', */
/*        k*(f + 7*N) + f,                  if FORM =  'F', */
/*     where k is the number of CG iterations performed, and f is the */
/*     number of floating point operations required by the subroutine F. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     March, 2002. */

/*     KEYWORDS */

/*     Conjugate gradients, convergence, linear system of equations, */
/*     matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    --ipar;
    --dpar;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --x;
    --dwork;

    /* Function Body */
    mat = lsame_(form, "U", (ftnlen)1, (ftnlen)1) || lsame_(form, "L", (
	    ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    *iwarn = 0;
    *info = 0;
    if (! (mat || lsame_(form, "F", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -3;
    } else if (! mat && *lipar < 0) {
	*info = -5;
    } else if (! mat && *ldpar < 0) {
	*info = -7;
    } else if (*itmax < 0) {
	*info = -8;
    } else if (*lda < 1 || mat && *lda < *n) {
	*info = -10;
    } else if (*incb <= 0) {
	*info = -12;
    } else if (*incx <= 0) {
	*info = -14;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 2, i__2 = *n * 3;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -17;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02WD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 0.;
	dwork[2] = 0.;
	return 0;
    }

    if (*itmax == 0) {
	dwork[1] = 0.;
	*iwarn = 2;
	return 0;
    }

/*     Set default tolerance, if needed. */

    toldef = *tol;
    if (toldef <= 0.) {
	toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7) * dnrm2_(n,
		 &b[1], incb);
    }

/*     Initialize local variables. */

    k = 0;

/*     Vector q is stored in DWORK(1), A*q or f(A, q) in DWORK(AQ), */
/*     and r in DWORK(R). The workspace for F starts in DWORK(DWLEFT). */

    aq = *n + 1;
    r__ = *n + aq;
    dwleft = *n + r__;

/*     Prepare the first iteration, initialize r and q. */

    if (mat) {
	dcopy_(n, &b[1], incb, &dwork[r__], &c__1);
	dsymv_(form, n, &c_b8, &a[a_offset], lda, &x[1], incx, &c_b9, &dwork[
		r__], &c__1, (ftnlen)1);
    } else {
	dcopy_(n, &x[1], incx, &dwork[r__], &c__1);
	i__1 = *ldwork - dwleft + 1;
	(*f)(n, &ipar[1], lipar, &dpar[1], ldpar, &a[a_offset], lda, &dwork[
		r__], &c__1, &dwork[dwleft], &i__1, info);
	if (*info != 0) {
	    return 0;
	}
	daxpy_(n, &c_b9, &b[1], incb, &dwork[r__], &c__1);
    }
    dcopy_(n, &dwork[r__], &c__1, &dwork[1], &c__1);

    res = dnrm2_(n, &dwork[r__], &c__1);

/*     Do nothing if x is already the solution. */

    if (res <= toldef) {
	goto L20;
    }

/*     Begin of the iteration loop. */

/*     WHILE ( RES.GT.TOLDEF .AND. K.LE.ITMAX ) DO */
L10:

/*        Calculate A*q or f(A, q). */

    if (mat) {
	dsymv_(form, n, &c_b8, &a[a_offset], lda, &dwork[1], &c__1, &c_b22, &
		dwork[aq], &c__1, (ftnlen)1);
    } else {
	dcopy_(n, &dwork[1], &c__1, &dwork[aq], &c__1);
	i__1 = *ldwork - dwleft + 1;
	(*f)(n, &ipar[1], lipar, &dpar[1], ldpar, &a[a_offset], lda, &dwork[
		aq], &c__1, &dwork[dwleft], &i__1, info);
	if (*info != 0) {
	    return 0;
	}
    }

/*        Calculate ALPHA(k). */

    alpha = ddot_(n, &dwork[1], &c__1, &dwork[r__], &c__1) / ddot_(n, &dwork[
	    1], &c__1, &dwork[aq], &c__1);

/*        x(k+1) = x(k) - ALPHA(k)*q(k). */

    d__1 = -alpha;
    daxpy_(n, &d__1, &dwork[1], &c__1, &x[1], incx);

/*        r(k+1) = r(k) - ALPHA(k)*(A*q(k)). */

    d__1 = -alpha;
    daxpy_(n, &d__1, &dwork[aq], &c__1, &dwork[r__], &c__1);

/*        Save RES and calculate a new RES. */

    resold = res;
    res = dnrm2_(n, &dwork[r__], &c__1);

/*        Exit if tolerance is reached. */

    if (res <= toldef) {
	goto L20;
    }

/*        Calculate BETA(k). */

/* Computing 2nd power */
    d__1 = res / resold;
    beta = d__1 * d__1;

/*        q(k+1) = r(k+1) + BETA(k)*q(k). */

    dscal_(n, &beta, &dwork[1], &c__1);
    daxpy_(n, &c_b8, &dwork[r__], &c__1, &dwork[1], &c__1);

/*        End of the iteration loop. */

    ++k;
    if (k < *itmax) {
	goto L10;
    }
/*     END WHILE 10 */

/*     Tolerance was not reached! */

    *iwarn = 1;

L20:

    dwork[1] = (doublereal) k;
    dwork[2] = res;

/* *** Last line of MB02WD *** */
    return 0;
} /* mb02wd_ */

