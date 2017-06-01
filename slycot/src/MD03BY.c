/* MD03BY.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int md03by_(char *cond, integer *n, doublereal *r__, integer 
	*ldr, integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *
	delta, doublereal *par, integer *rank, doublereal *x, doublereal *rx, 
	doublereal *tol, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, l, n2;
    static doublereal fp, dum[3], parc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal parl;
    static logical sing;
    static integer iter;
    static doublereal temp, paru;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd;
    extern /* Subroutine */ int mb02yd_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static char condl[1];
    static logical ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwarf, dmino;
    static logical ucond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal gnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dtrsv_(char *, char *, char *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen,
	     ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal dxnorm;


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

/*     To determine a value for the parameter PAR such that if x solves */
/*     the system */

/*           A*x = b ,     sqrt(PAR)*D*x = 0 , */

/*     in the least squares sense, where A is an m-by-n matrix, D is an */
/*     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if */
/*     DELTA is a positive number, DXNORM is the Euclidean norm of D*x, */
/*     then either PAR is zero and */

/*           ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*     It is assumed that a QR factorization, with column pivoting, of A */
/*     is available, that is, A*P = Q*R, where P is a permutation matrix, */
/*     Q has orthogonal columns, and R is an upper triangular matrix */
/*     with diagonal elements of nonincreasing magnitude. */
/*     The routine needs the full upper triangle of R, the permutation */
/*     matrix P, and the first n components of Q'*b (' denotes the */
/*     transpose). On output, MD03BY also provides an upper triangular */
/*     matrix S such that */

/*           P'*(A'*A + PAR*D*D)*P = S'*S . */

/*     Matrix S is used in the solution process. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the matrices R and S */
/*             should be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation for R and S; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R and S for zero values; */
/*             = 'U' :  use the rank already stored in RANK (for R). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N >= 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR, N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the full upper triangle is unaltered, and the */
/*             strict lower triangle contains the strict upper triangle */
/*             (transposed) of the upper triangular matrix S. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             A*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             An upper bound on the Euclidean norm of D*x.  DELTA > 0. */

/*     PAR     (input/output) DOUBLE PRECISION */
/*             On entry, PAR must contain an initial estimate of the */
/*             Levenberg-Marquardt parameter.  PAR >= 0. */
/*             On exit, it contains the final estimate of this parameter. */

/*     RANK    (input or output) INTEGER */
/*             On entry, if COND = 'U', this parameter must contain the */
/*             (numerical) rank of the matrix R. */
/*             On exit, this parameter contains the numerical rank of */
/*             the matrix S. */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system A*x = b, sqrt(PAR)*D*x = 0. */

/*     RX      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the matrix-vector product -R*P'*x. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             rank of the matrices R and S. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 4*N, if COND =  'E'; */
/*             LDWORK >= 2*N, if COND <> 'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm computes the Gauss-Newton direction. A least squares */
/*     solution is found if the Jacobian is rank deficient. If the Gauss- */
/*     Newton direction is not acceptable, then an iterative algorithm */
/*     obtains improved lower and upper bounds for the parameter PAR. */
/*     Only a few iterations are generally needed for convergence of the */
/*     algorithm. If, however, the limit of ITMAX = 10 iterations is */
/*     reached, then the output PAR will contain the best value obtained */
/*     so far. If the Gauss-Newton step is acceptable, it is stored in x, */
/*     and PAR is set to zero, hence S = R. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     This routine is a LAPACK-based modification of LMPAR from the */
/*     MINPACK package [1], and with optional condition estimation. */
/*     The option COND = 'U' is useful when dealing with several */
/*     right-hand side vectors, but RANK should be reset. */
/*     If COND = 'E', but the matrix S is guaranteed to be nonsingular */
/*     and well conditioned relative to TOL, i.e., rank(R) = N, and */
/*     min(DIAG) > 0, then its condition is not estimated. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

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
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --ipvt;
    --diag;
    --qtb;
    --x;
    --rx;
    --dwork;

    /* Function Body */
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
    ucond = lsame_(cond, "U", (ftnlen)1, (ftnlen)1);
    *info = 0;
    if (! (econd || ncond || ucond)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldr < max(1,*n)) {
	*info = -4;
    } else if (*delta <= 0.) {
	*info = -8;
    } else if (*par < 0.) {
	*info = -9;
    } else if (ucond && (*rank < 0 || *rank > *n)) {
	*info = -10;
    } else if (*ldwork < *n << 1 || econd && *ldwork < *n << 2) {
	*info = -15;
    } else if (*n > 0) {
	dmino = diag[1];
	sing = FALSE_;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (diag[j] < dmino) {
		dmino = diag[j];
	    }
	    sing = sing || diag[j] == 0.;
/* L10: */
	}

	if (sing) {
	    *info = -6;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MD03BY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*par = 0.;
	*rank = 0;
	return 0;
    }

/*     DWARF is the smallest positive magnitude. */

    dwarf = dlamch_("Underflow", (ftnlen)9);
    n2 = *n;

/*     Estimate the rank of R, if required. */

    if (econd) {
	n2 = *n << 1;
	temp = *tol;
	if (temp <= 0.) {

/*           Use the default tolerance in rank determination. */

	    temp = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	}

/*        Estimate the reciprocal condition number of R and set the rank. */
/*        Workspace: 2*N. */

	mb03od_("No QR", n, n, &r__[r_offset], ldr, &ipvt[1], &temp, &c_b10, &
		dwork[1], rank, dum, &dwork[1], ldwork, info, (ftnlen)5);

    } else if (ncond) {
	j = 1;

L20:
	if (r__[j + j * r_dim1] != 0.) {
	    ++j;
	    if (j <= *n) {
		goto L20;
	    }
	}

	*rank = j - 1;
    }

/*     Compute and store in x the Gauss-Newton direction. If the */
/*     Jacobian is rank-deficient, obtain a least squares solution. */
/*     The array RX is used as workspace. */

    dcopy_(rank, &qtb[1], &c__1, &rx[1], &c__1);
    dum[0] = 0.;
    if (*rank < *n) {
	i__1 = *n - *rank;
	dcopy_(&i__1, dum, &c__0, &rx[*rank + 1], &c__1);
    }
    dtrsv_("Upper", "No transpose", "Non unit", rank, &r__[r_offset], ldr, &
	    rx[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = rx[j];
/* L30: */
    }

/*     Initialize the iteration counter. */
/*     Evaluate the function at the origin, and test */
/*     for acceptance of the Gauss-Newton direction. */

    iter = 0;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dwork[j] = diag[j] * x[j];
/* L40: */
    }

    dxnorm = dnrm2_(n, &dwork[1], &c__1);
    fp = dxnorm - *delta;
    if (fp > *delta * .1) {

/*        Set an appropriate option for estimating the condition of */
/*        the matrix S. */

	if (ucond) {
	    if (*ldwork >= *n << 2) {
		*(unsigned char *)condl = 'E';
		toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	    } else {
		*(unsigned char *)condl = 'N';
		toldef = *tol;
	    }
	} else {
	    *(unsigned char *)condl = *(unsigned char *)cond;
	    toldef = *tol;
	}

/*        If the Jacobian is not rank deficient, the Newton */
/*        step provides a lower bound, PARL, for the zero of */
/*        the function. Otherwise set this bound to zero. */

	if (*rank == *n) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		rx[j] = diag[l] * (dwork[l] / dxnorm);
/* L50: */
	    }

	    dtrsv_("Upper", "Transpose", "Non unit", n, &r__[r_offset], ldr, &
		    rx[1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
	    temp = dnrm2_(n, &rx[1], &c__1);
	    parl = fp / *delta / temp / temp;

/*           For efficiency, use CONDL = 'U', if possible. */

	    if (! lsame_(condl, "U", (ftnlen)1, (ftnlen)1) && dmino > 0.) {
		*(unsigned char *)condl = 'U';
	    }
	} else {
	    parl = 0.;
	}

/*        Calculate an upper bound, PARU, for the zero of the function. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    l = ipvt[j];
	    rx[j] = ddot_(&j, &r__[j * r_dim1 + 1], &c__1, &qtb[1], &c__1) / 
		    diag[l];
/* L60: */
	}

	gnorm = dnrm2_(n, &rx[1], &c__1);
	paru = gnorm / *delta;
	if (paru == 0.) {
	    paru = dwarf / min(*delta,.1) / .001;
	}

/*        If the input PAR lies outside of the interval (PARL,PARU), */
/*        set PAR to the closer endpoint. */

	*par = max(*par,parl);
	*par = min(*par,paru);
	if (*par == 0.) {
	    *par = gnorm / dxnorm;
	}

/*        Beginning of an iteration. */

L70:
	++iter;

/*           Evaluate the function at the current value of PAR. */

	if (*par == 0.) {
/* Computing MAX */
	    d__1 = dwarf, d__2 = paru * .001;
	    *par = max(d__1,d__2);
	}
	temp = sqrt(*par);

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    rx[j] = temp * diag[j];
/* L80: */
	}

/*           Solve the system A*x = b , sqrt(PAR)*D*x = 0 , in a least */
/*           square sense. The first N elements of DWORK contain the */
/*           diagonal elements of the upper triangular matrix S, and */
/*           the next N elements contain the vector z, so that x = P*z. */
/*           The vector z is preserved if COND = 'E'. */
/*           Workspace:   4*N, if CONDL =  'E'; */
/*                        2*N, if CONDL <> 'E'. */

	mb02yd_(condl, n, &r__[r_offset], ldr, &ipvt[1], &rx[1], &qtb[1], 
		rank, &x[1], &toldef, &dwork[1], ldwork, info, (ftnlen)1);

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dwork[n2 + j] = diag[j] * x[j];
/* L90: */
	}

	dxnorm = dnrm2_(n, &dwork[n2 + 1], &c__1);
	temp = fp;
	fp = dxnorm - *delta;

/*           If the function is small enough, accept the current value */
/*           of PAR. Also test for the exceptional cases where PARL */
/*           is zero or the number of iterations has reached ITMAX. */

	if (abs(fp) > *delta * .1 && (parl != 0. || fp > temp || temp >= 0.) 
		&& iter < 10) {

/*              Compute the Newton correction. */

	    i__1 = *rank;
	    for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		rx[j] = diag[l] * (dwork[n2 + l] / dxnorm);
/* L100: */
	    }

	    if (*rank < *n) {
		i__1 = *n - *rank;
		dcopy_(&i__1, dum, &c__0, &rx[*rank + 1], &c__1);
	    }
	    i__1 = *ldr + 1;
	    dswap_(n, &r__[r_offset], &i__1, &dwork[1], &c__1);
	    dtrsv_("Lower", "No transpose", "Non Unit", rank, &r__[r_offset], 
		    ldr, &rx[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	    i__1 = *ldr + 1;
	    dswap_(n, &r__[r_offset], &i__1, &dwork[1], &c__1);
	    temp = dnrm2_(rank, &rx[1], &c__1);
	    parc = fp / *delta / temp / temp;

/*              Depending on the sign of the function, update PARL */
/*              or PARU. */

	    if (fp > 0.) {
		parl = max(parl,*par);
	    } else if (fp < 0.) {
		paru = min(paru,*par);
	    }

/*              Compute an improved estimate for PAR. */

/* Computing MAX */
	    d__1 = parl, d__2 = *par + parc;
	    *par = max(d__1,d__2);

/*              End of an iteration. */

	    goto L70;
	}
    }

/*     Compute -R*P'*x = -R*z. */

    if (econd && iter > 0) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    rx[j] = -dwork[*n + j];
/* L110: */
	}

	dtrmv_("Upper", "NoTranspose", "NonUnit", n, &r__[r_offset], ldr, &rx[
		1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
    } else {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    rx[j] = 0.;
	    l = ipvt[j];
	    d__1 = -x[l];
	    daxpy_(&j, &d__1, &r__[j * r_dim1 + 1], &c__1, &rx[1], &c__1);
/* L120: */
	}

    }

/*     Termination. If PAR = 0, set S. */

    if (iter == 0) {
	*par = 0.;

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    dwork[j] = r__[j + j * r_dim1];
	    i__2 = *n - j;
	    dcopy_(&i__2, &r__[j + (j + 1) * r_dim1], ldr, &r__[j + 1 + j * 
		    r_dim1], &c__1);
/* L130: */
	}

	dwork[*n] = r__[*n + *n * r_dim1];
    }

    return 0;

/* *** Last line of MD03BY *** */
} /* md03by_ */

