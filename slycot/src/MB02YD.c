/* MB02YD.f -- translated by f2c (version 20100827).
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
static doublereal c_b19 = 0.;
static integer c__0 = 0;

/* Subroutine */ int mb02yd_(char *cond, integer *n, doublereal *r__, integer 
	*ldr, integer *ipvt, doublereal *diag, doublereal *qtb, integer *rank,
	 doublereal *x, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal cs, sn, dum[3], temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), mb03od_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen);
    static logical econd, ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ucond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal qtbpj;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen);


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

/*     To determine a vector x which solves the system of linear */
/*     equations */

/*           A*x = b ,     D*x = 0 , */

/*     in the least squares sense, where A is an m-by-n matrix, */
/*     D is an n-by-n diagonal matrix, and b is an m-vector. */
/*     It is assumed that a QR factorization, with column pivoting, of A */
/*     is available, that is, A*P = Q*R, where P is a permutation matrix, */
/*     Q has orthogonal columns, and R is an upper triangular matrix */
/*     with diagonal elements of nonincreasing magnitude. */
/*     The routine needs the full upper triangle of R, the permutation */
/*     matrix P, and the first n components of Q'*b (' denotes the */
/*     transpose). The system A*x = b, D*x = 0, is then equivalent to */

/*           R*z = Q'*b ,  P'*D*P*z = 0 ,                             (1) */

/*     where x = P*z. If this system does not have full rank, then a */
/*     least squares solution is obtained. On output, MB02YD also */
/*     provides an upper triangular matrix S such that */

/*           P'*(A'*A + D*D)*P = S'*S . */

/*     The system (1) is equivalent to S*z = c , where c contains the */
/*     first n components of the vector obtained by applying to */
/*     [ (Q'*b)'  0 ]' the transformations which triangularized */
/*     [ R'  P'*D*P ]', getting S. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the matrix S should be */
/*             estimated, as follows: */
/*             = 'E' :  use incremental condition estimation and store */
/*                      the numerical rank of S in RANK; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of S for zero values; */
/*             = 'U' :  use the rank already stored in RANK. */

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
/*             matrix D. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     RANK    (input or output) INTEGER */
/*             On entry, if COND = 'U', this parameter must contain the */
/*             (numerical) rank of the matrix S. */
/*             On exit, if COND = 'E' or 'N', this parameter contains */
/*             the numerical rank of the matrix S, estimated according */
/*             to the value of COND. */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system A*x = b, D*x = 0. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             rank of the matrix S. If the user sets TOL > 0, then the */
/*             given value of TOL is used as a lower bound for the */
/*             reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S, and */
/*             the next N elements contain the solution z. */

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

/*     Standard plane rotations are used to annihilate the elements of */
/*     the diagonal matrix D, updating the upper triangular matrix R */
/*     and the first n elements of the vector Q'*b. A basic least squares */
/*     solution is computed. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     This routine is a LAPACK-based modification of QRSOLV from the */
/*     MINPACK package [1], and with optional condition estimation. */
/*     The option COND = 'U' is useful when dealing with several */
/*     right-hand side vectors. */

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
    } else if (ucond && (*rank < 0 || *rank > *n)) {
	*info = -8;
    } else if (*ldwork < *n << 1 || econd && *ldwork < *n << 2) {
	*info = -12;
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02YD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (! ucond) {
	    *rank = 0;
	}
	return 0;
    }

/*     Copy R and Q'*b to preserve input and initialize S. */
/*     In particular, save the diagonal elements of R in X. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = r__[j + j * r_dim1];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
/* L10: */
	}
/* L20: */
    }

    dcopy_(n, &qtb[1], &c__1, &dwork[*n + 1], &c__1);

/*     Eliminate the diagonal matrix D using Givens rotations. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*        Prepare the row of D to be eliminated, locating the */
/*        diagonal element using P from the QR factorization. */

	l = ipvt[j];
	if (diag[l] != 0.) {
	    qtbpj = 0.;
	    dwork[j] = diag[l];

	    i__2 = *n;
	    for (k = j + 1; k <= i__2; ++k) {
		dwork[k] = 0.;
/* L30: */
	    }

/*           The transformations to eliminate the row of D modify only */
/*           a single element of Q'*b beyond the first n, which is */
/*           initially zero. */

	    i__2 = *n;
	    for (k = j; k <= i__2; ++k) {

/*              Determine a Givens rotation which eliminates the */
/*              appropriate element in the current row of D. */

		if (dwork[k] != 0.) {

		    dlartg_(&r__[k + k * r_dim1], &dwork[k], &cs, &sn, &temp);

/*                 Compute the modified diagonal element of R and */
/*                 the modified elements of (Q'*b,0). */
/*                 Accumulate the tranformation in the row of S. */

		    temp = cs * dwork[*n + k] + sn * qtbpj;
		    qtbpj = -sn * dwork[*n + k] + cs * qtbpj;
		    dwork[*n + k] = temp;
		    i__3 = *n - k + 1;
		    drot_(&i__3, &r__[k + k * r_dim1], &c__1, &dwork[k], &
			    c__1, &cs, &sn);

		}
/* L40: */
	    }

	}

/*        Store the diagonal element of S and, if COND <> 'E', restore */
/*        the corresponding diagonal element of R. */

	dwork[j] = r__[j + j * r_dim1];
	if (! econd) {
	    r__[j + j * r_dim1] = x[j];
	}
/* L50: */
    }

/*     Solve the triangular system for z. If the system is singular, */
/*     then obtain a least squares solution. */

    if (econd) {
	toldef = *tol;
	if (toldef <= 0.) {

/*           Use the default tolerance in rank determination. */

	    toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	}

/*        Interchange the strict upper and lower triangular parts of R. */

	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    dswap_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], ldr);
/* L60: */
	}

/*        Estimate the reciprocal condition number of S and set the rank. */
/*        Additional workspace: 2*N. */

	i__1 = *ldwork - (*n << 1);
	mb03od_("No QR", n, n, &r__[r_offset], ldr, &ipvt[1], &toldef, &c_b19,
		 &dwork[1], rank, dum, &dwork[(*n << 1) + 1], &i__1, info, (
		ftnlen)5);
	r__[r_dim1 + 1] = x[1];

/*        Restore the strict upper and lower triangular parts of R. */

	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    dswap_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], ldr);
	    r__[j + j * r_dim1] = x[j];
/* L70: */
	}

    } else if (ncond) {

/*        Determine rank(S) by checking zero diagonal entries. */

	*rank = *n;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (dwork[j] == 0. && *rank == *n) {
		*rank = j - 1;
	    }
/* L80: */
	}

    }

    dum[0] = 0.;
    if (*rank < *n) {
	i__1 = *n - *rank;
	dcopy_(&i__1, dum, &c__0, &dwork[*n + *rank + 1], &c__1);
    }

/*     Solve S*z = c using back substitution. */

    for (j = *rank; j >= 1; --j) {
	temp = 0.;

	i__1 = *rank;
	for (i__ = j + 1; i__ <= i__1; ++i__) {
	    temp += r__[i__ + j * r_dim1] * dwork[*n + i__];
/* L90: */
	}

	dwork[*n + j] = (dwork[*n + j] - temp) / dwork[j];
/* L100: */
    }

/*     Permute the components of z back to components of x. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = dwork[*n + j];
/* L110: */
    }

    return 0;

/* *** Last line of MB02YD *** */
} /* mb02yd_ */

