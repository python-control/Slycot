/* TD04AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b8 = 0.;
static integer c__1 = 1;

/* Subroutine */ int td04ad_(char *rowcol, integer *m, integer *p, integer *
	index, doublereal *dcoeff, integer *lddcoe, doublereal *ucoeff, 
	integer *lduco1, integer *lduco2, integer *nr, doublereal *a, integer 
	*lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen rowcol_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, dcoeff_dim1, dcoeff_offset, ucoeff_dim1, ucoeff_dim2, 
	    ucoeff_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, n;
    extern /* Subroutine */ int tb01pd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), td03ay_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    static integer mplim;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jstop, mwork, pwork, kdcoef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lrococ, lrocor;


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

/*     To find a minimal state-space representation (A,B,C,D) for a */
/*     proper transfer matrix T(s) given as either row or column */
/*     polynomial vectors over denominator polynomials, possibly with */
/*     uncancelled common terms. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ROWCOL  CHARACTER*1 */
/*             Indicates whether the transfer matrix T(s) is given as */
/*             rows or columns over common denominators as follows: */
/*             = 'R':  T(s) is given as rows over common denominators; */
/*             = 'C':  T(s) is given as columns over common denominators. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     INDEX   (input) INTEGER array, dimension (porm), where porm = P, */
/*             if ROWCOL = 'R', and porm = M, if ROWCOL = 'C'. */
/*             This array must contain the degrees of the denominator */
/*             polynomials in D(s). */

/*     DCOEFF  (input) DOUBLE PRECISION array, dimension (LDDCOE,kdcoef), */
/*             where kdcoef = MAX(INDEX(I)) + 1. */
/*             The leading porm-by-kdcoef part of this array must contain */
/*             the coefficients of each denominator polynomial. */
/*             DCOEFF(I,K) is the coefficient in s**(INDEX(I)-K+1) of the */
/*             I-th denominator polynomial in D(s), where */
/*             K = 1,2,...,kdcoef. */

/*     LDDCOE  INTEGER */
/*             The leading dimension of array DCOEFF. */
/*             LDDCOE >= MAX(1,P) if ROWCOL = 'R'; */
/*             LDDCOE >= MAX(1,M) if ROWCOL = 'C'. */

/*     UCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDUCO1,LDUCO2,kdcoef) */
/*             The leading P-by-M-by-kdcoef part of this array must */
/*             contain the numerator matrix U(s); if ROWCOL = 'C', this */
/*             array is modified internally but restored on exit, and the */
/*             remainder of the leading MAX(M,P)-by-MAX(M,P)-by-kdcoef */
/*             part is used as internal workspace. */
/*             UCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef; */
/*             if ROWCOL = 'R' then iorj = I, otherwise iorj = J. */
/*             Thus for ROWCOL = 'R', U(s) = */
/*             diag(s**INDEX(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...). */

/*     LDUCO1  INTEGER */
/*             The leading dimension of array UCOEFF. */
/*             LDUCO1 >= MAX(1,P)   if ROWCOL = 'R'; */
/*             LDUCO1 >= MAX(1,M,P) if ROWCOL = 'C'. */

/*     LDUCO2  INTEGER */
/*             The second dimension of array UCOEFF. */
/*             LDUCO2 >= MAX(1,M)   if ROWCOL = 'R'; */
/*             LDUCO2 >= MAX(1,M,P) if ROWCOL = 'C'. */

/*     NR      (output) INTEGER */
/*             The order of the resulting minimal realization, i.e. the */
/*             order of the state dynamics matrix A. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N), */
/*                       porm */
/*             where N = SUM INDEX(I). */
/*                       I=1 */
/*             The leading NR-by-NR part of this array contains the upper */
/*             block Hessenberg state dynamics matrix A of a minimal */
/*             realization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P)) */
/*             The leading NR-by-M part of this array contains the */
/*             input/state matrix B of a minimal realization; the */
/*             remainder of the leading N-by-MAX(M,P) part is used as */
/*             internal workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-NR part of this array contains the */
/*             state/output matrix C of a minimal realization; the */
/*             remainder of the leading MAX(M,P)-by-N part is used as */
/*             internal workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M,P). */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M), */
/*             if ROWCOL = 'R', and (LDD,MAX(M,P)) if ROWCOL = 'C'. */
/*             The leading P-by-M part of this array contains the direct */
/*             transmission matrix D; if ROWCOL = 'C', the remainder of */
/*             the leading MAX(M,P)-by-MAX(M,P) part is used as internal */
/*             workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P)   if ROWCOL = 'R'; */
/*             LDD >= MAX(1,M,P) if ROWCOL = 'C'. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B, C). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance */
/*             (determined by the SLICOT routine TB01UD) is used instead. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+MAX(M,P)) */
/*             On exit, if INFO = 0, the first nonzero elements of */
/*             IWORK(1:N) return the orders of the diagonal blocks of A. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then i is the first integer for which */
/*                   ABS( DCOEFF(I,1) ) is so small that the calculations */
/*                   would overflow (see SLICOT Library routine TD03AY); */
/*                   that is, the leading coefficient of a polynomial is */
/*                   nearly zero; no state-space representation is */
/*                   calculated. */

/*     METHOD */

/*     The method for transfer matrices factorized by rows will be */
/*     described here: T(s) factorized by columns is dealt with by */
/*     operating on the dual T'(s). This description for T(s) is */
/*     actually the left polynomial matrix representation */

/*          T(s) = inv(D(s))*U(s), */

/*     where D(s) is diagonal with its (I,I)-th polynomial element of */
/*     degree INDEX(I). The first step is to check whether the leading */
/*     coefficient of any polynomial element of D(s) is approximately */
/*     zero; if so the routine returns with INFO > 0. Otherwise, */
/*     Wolovich's Observable Structure Theorem is used to construct a */
/*     state-space representation in observable companion form which */
/*     is equivalent to the above polynomial matrix representation. */
/*     The method is particularly easy here due to the diagonal form */
/*     of D(s). This state-space representation is not necessarily */
/*     controllable (as D(s) and U(s) are not necessarily relatively */
/*     left prime), but it is in theory completely observable; however, */
/*     its observability matrix may be poorly conditioned, so it is */
/*     treated as a general state-space representation and SLICOT */
/*     Library routine TB01PD is then called to separate out a minimal */
/*     realization from this general state-space representation by means */
/*     of orthogonal similarity transformations. */

/*     REFERENCES */

/*     [1] Patel, R.V. */
/*         Computation of Minimal-Order State-Space Realizations and */
/*         Observability Indices using Orthogonal Transformations. */
/*         Int. J. Control, 33, pp. 227-246, 1981. */

/*     [2] Wolovich, W.A. */
/*         Linear Multivariable Systems, (Theorem 4.3.3). */
/*         Springer-Verlag, 1974. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998. */
/*     Supersedes Release 3.0 routine TD01OD. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Controllability, elementary polynomial operations, minimal */
/*     realization, polynomial matrix, state-space representation, */
/*     transfer matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --index;
    dcoeff_dim1 = *lddcoe;
    dcoeff_offset = 1 + dcoeff_dim1;
    dcoeff -= dcoeff_offset;
    ucoeff_dim1 = *lduco1;
    ucoeff_dim2 = *lduco2;
    ucoeff_offset = 1 + ucoeff_dim1 * (1 + ucoeff_dim2);
    ucoeff -= ucoeff_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    lrocor = lsame_(rowcol, "R", (ftnlen)1, (ftnlen)1);
    lrococ = lsame_(rowcol, "C", (ftnlen)1, (ftnlen)1);
/* Computing MAX */
    i__1 = max(1,*m);
    mplim = max(i__1,*p);

/*     Test the input scalar arguments. */

    if (! lrocor && ! lrococ) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (lrocor && *lddcoe < max(1,*p) || lrococ && *lddcoe < max(1,*m))
	     {
	*info = -6;
    } else if (lrocor && *lduco1 < max(1,*p) || lrococ && *lduco1 < mplim) {
	*info = -8;
    } else if (lrocor && *lduco2 < max(1,*m) || lrococ && *lduco2 < mplim) {
	*info = -9;
    }

    n = 0;
    if (*info == 0) {
	if (lrocor) {

/*           Initialization for T(s) given as rows over common */
/*           denominators. */

	    pwork = *p;
	    mwork = *m;
	} else {

/*           Initialization for T(s) given as columns over common */
/*           denominators. */

	    pwork = *m;
	    mwork = *p;
	}

/*        Calculate N, the order of the resulting state-space */
/*        representation. */

	kdcoef = 0;

	i__1 = pwork;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    i__2 = kdcoef, i__3 = index[i__];
	    kdcoef = max(i__2,i__3);
	    n += index[i__];
/* L10: */
	}

	++kdcoef;

	if (*lda < max(1,n)) {
	    *info = -12;
	} else if (*ldb < max(1,n)) {
	    *info = -14;
	} else if (*ldc < mplim) {
	    *info = -16;
	} else if (lrocor && *ldd < max(1,*p) || lrococ && *ldd < mplim) {
	    *info = -18;
	} else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	    i__3 = n, i__4 = *m * 3, i__3 = max(i__3,i__4), i__4 = *p * 3;
	    i__1 = 1, i__2 = n + max(i__3,i__4);
	    if (*ldwork < max(i__1,i__2)) {
		*info = -22;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TD04AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(n,*m);
    if (max(i__1,*p) == 0) {
	*nr = 0;
	dwork[1] = 1.;
	return 0;
    }

    if (lrococ) {

/*        Initialize the remainder of the leading */
/*        MPLIM-by-MPLIM-by-KDCOEF part of U(s) to zero. */

	if (*p < *m) {

	    i__1 = kdcoef;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *m - *p;
		dlaset_("Full", &i__2, &mplim, &c_b8, &c_b8, &ucoeff[*p + 1 + 
			(k * ucoeff_dim2 + 1) * ucoeff_dim1], lduco1, (ftnlen)
			4);
/* L20: */
	    }

	} else if (*p > *m) {

	    i__1 = kdcoef;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *p - *m;
		dlaset_("Full", &mplim, &i__2, &c_b8, &c_b8, &ucoeff[(*m + 1 
			+ k * ucoeff_dim2) * ucoeff_dim1 + 1], lduco1, (
			ftnlen)4);
/* L30: */
	    }

	}

	if (mplim != 1) {

/*           Non-scalar T(s) factorized by columns: transpose it (i.e. */
/*           U(s)). */

	    jstop = mplim - 1;

	    i__1 = kdcoef;
	    for (k = 1; k <= i__1; ++k) {

		i__2 = jstop;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = mplim - j;
		    dswap_(&i__3, &ucoeff[j + 1 + (j + k * ucoeff_dim2) * 
			    ucoeff_dim1], &c__1, &ucoeff[j + (j + 1 + k * 
			    ucoeff_dim2) * ucoeff_dim1], lduco1);
/* L40: */
		}

/* L50: */
	    }

	}
    }

/*     Construct non-minimal state-space representation (by Wolovich's */
/*     Structure Theorem) which has transfer matrix T(s) or T'(s) as */
/*     appropriate ... */

    td03ay_(&mwork, &pwork, &index[1], &dcoeff[dcoeff_offset], lddcoe, &
	    ucoeff[ucoeff_offset], lduco1, lduco2, &n, &a[a_offset], lda, &b[
	    b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, info);
    if (*info > 0) {
	return 0;
    }

/*     and then separate out a minimal realization from this. */

/*     Workspace: need  N + MAX(N, 3*MWORK, 3*PWORK). */

    tb01pd_("Minimal", "Scale", &n, &mwork, &pwork, &a[a_offset], lda, &b[
	    b_offset], ldb, &c__[c_offset], ldc, nr, tol, &iwork[1], &dwork[1]
	    , ldwork, info, (ftnlen)7, (ftnlen)5);

    if (lrococ) {

/*        If T(s) originally factorized by columns, find dual of minimal */
/*        state-space representation, and reorder the rows and columns */
/*        to get an upper block Hessenberg state dynamics matrix. */

	k = iwork[1] + iwork[2] - 1;
	i__1 = *nr - 1;
	tb01xd_("D", nr, &mwork, &pwork, &k, &i__1, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, 
		info, (ftnlen)1);
	if (mplim != 1) {

/*           Also, retranspose U(s) if this is non-scalar. */

	    i__1 = kdcoef;
	    for (k = 1; k <= i__1; ++k) {

		i__2 = jstop;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = mplim - j;
		    dswap_(&i__3, &ucoeff[j + 1 + (j + k * ucoeff_dim2) * 
			    ucoeff_dim1], &c__1, &ucoeff[j + (j + 1 + k * 
			    ucoeff_dim2) * ucoeff_dim1], lduco1);
/* L60: */
		}

/* L70: */
	    }

	}
    }

    return 0;
/* *** Last line of TD04AD *** */
} /* td04ad_ */

