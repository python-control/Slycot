/* TB04AD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb04ad_(char *rowcol, integer *n, integer *m, integer *p,
	 doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal 
	*c__, integer *ldc, doublereal *d__, integer *ldd, integer *nr, 
	integer *index, doublereal *dcoeff, integer *lddcoe, doublereal *
	ucoeff, integer *lduco1, integer *lduco2, doublereal *tol1, 
	doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen rowcol_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, dcoeff_dim1, dcoeff_offset, ucoeff_dim1, ucoeff_dim2, 
	    ucoeff_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, n1, ia;
    static char jobd[1];
    static integer itau;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), tb04ay_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer mplim;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxmp, jwork, mwork, pwork, kdcoef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lrococ, lrocor;
    static integer maxmpn;


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

/*     To find the transfer matrix T(s) of a given state-space */
/*     representation (A,B,C,D). T(s) is expressed as either row or */
/*     column polynomial vectors over monic least common denominator */
/*     polynomials. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ROWCOL  CHARACTER*1 */
/*             Indicates whether the transfer matrix T(s) is required */
/*             as rows or columns over common denominators as follows: */
/*             = 'R':  T(s) is required as rows over common denominators; */
/*             = 'C':  T(s) is required as columns over common */
/*                     denominators. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, i.e. the */
/*             order of the original state dynamics matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the upper block Hessenberg state dynamics matrix A of a */
/*             transformed representation for the original system: this */
/*             is completely controllable if ROWCOL = 'R', or completely */
/*             observable if ROWCOL = 'C'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M), */
/*             if ROWCOL = 'R', and (LDB,MAX(M,P)) if ROWCOL = 'C'. */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B; if */
/*             ROWCOL = 'C', the remainder of the leading N-by-MAX(M,P) */
/*             part is used as internal workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the transformed input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C; if */
/*             ROWCOL = 'C', the remainder of the leading MAX(M,P)-by-N */
/*             part is used as internal workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P)   if ROWCOL = 'R'; */
/*             LDC >= MAX(1,M,P) if ROWCOL = 'C'. */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M), */
/*             if ROWCOL = 'R', and (LDD,MAX(M,P)) if ROWCOL = 'C'. */
/*             The leading P-by-M part of this array must contain the */
/*             original direct transmission matrix D; if ROWCOL = 'C', */
/*             this array is modified internally, but restored on exit, */
/*             and the remainder of the leading MAX(M,P)-by-MAX(M,P) */
/*             part is used as internal workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P)   if ROWCOL = 'R'; */
/*             LDD >= MAX(1,M,P) if ROWCOL = 'C'. */

/*     NR      (output) INTEGER */
/*             The order of the transformed state-space representation. */

/*     INDEX   (output) INTEGER array, dimension (porm), where porm = P, */
/*             if ROWCOL = 'R', and porm = M, if ROWCOL = 'C'. */
/*             The degrees of the denominator polynomials. */

/*     DCOEFF  (output) DOUBLE PRECISION array, dimension (LDDCOE,N+1) */
/*             The leading porm-by-kdcoef part of this array contains */
/*             the coefficients of each denominator polynomial, where */
/*             kdcoef = MAX(INDEX(I)) + 1. */
/*             DCOEFF(I,K) is the coefficient in s**(INDEX(I)-K+1) of */
/*             the I-th denominator polynomial, where K = 1,2,...,kdcoef. */

/*     LDDCOE  INTEGER */
/*             The leading dimension of array DCOEFF. */
/*             LDDCOE >= MAX(1,P) if ROWCOL = 'R'; */
/*             LDDCOE >= MAX(1,M) if ROWCOL = 'C'. */

/*     UCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDUCO1,LDUCO2,N+1) */
/*             If ROWCOL = 'R' then porp = M, otherwise porp = P. */
/*             The leading porm-by-porp-by-kdcoef part of this array */
/*             contains the coefficients of the numerator matrix U(s). */
/*             UCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef; */
/*             if ROWCOL = 'R' then iorj = I, otherwise iorj = J. */
/*             Thus for ROWCOL = 'R', U(s) = */
/*             diag(s**INDEX(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...). */

/*     LDUCO1  INTEGER */
/*             The leading dimension of array UCOEFF. */
/*             LDUCO1 >= MAX(1,P) if ROWCOL = 'R'; */
/*             LDUCO1 >= MAX(1,M) if ROWCOL = 'C'. */

/*     LDUCO2  INTEGER */
/*             The second dimension of array UCOEFF. */
/*             LDUCO2 >= MAX(1,M) if ROWCOL = 'R'; */
/*             LDUCO2 >= MAX(1,P) if ROWCOL = 'C'. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             The tolerance to be used in determining the i-th row of */
/*             T(s), where i = 1,2,...,porm. If the user sets TOL1 > 0, */
/*             then the given value of TOL1 is used as an absolute */
/*             tolerance; elements with absolute value less than TOL1 are */
/*             considered neglijible. If the user sets TOL1 <= 0, then */
/*             an implicitly computed, default tolerance, defined in */
/*             the SLICOT Library routine TB01ZD, is used instead. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance to be used to separate out a controllable */
/*             subsystem of (A,B,C). If the user sets TOL2 > 0, then */
/*             the given value of TOL2 is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL2 is considered to be of full rank.  If the user sets */
/*             TOL2 <= 0, then an implicitly computed, default tolerance, */
/*             defined in the SLICOT Library routine TB01UD, is used */
/*             instead. */

/*     Workspace */

/*     IWORK   DOUBLE PRECISION array, dimension (N+MAX(M,P)) */
/*             On exit, if INFO = 0, the first nonzero elements of */
/*             IWORK(1:N) return the orders of the diagonal blocks of A. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*(N + 1) + MAX(N*MP + 2*N + MAX(N,MP), */
/*                                       3*MP, PM)), */
/*             where MP = M, PM = P, if ROWCOL = 'R'; */
/*                   MP = P, PM = M, if ROWCOL = 'C'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The method for transfer matrices factorized by rows will be */
/*     described here: T(s) factorized by columns is dealt with by */
/*     operating on the dual of the original system.  Each row of */
/*     T(s) is simply a single-output relatively left prime polynomial */
/*     matrix representation, so can be calculated by applying a */
/*     simplified version of the Orthogonal Structure Theorem to a */
/*     minimal state-space representation for the corresponding row of */
/*     the given system. A minimal state-space representation is obtained */
/*     using the Orthogonal Canonical Form to first separate out a */
/*     completely controllable one for the overall system and then, for */
/*     each row in turn, applying it again to the resulting dual SIMO */
/*     (single-input multi-output) system. Note that the elements of the */
/*     transformed matrix A so calculated are individually scaled in a */
/*     way which guarantees a monic denominator polynomial. */

/*     REFERENCES */

/*     [1] Williams, T.W.C. */
/*         An Orthogonal Structure Theorem for Linear Systems. */
/*         Control Systems Research Group, Kingston Polytechnic, */
/*         Internal Report 82/2, 1982. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998. */
/*     Supersedes Release 3.0 routine TB01QD. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Controllability, dual system, minimal realization, orthogonal */
/*     canonical form, orthogonal transformation, polynomial matrix, */
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
    --index;
    dcoeff_dim1 = *lddcoe;
    dcoeff_offset = 1 + dcoeff_dim1;
    dcoeff -= dcoeff_offset;
    ucoeff_dim1 = *lduco1;
    ucoeff_dim2 = *lduco2;
    ucoeff_offset = 1 + ucoeff_dim1 * (1 + ucoeff_dim2);
    ucoeff -= ucoeff_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    lrocor = lsame_(rowcol, "R", (ftnlen)1, (ftnlen)1);
    lrococ = lsame_(rowcol, "C", (ftnlen)1, (ftnlen)1);
    maxmp = max(*m,*p);
    mplim = max(1,maxmp);
    maxmpn = max(maxmp,*n);
    n1 = max(1,*n);
    if (lrocor) {

/*        T(s) given as rows over common denominators. */

	pwork = *p;
	mwork = *m;
    } else {

/*        T(s) given as columns over common denominators. */

	pwork = *m;
	mwork = *p;
    }

/*     Test the input scalar arguments. */

    if (! lrocor && ! lrococ) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*lda < n1) {
	*info = -6;
    } else if (*ldb < n1) {
	*info = -8;
    } else if (lrococ && *ldc < mplim || *ldc < max(1,*p)) {
	*info = -10;
    } else if (lrococ && *ldd < mplim || *ldd < max(1,*p)) {
	*info = -12;
    } else if (*lddcoe < max(1,pwork)) {
	*info = -16;
    } else if (*lduco1 < max(1,pwork)) {
	*info = -18;
    } else if (*lduco2 < max(1,mwork)) {
	*info = -19;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n * mwork + (*n << 1) + max(*n,mwork), i__4 = mwork * 3, i__3 
		= max(i__3,i__4);
	i__1 = 1, i__2 = *n * (*n + 1) + max(i__3,pwork);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -24;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB04AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (maxmpn == 0) {
	return 0;
    }

    *(unsigned char *)jobd = 'D';
    ia = 1;
    itau = ia + *n * *n;
    jwork = itau + *n;

    if (lrococ) {

/*        Initialization for T(s) given as columns over common */
/*        denominators. */

	ab07md_(jobd, n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)1);
    }

/*     Initialize polynomial matrix U(s) to zero. */

    i__1 = *n + 1;
    for (k = 1; k <= i__1; ++k) {
	dlaset_("Full", &pwork, &mwork, &c_b8, &c_b8, &ucoeff[(k * 
		ucoeff_dim2 + 1) * ucoeff_dim1 + 1], lduco1, (ftnlen)4);
/* L10: */
    }

/*     Calculate T(s) by applying the Orthogonal Structure Theorem to */
/*     each of the PWORK MISO subsystems (A,B,C:I,D:I) in turn. */

    i__1 = *ldwork - jwork + 1;
    tb04ay_(n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, nr, &index[1], &dcoeff[
	    dcoeff_offset], lddcoe, &ucoeff[ucoeff_offset], lduco1, lduco2, &
	    dwork[ia], &n1, &dwork[itau], tol1, tol2, &iwork[1], &dwork[jwork]
	    , &i__1, info);
    dwork[1] = dwork[jwork] + (doublereal) (jwork - 1);

    if (lrococ) {

/*        For T(s) factorized by columns, return to original (dual of */
/*        dual) system, and reorder the rows and columns to get an upper */
/*        block Hessenberg state dynamics matrix. */

	i__1 = iwork[1] + iwork[2] - 1;
	i__2 = *n - 1;
	tb01xd_(jobd, n, &mwork, &pwork, &i__1, &i__2, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, 
		info, (ftnlen)1);

	if (mplim != 1) {

/*           Also, transpose U(s) (not 1-by-1). */

	    kdcoef = 0;

	    i__1 = pwork;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		i__2 = kdcoef, i__3 = index[i__];
		kdcoef = max(i__2,i__3);
/* L20: */
	    }

	    ++kdcoef;

	    i__1 = kdcoef;
	    for (k = 1; k <= i__1; ++k) {

		i__2 = mplim - 1;
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

    return 0;
/* *** Last line of TB04AD *** */
} /* tb04ad_ */

