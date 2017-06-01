/* TB03AD.f -- translated by f2c (version 20100827).
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
static doublereal c_b20 = 0.;
static doublereal c_b53 = 1.;
static integer c_n1 = -1;

/* Subroutine */ int tb03ad_(char *leri, char *equil, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, integer 
	*nr, integer *index, doublereal *pcoeff, integer *ldpco1, integer *
	ldpco2, doublereal *qcoeff, integer *ldqco1, integer *ldqco2, 
	doublereal *vcoeff, integer *ldvco1, integer *ldvco2, doublereal *tol,
	 integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen leri_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, 
	    qcoeff_dim2, qcoeff_offset, vcoeff_dim1, vcoeff_dim2, 
	    vcoeff_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ic, iz, ioff, joff, ncol, kmax, itau, nrow;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), ma02gd_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *), tb01id_(char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), tc01od_(char *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), tb01ud_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb03ay_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), tb01yd_(integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static integer mplim, ncont, maxmp;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork, kwork, istop, kplus, mwork, pwork, indblk, irankc, 
	    kpcoef, nreflc;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical lleril;
    static integer ldwric;
    static logical llerir, lequil;
    static integer ifirst;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer istart, inplus, wrkopt;


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

/*     To find a relatively prime left polynomial matrix representation */
/*     inv(P(s))*Q(s) or right polynomial matrix representation */
/*     Q(s)*inv(P(s)) with the same transfer matrix T(s) as that of a */
/*     given state-space representation, i.e. */

/*        inv(P(s))*Q(s) = Q(s)*inv(P(s)) = T(s) = C*inv(s*I-A)*B + D. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether the left polynomial matrix */
/*             representation or the right polynomial matrix */
/*             representation is required as follows: */
/*             = 'L':  A left matrix fraction is required; */
/*             = 'R':  A right matrix fraction is required. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the triplet */
/*             (A,B,C), before computing a minimal state-space */
/*             representation, as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

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
/*             the upper block Hessenberg state dynamics matrix Amin of a */
/*             minimal realization for the original system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B; the remainder */
/*             of the leading N-by-MAX(M,P) part is used as internal */
/*             workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the transformed input/state matrix Bmin. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C; the remainder */
/*             of the leading MAX(M,P)-by-N part is used as internal */
/*             workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix Cmin. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,MAX(M,P)) */
/*             The leading P-by-M part of this array must contain the */
/*             original direct transmission matrix D; the remainder of */
/*             the leading MAX(M,P)-by-MAX(M,P) part is used as internal */
/*             workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M,P). */

/*     NR      (output) INTEGER */
/*             The order of the minimal state-space representation */
/*             (Amin,Bmin,Cmin). */

/*     INDEX   (output) INTEGER array, dimension (P), if LERI = 'L', or */
/*                                     dimension (M), if LERI = 'R'. */
/*             If LERI = 'L', INDEX(I), I = 1,2,...,P, contains the */
/*             maximum degree of the polynomials in the I-th row of the */
/*             denominator matrix P(s) of the left polynomial matrix */
/*             representation. */
/*             These elements are ordered so that */
/*             INDEX(1) >= INDEX(2) >= ... >= INDEX(P). */
/*             If LERI = 'R', INDEX(I), I = 1,2,...,M, contains the */
/*             maximum degree of the polynomials in the I-th column of */
/*             the denominator matrix P(s) of the right polynomial */
/*             matrix representation. */
/*             These elements are ordered so that */
/*             INDEX(1) >= INDEX(2) >= ... >= INDEX(M). */

/*     PCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,N+1) */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             The leading porm-by-porm-by-kpcoef part of this array */
/*             contains the coefficients of the denominator matrix P(s), */
/*             where kpcoef = MAX(INDEX(I)) + 1. */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if */
/*             LERI = 'L' then iorj = I, otherwise iorj = J. */
/*             Thus for LERI = 'L', P(s) = */
/*             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...). */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P), if LERI = 'L'; */
/*             LDPCO1 >= MAX(1,M), if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P), if LERI = 'L'; */
/*             LDPCO2 >= MAX(1,M), if LERI = 'R'. */

/*     QCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,N+1) */
/*             If LERI = 'L' then porp = M, otherwise porp = P. */
/*             If LERI = 'L', the leading porm-by-porp-by-kpcoef part */
/*             of this array contains the coefficients of the numerator */
/*             matrix Q(s). */
/*             If LERI = 'R', the leading porp-by-porm-by-kpcoef part */
/*             of this array contains the coefficients of the numerator */
/*             matrix Q(s). */
/*             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,P),   if LERI = 'L'; */
/*             LDQCO1 >= MAX(1,M,P), if LERI = 'R'. */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M),   if LERI = 'L'; */
/*             LDQCO2 >= MAX(1,M,P), if LERI = 'R'. */

/*     VCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDVCO1,LDVCO2,N+1) */
/*             The leading porm-by-NR-by-kpcoef part of this array */
/*             contains the coefficients of the intermediate matrix V(s). */
/*             VCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */

/*     LDVCO1  INTEGER */
/*             The leading dimension of array VCOEFF. */
/*             LDVCO1 >= MAX(1,P), if LERI = 'L'; */
/*             LDVCO1 >= MAX(1,M), if LERI = 'R'. */

/*     LDVCO2  INTEGER */
/*             The second dimension of array VCOEFF.  LDVCO2 >= MAX(1,N). */

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
/*             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P), PM*(PM + 2)) */
/*             where  PM = P, if LERI = 'L'; */
/*                    PM = M, if LERI = 'R'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if a singular matrix was encountered during the */
/*                   computation of V(s); */
/*             = 2:  if a singular matrix was encountered during the */
/*                   computation of P(s). */

/*     METHOD */

/*     The method for a left matrix fraction will be described here: */
/*     right matrix fractions are dealt with by constructing a left */
/*     fraction for the dual of the original system. The first step is to */
/*     obtain, by means of orthogonal similarity transformations, a */
/*     minimal state-space representation (Amin,Bmin,Cmin,D) for the */
/*     original system (A,B,C,D), where Amin is lower block Hessenberg */
/*     with all its superdiagonal blocks upper triangular and Cmin has */
/*     all but its first rank(C) columns zero.  The number and dimensions */
/*     of the blocks of Amin now immediately yield the row degrees of */
/*     P(s) with P(s) row proper: furthermore, the P-by-NR polynomial */
/*     matrix V(s) (playing a similar role to S(s) in Wolovich's */
/*     Structure Theorem) can be calculated a column block at a time, in */
/*     reverse order, from Amin. P(s) is then found as if it were the */
/*     O-th column block of V(s) (using Cmin as well as Amin), while */
/*     Q(s) = (V(s) * Bmin) + (P(s) * D). Finally, a special similarity */
/*     transformation is used to put Amin in an upper block Hessenberg */
/*     form. */

/*     REFERENCES */

/*     [1] Williams, T.W.C. */
/*         An Orthogonal Structure Theorem for Linear Systems. */
/*         Kingston Polytechnic Control Systems Research Group, */
/*         Internal Report 82/2, July 1982. */

/*     [2] Patel, R.V. */
/*         On Computing Matrix Fraction Descriptions and Canonical */
/*         Forms of Linear Time-Invariant Systems. */
/*         UMIST Control Systems Centre Report 489, 1980. */
/*         (Algorithms 1 and 2, extensively modified). */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998. */
/*     Supersedes Release 3.0 routine TB01SD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000. */

/*     KEYWORDS */

/*     Canonical form, coprime matrix fraction, dual system, elementary */
/*     polynomial operations, Hessenberg form, minimal realization, */
/*     orthogonal transformation, polynomial matrix, state-space */
/*     representation, transfer matrix. */

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
    pcoeff_dim1 = *ldpco1;
    pcoeff_dim2 = *ldpco2;
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
    pcoeff -= pcoeff_offset;
    qcoeff_dim1 = *ldqco1;
    qcoeff_dim2 = *ldqco2;
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
    qcoeff -= qcoeff_offset;
    vcoeff_dim1 = *ldvco1;
    vcoeff_dim2 = *ldvco2;
    vcoeff_offset = 1 + vcoeff_dim1 * (1 + vcoeff_dim2);
    vcoeff -= vcoeff_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    lleril = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
    llerir = lsame_(leri, "R", (ftnlen)1, (ftnlen)1);
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
    maxmp = max(*m,*p);
    mplim = max(1,maxmp);
    if (llerir) {

/*        Initialization for right matrix fraction. */

	pwork = *m;
	mwork = *p;
    } else {

/*        Initialization for left matrix fraction. */

	pwork = *p;
	mwork = *m;
    }

/*     Test the input scalar arguments. */

    if (! lleril && ! llerir) {
	*info = -1;
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < mplim) {
	*info = -11;
    } else if (*ldd < mplim) {
	*info = -13;
    } else if (*ldpco1 < max(1,pwork)) {
	*info = -17;
    } else if (*ldpco2 < max(1,pwork)) {
	*info = -18;
    } else if (*ldqco1 < max(1,pwork) || llerir && *ldqco1 < mplim) {
	*info = -20;
    } else if (*ldqco2 < max(1,mwork) || llerir && *ldqco2 < mplim) {
	*info = -21;
    } else if (*ldvco1 < max(1,pwork)) {
	*info = -23;
    } else if (*ldvco2 < max(1,*n)) {
	*info = -24;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n, i__4 = maxmp * 3;
	i__1 = 1, i__2 = *n + max(i__3,i__4), i__1 = max(i__1,i__2), i__2 = 
		pwork * (pwork + 2);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -28;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB03AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*n,*m);
    if (max(i__1,*p) == 0) {
	*nr = 0;
	dwork[1] = 1.;
	return 0;
    }

    if (llerir) {

/*        For right matrix fraction, obtain dual system. */

	ab07md_("D", n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)1);
    }

/*     Obtain minimal realization, in canonical form, for this system. */
/*     Part of the code in SLICOT routine TB01PD is included in-line */
/*     here. (TB01PD cannot be directly used.) */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/*     If required, balance the triplet (A,B,C) (default MAXRED). */
/*     Workspace: need N. */

    if (lequil) {
	maxred = 0.;
	tb01id_("A", n, &mwork, &pwork, &maxred, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &dwork[1], info, (ftnlen)
		1);
    }

    iz = 1;
    itau = 1;
    jwork = itau + *n;

/*     Separate out controllable subsystem (of order NCONT): */
/*     A <-- Z'*A*Z,  B <-- Z'*B,  C <-- C*Z. */

/*     Workspace: need   N + MAX(N, 3*MWORK, PWORK). */
/*                prefer larger. */

    i__1 = *ldwork - jwork + 1;
    tb01ud_("No Z", n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &ncont, &indblk, &iwork[1], &dwork[iz], &c__1,
	     &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, info, (
	    ftnlen)4);

    wrkopt = (integer) dwork[jwork] + jwork - 1;

/*     Separate out the observable subsystem (of order NR): */
/*     Form the dual of the subsystem of order NCONT (which is */
/*     controllable), leaving rest as it is. */

    ab07md_("Z", &ncont, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb,
	     &c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*     And separate out the controllable part of this dual subsystem. */

/*     Workspace: need   NCONT + MAX(NCONT, 3*PWORK, MWORK). */
/*                prefer larger. */

    i__1 = *ldwork - jwork + 1;
    tb01ud_("No Z", &ncont, &pwork, &mwork, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, nr, &indblk, &iwork[1], &dwork[iz], &
	    c__1, &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, 
	    info, (ftnlen)4);

/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Retranspose, giving controllable and observable (i.e. minimal) */
/*     part of original system. */

    ab07md_("Z", nr, &pwork, &mwork, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*     Annihilate the trailing components of IWORK(1:N). */

    i__1 = *n;
    for (i__ = indblk + 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
/* L10: */
    }

/*     Initialize polynomial matrices P(s), Q(s) and V(s) to zero. */

    i__1 = *n + 1;
    for (k = 1; k <= i__1; ++k) {
	dlaset_("Full", &pwork, &pwork, &c_b20, &c_b20, &pcoeff[(k * 
		pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (ftnlen)4);
	dlaset_("Full", &pwork, &mwork, &c_b20, &c_b20, &qcoeff[(k * 
		qcoeff_dim2 + 1) * qcoeff_dim1 + 1], ldqco1, (ftnlen)4);
	dlaset_("Full", &pwork, nr, &c_b20, &c_b20, &vcoeff[(k * vcoeff_dim2 
		+ 1) * vcoeff_dim1 + 1], ldvco1, (ftnlen)4);
/* L20: */
    }

/*     Finish initializing V(s), and set up row degrees of P(s). */

    inplus = indblk + 1;
    istart = 1;
    joff = *nr;

    i__1 = indblk;
    for (k = 1; k <= i__1; ++k) {
	kwork = inplus - k;
	kplus = kwork + 1;
	istop = iwork[kwork];
	joff -= istop;

	i__2 = istop;
	for (i__ = istart; i__ <= i__2; ++i__) {
	    index[i__] = kwork;
	    vcoeff[i__ + (joff + i__ + kplus * vcoeff_dim2) * vcoeff_dim1] = 
		    1.;
/* L30: */
	}

	istart = istop + 1;
/* L40: */
    }

/*     ISTART = IWORK(1)+1 now: if .LE. PWORK, set up final rows of P(s). */

    i__1 = pwork;
    for (i__ = istart; i__ <= i__1; ++i__) {
	index[i__] = 0;
	pcoeff[i__ + (i__ + pcoeff_dim2) * pcoeff_dim1] = 1.;
/* L50: */
    }

/*     Triangularize the superdiagonal blocks of Amin. */

    nrow = iwork[indblk];
    ioff = *nr - nrow;
    kmax = indblk - 1;
    itau = 1;
    ifirst = 0;
    if (indblk > 2) {
	ifirst = ioff - iwork[kmax];
    }

/*     QR decomposition of each superdiagonal block of A in turn */
/*     (done in reverse order to preserve upper triangular blocks in A). */

    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {

/*        Calculate dimensions of new block & its position in A. */

	kwork = indblk - k;
	ncol = nrow;
	nrow = iwork[kwork];
	joff = ioff;
	ioff -= nrow;
	nreflc = min(nrow,ncol);
	jwork = itau + nreflc;
	if (kwork >= 2) {
	    ifirst -= iwork[kwork - 1];
	}

/*        Find QR decomposition of this (full rank) block: */
/*        block = QR.  No pivoting is needed. */

/*        Workspace: need   MIN(NROW,NCOL) + NCOL; */
/*                   prefer MIN(NROW,NCOL) + NCOL*NB. */

	i__2 = *ldwork - jwork + 1;
	dgeqrf_(&nrow, &ncol, &a[ioff + 1 + (joff + 1) * a_dim1], lda, &dwork[
		itau], &dwork[jwork], &i__2, info);

/* Computing MAX */
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__2,i__3);

/*        Premultiply appropriate row block of A by Q'. */

/*        Workspace: need   MIN(NROW,NCOL) + JOFF; */
/*                   prefer MIN(NROW,NCOL) + JOFF*NB. */

	i__2 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &nrow, &joff, &nreflc, &a[ioff + 1 + (
		joff + 1) * a_dim1], lda, &dwork[itau], &a[ioff + 1 + a_dim1],
		 lda, &dwork[jwork], &i__2, info, (ftnlen)4, (ftnlen)9);

/* Computing MAX */
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__2,i__3);

/*        Premultiply appropriate row block of B by Q' also. */

/*        Workspace: need   MIN(NROW,NCOL) + MWORK; */
/*                   prefer MIN(NROW,NCOL) + MWORK*NB. */

	i__2 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &nrow, &mwork, &nreflc, &a[ioff + 1 + (
		joff + 1) * a_dim1], lda, &dwork[itau], &b[ioff + 1 + b_dim1],
		 ldb, &dwork[jwork], &i__2, info, (ftnlen)4, (ftnlen)9);

/* Computing MAX */
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__2,i__3);

/*        And postmultiply the non-zero part of appropriate column */
/*        block of A by Q. */

/*        Workspace: need   MIN(NROW,NCOL) + NR; */
/*                   prefer MIN(NROW,NCOL) + NR*NB. */

	i__2 = *nr - ifirst;
	i__3 = *ldwork - jwork + 1;
	dormqr_("Right", "No Transpose", &i__2, &nrow, &nreflc, &a[ioff + 1 + 
		(joff + 1) * a_dim1], lda, &dwork[itau], &a[ifirst + 1 + (
		ioff + 1) * a_dim1], lda, &dwork[jwork], &i__3, info, (ftnlen)
		5, (ftnlen)12);

/* Computing MAX */
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__2,i__3);

/*        Annihilate the lower triangular part of the block in A. */

	if (k != kmax && nrow > 1) {
	    i__2 = nrow - 1;
	    dlaset_("Lower", &i__2, &ncol, &c_b20, &c_b20, &a[ioff + 2 + (
		    joff + 1) * a_dim1], lda, (ftnlen)5);
	}

/* L60: */
    }

/*     Finally: postmultiply non-zero columns of C by Q (K = KMAX). */

/*     Workspace: need   MIN(NROW,NCOL) + PWORK; */
/*                prefer MIN(NROW,NCOL) + PWORK*NB. */

    i__1 = *ldwork - jwork + 1;
    dormqr_("Right", "No Transpose", &pwork, &nrow, &nreflc, &a[ioff + 1 + (
	    joff + 1) * a_dim1], lda, &dwork[itau], &c__[c_offset], ldc, &
	    dwork[jwork], &i__1, info, (ftnlen)5, (ftnlen)12);

/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Annihilate the lower triangular part of the block in A. */

    if (nrow > 1) {
	i__1 = nrow - 1;
	dlaset_("Lower", &i__1, &ncol, &c_b20, &c_b20, &a[ioff + 2 + (joff + 
		1) * a_dim1], lda, (ftnlen)5);
    }

/*     Calculate the (PWORK x NR) polynomial matrix V(s) ... */

    tb03ay_(nr, &a[a_offset], lda, &indblk, &iwork[1], &vcoeff[vcoeff_offset],
	     ldvco1, ldvco2, &pcoeff[pcoeff_offset], ldpco1, ldpco2, info);

    if (*info != 0) {
	*info = 1;
	return 0;
    } else {

/*        And then use this matrix to calculate P(s): first store */
/*        C1 from C. */

	ic = 1;
	irankc = iwork[1];
	ldwric = max(1,pwork);
	dlacpy_("Full", &pwork, &irankc, &c__[c_offset], ldc, &dwork[ic], &
		ldwric, (ftnlen)4);

	if (irankc < pwork) {

/*           rank(C) .LT. PWORK: obtain QR decomposition of C1, */
/*           giving R and Q. */

/*           Workspace: need   PWORK*IRANKC + 2*IRANKC; */
/*                      prefer PWORK*IRANKC +   IRANKC + IRANKC*NB. */

	    itau = ic + ldwric * irankc;
	    jwork = itau + irankc;

	    i__1 = *ldwork - jwork + 1;
	    dgeqrf_(&pwork, &irankc, &dwork[ic], &ldwric, &dwork[itau], &
		    dwork[jwork], &i__1, info);

/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__1,i__2);

/*           First IRANKC rows of Pbar(s) are given by Wbar(s) * inv(R). */
/*           Check for zero diagonal elements of R. */

	    i__1 = irankc;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (dwork[ic + (i__ - 1) * ldwric + i__ - 1] == 0.) {

/*                 Error return. */

		    *info = 2;
		    return 0;
		}
/* L70: */
	    }

	    nrow = irankc;

	    i__1 = inplus;
	    for (k = 1; k <= i__1; ++k) {
		dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &nrow, &
			irankc, &c_b53, &dwork[ic], &ldwric, &pcoeff[(k * 
			pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (ftnlen)
			5, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		nrow = iwork[k];
/* L80: */
	    }

/*           P(s) itself is now given by Pbar(s) * Q'. */

	    nrow = pwork;

	    i__1 = inplus;
	    for (k = 1; k <= i__1; ++k) {

/*              Workspace: need   PWORK*IRANKC + IRANKC + NROW; */
/*                         prefer PWORK*IRANKC + IRANKC + NROW*NB. */

		i__2 = *ldwork - jwork + 1;
		dormqr_("Right", "Transpose", &nrow, &pwork, &irankc, &dwork[
			ic], &ldwric, &dwork[itau], &pcoeff[(k * pcoeff_dim2 
			+ 1) * pcoeff_dim1 + 1], ldpco1, &dwork[jwork], &i__2,
			 info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
		i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
		wrkopt = max(i__2,i__3);
		nrow = iwork[k];
/* L90: */
	    }

	} else {

/*           Special case rank(C) = PWORK, full: */
/*           no QR decomposition (P(s)=Wbar(s)*inv(C1)). */

	    dgetrf_(&pwork, &pwork, &dwork[ic], &ldwric, &iwork[*n + 1], info)
		    ;

	    if (*info != 0) {

/*              Error return. */

		*info = 2;
		return 0;
	    } else {

		nrow = irankc;

/*              Workspace: need   PWORK*IRANKC + N. */

		i__1 = inplus;
		for (k = 1; k <= i__1; ++k) {
		    dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &
			    nrow, &pwork, &c_b53, &dwork[ic], &ldwric, &
			    pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], 
			    ldpco1, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)
			    8);
		    dtrsm_("Right", "Lower", "No Transpose", "Unit", &nrow, &
			    pwork, &c_b53, &dwork[ic], &ldwric, &pcoeff[(k * 
			    pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (
			    ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		    ma02gd_(&nrow, &pcoeff[(k * pcoeff_dim2 + 1) * 
			    pcoeff_dim1 + 1], ldpco1, &c__1, &pwork, &iwork[*
			    n + 1], &c_n1);
		    nrow = iwork[k];
/* L100: */
		}
	    }
	}

/*        Finally, Q(s) = V(s) * B + P(s) * D can now be evaluated. */

	nrow = pwork;

	i__1 = inplus;
	for (k = 1; k <= i__1; ++k) {
	    dgemm_("No transpose", "No transpose", &nrow, &mwork, nr, &c_b53, 
		    &vcoeff[(k * vcoeff_dim2 + 1) * vcoeff_dim1 + 1], ldvco1, 
		    &b[b_offset], ldb, &c_b20, &qcoeff[(k * qcoeff_dim2 + 1) *
		     qcoeff_dim1 + 1], ldqco1, (ftnlen)12, (ftnlen)12);
	    dgemm_("No transpose", "No transpose", &nrow, &mwork, &pwork, &
		    c_b53, &pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], 
		    ldpco1, &d__[d_offset], ldd, &c_b53, &qcoeff[(k * 
		    qcoeff_dim2 + 1) * qcoeff_dim1 + 1], ldqco1, (ftnlen)12, (
		    ftnlen)12);
	    nrow = iwork[k];
/* L110: */
	}

    }

    if (llerir) {

/*        For right matrix fraction, return to original (dual of dual) */
/*        system. */

	ab07md_("Z", nr, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*        Also, obtain the dual of the polynomial matrix representation. */

	kpcoef = 0;

	i__1 = pwork;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    i__2 = kpcoef, i__3 = index[i__];
	    kpcoef = max(i__2,i__3);
/* L120: */
	}

	++kpcoef;
	tc01od_("L", &mwork, &pwork, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (ftnlen)
		1);
    } else {

/*        Reorder the rows and columns of the system, to get an upper */
/*        block Hessenberg matrix A of the minimal system. */

	tb01yd_(nr, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset]
		, ldc, info);
    }

/*     Set optimal workspace dimension. */

    dwork[1] = (doublereal) wrkopt;
    return 0;
/* *** Last line of TB03AD *** */
} /* tb03ad_ */

