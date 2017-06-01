/* AG8BYZ.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static logical c_true = TRUE_;
static integer c__2 = 2;

/* Subroutine */ int ag8byz_(logical *first, integer *n, integer *m, integer *
	p, doublereal *svlmax, doublecomplex *abcd, integer *ldabcd, 
	doublecomplex *e, integer *lde, integer *nr, integer *pr, integer *
	ninfz, integer *dinfz, integer *nkronl, integer *infz, integer *kronl,
	 doublereal *tol, integer *iwork, doublereal *dwork, doublecomplex *
	zwork, integer *lzwork, integer *info)
{
    /* System generated locals */
    integer abcd_dim1, abcd_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublecomplex s;
    static doublereal t;
    static doublecomplex c1, c2;
    static integer n1;
    static doublecomplex s1, s2;
    static integer nb;
    static doublecomplex tc;
    static integer mn, pn, ro;
    static doublereal tt;
    static integer mn1, mp1, ro1, irc;
    static doublecomplex dum[1];
    static integer mpm, mui, mnr, icol, rank, itau, taui;
    static doublereal sval[3], smin, smax;
    static integer irow;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer muim1, sigma;
    static doublereal rcond;
    static integer ilast, jlast, ismin, ismax, mntau;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zlaic1_(integer *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static integer jwork1, jwork2;
    extern /* Subroutine */ int mb3oyz_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zlaset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zlartg_(doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *), 
	    zlapmt_(logical *, integer *, integer *, doublecomplex *, integer 
	    *, integer *);
    static doublereal sminpr, smaxpr;
    static logical lquery;
    extern /* Subroutine */ int zlatzm_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static integer wrkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

/*     To extract from the (N+P)-by-(M+N) descriptor system pencil */

/*        S(lambda) = ( B   A - lambda*E  ) */
/*                    ( D        C        ) */

/*     with E nonsingular and upper triangular a */
/*     (NR+PR)-by-(M+NR) "reduced" descriptor system pencil */

/*                           ( Br  Ar-lambda*Er ) */
/*              Sr(lambda) = (                  ) */
/*                           ( Dr     Cr        ) */

/*     having the same finite Smith zeros as the pencil */
/*     S(lambda) but with Dr, a PR-by-M full row rank */
/*     left upper trapezoidal matrix, and Er, an NR-by-NR */
/*     upper triangular nonsingular matrix. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FIRST   LOGICAL */
/*             Specifies if AG8BYZ is called first time or it is called */
/*             for an already reduced system, with D full column rank */
/*             with the last M rows in upper triangular form: */
/*             FIRST = .TRUE.,  first time called; */
/*             FIRST = .FALSE., not first time called. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of matrix B, the number of columns of */
/*             matrix C and the order of square matrices A and E. */
/*             N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrices B and D.  M >= 0. */
/*             M <= P if FIRST = .FALSE. . */

/*     P       (input) INTEGER */
/*             The number of rows of matrices C and D.  P >= 0. */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             During each reduction step, the rank-revealing QR */
/*             factorization of a matrix stops when the estimated minimum */
/*             singular value is smaller than TOL * MAX(SVLMAX,EMSV), */
/*             where EMSV is the estimated maximum singular value. */
/*             SVLMAX >= 0. */

/*     ABCD    (input/output) COMPLEX*16 array, dimension (LDABCD,M+N) */
/*             On entry, the leading (N+P)-by-(M+N) part of this array */
/*             must contain the compound matrix */
/*                      (  B   A  ) , */
/*                      (  D   C  ) */
/*             where A is an N-by-N matrix, B is an N-by-M matrix, */
/*             C is a P-by-N matrix and D is a P-by-M matrix. */
/*             If FIRST = .FALSE., then D must be a full column */
/*             rank matrix with the last M rows in upper triangular form. */
/*             On exit, the leading (NR+PR)-by-(M+NR) part of ABCD */
/*             contains the reduced compound matrix */
/*                       (  Br  Ar ) , */
/*                       (  Dr  Cr ) */
/*             where Ar is an NR-by-NR matrix, Br is an NR-by-M matrix, */
/*             Cr is a PR-by-NR matrix, Dr is a PR-by-M full row rank */
/*             left upper trapezoidal matrix with the first PR columns */
/*             in upper triangular form. */

/*     LDABCD  INTEGER */
/*             The leading dimension of array ABCD. */
/*             LDABCD >= MAX(1,N+P). */

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular nonsingular matrix E. */
/*             On exit, the leading NR-by-NR part contains the reduced */
/*             upper triangular nonsingular matrix Er. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     NR      (output) INTEGER */
/*             The order of the reduced matrices Ar and Er; also the */
/*             number of rows of the reduced matrix Br and the number */
/*             of columns of the reduced matrix Cr. */
/*             If Dr is invertible, NR is also the number of finite */
/*             Smith zeros. */

/*     PR      (output) INTEGER */
/*             The rank of the resulting matrix Dr; also the number of */
/*             rows of reduced matrices Cr and Dr. */

/*     NINFZ   (output) INTEGER */
/*             Number of infinite zeros.  NINFZ = 0 if FIRST = .FALSE. . */

/*     DINFZ   (output) INTEGER */
/*             The maximal multiplicity of infinite zeros. */
/*             DINFZ = 0 if FIRST = .FALSE. . */

/*     NKRONL  (output) INTEGER */
/*             The maximal dimension of left elementary Kronecker blocks. */

/*     INFZ    (output) INTEGER array, dimension (N) */
/*             INFZ(i) contains the number of infinite zeros of */
/*             degree i, where i = 1,2,...,DINFZ. */
/*             INFZ is not referenced if FIRST = .FALSE. . */

/*     KRONL   (output) INTEGER array, dimension (N+1) */
/*             KRONL(i) contains the number of left elementary Kronecker */
/*             blocks of dimension i-by-(i-1), where i = 1,2,...,NKRONL. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             If the user sets TOL <= 0, then an implicitly computed, */
/*             default tolerance TOLDEF = (N+P)*(N+M)*EPS,  is used */
/*             instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH). */
/*             NOTE that when SVLMAX > 0, the estimated ranks could be */
/*             less than those defined above (see SVLMAX).  TOL <= 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */
/*             If FIRST = .FALSE., IWORK is not referenced. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             LDWORK >= 2*MAX(M,P), if FIRST = .TRUE.; */
/*             LDWORK >= 2*P,        if FIRST = .FALSE. . */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= 1, if P = 0; otherwise */
/*             LZWORK >= MAX( 1, N+M-1, MIN(P,M) + MAX(3*M-1,N), 3*P ), */
/*                                             if FIRST = .TRUE.; */
/*             LZWORK >= MAX( 1, N+M-1, 3*P ), if FIRST = .FALSE. . */
/*             The second term is not needed if M = 0. */
/*             For optimum performance LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the reduction algorithm of [1]. */

/*     REFERENCES */

/*     [1] P. Misra, P. Van Dooren and A. Varga. */
/*         Computation of structural invariants of generalized */
/*         state-space systems. */
/*         Automatica, 30, pp. 1921-1936, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( (P+N)*(M+N)*N )  floating point operations. */

/*     FURTHER COMMENTS */

/*     The number of infinite zeros is computed as */

/*                   DINFZ */
/*        NINFZ =     Sum  (INFZ(i)*i) . */
/*                    i=1 */
/*     Note that each infinite zero of multiplicity k corresponds to */
/*     an infinite eigenvalue of multiplicity k+1. */
/*     The multiplicities of the infinite eigenvalues can be determined */
/*     from PR, DINFZ and INFZ(i), i = 1, ..., DINFZ, as follows: */

/*                     DINFZ */
/*     - there are PR - Sum (INFZ(i)) simple infinite eigenvalues; */
/*                      i=1 */

/*     - there are INFZ(i) infinite eigenvalues with multiplicity i+1, */
/*       for i = 1, ..., DINFZ. */

/*     The left Kronecker indices are: */

/*     [ 0  0 ...  0  | 1  1  ...  1 |  .... | NKRONL  ...  NKRONL ] */
/*     |<- KRONL(1) ->|<- KRONL(2) ->|       |<-  KRONL(NKRONL)  ->| */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     May 1999. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, unitary transformation, structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    abcd_dim1 = *ldabcd;
    abcd_offset = 1 + abcd_dim1;
    abcd -= abcd_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    --infz;
    --kronl;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    lquery = *lzwork == -1;
    *info = 0;
    pn = *p + *n;
    mn = *m + *n;
    mpm = min(*p,*m);
    if (*n < 0) {
	*info = -2;
    } else if (*m < 0 || ! (*first) && *m > *p) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*svlmax < 0.) {
	*info = -5;
    } else if (*ldabcd < max(1,pn)) {
	*info = -7;
    } else if (*lde < max(1,*n)) {
	*info = -9;
    } else if (*tol > 1.) {
	*info = -17;
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = *p * 3;
	wrkopt = max(i__1,i__2);
	if (*p > 0) {
	    if (*m > 0) {
/* Computing MAX */
		i__1 = wrkopt, i__2 = mn - 1;
		wrkopt = max(i__1,i__2);
		if (*first) {
/* Computing MAX */
/* Computing MAX */
		    i__3 = *m * 3 - 1;
		    i__1 = wrkopt, i__2 = mpm + max(i__3,*n);
		    wrkopt = max(i__1,i__2);
		    if (lquery) {
/* Computing MIN */
			i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", p, n,
				 &mpm, &c_n1, (ftnlen)6, (ftnlen)2);
			nb = min(i__1,i__2);
/* Computing MAX */
			i__1 = wrkopt, i__2 = mpm + max(1,*n) * nb;
			wrkopt = max(i__1,i__2);
		    }
		}
	    }
	}
	if (*lzwork < wrkopt && ! lquery) {
	    *info = -21;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("AG8BYZ", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
	return 0;
    }

/*     Initialize output variables. */

    *pr = *p;
    *nr = *n;
    *dinfz = 0;
    *ninfz = 0;
    *nkronl = 0;

/*     Quick return if possible. */

    if (*p == 0) {
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }
    if (*n == 0 && *m == 0) {
	*pr = 0;
	*nkronl = 1;
	kronl[1] = *p;
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }

    rcond = *tol;
    if (rcond <= 0.) {

/*        Use the default tolerance in rank determination. */

	rcond = (doublereal) (pn * mn) * dlamch_("EPSILON", (ftnlen)7);
    }

/*     The D matrix is (RO+SIGMA)-by-M, where RO = P - SIGMA and */
/*     SIGMA = 0 for FIRST = .TRUE. and SIGMA = M for FIRST = .FALSE.. */
/*     The leading (RO+SIGMA)-by-SIGMA submatrix of D has full column */
/*     rank, with the trailing SIGMA-by-SIGMA submatrix upper triangular. */

    if (*first) {
	sigma = 0;
    } else {
	sigma = *m;
    }
    ro = *p - sigma;
    mp1 = *m + 1;
    mui = 0;
    dum[0].r = 0., dum[0].i = 0.;

    itau = 1;
    jwork1 = itau + mpm;
    ismin = 1;
    ismax = ismin + *p;
    jwork2 = ismax + *p;
    nblcks = 0;
    wrkopt = 1;

L10:
    if (*pr == 0) {
	goto L90;
    }

/*     (NR+1,ICOL+1) points to the current position of matrix D. */

    ro1 = ro;
    mnr = *m + *nr;
    if (*m > 0) {

/*        Compress rows of D; first exploit the trapezoidal shape of the */
/*        (RO+SIGMA)-by-SIGMA matrix in the first SIGMA columns of D; */
/*        compress the first SIGMA columns without column pivoting: */

/*              ( x x x x x )       ( x x x x x ) */
/*              ( x x x x x )       ( 0 x x x x ) */
/*              ( x x x x x )  - >  ( 0 0 x x x ) */
/*              ( 0 x x x x )       ( 0 0 0 x x ) */
/*              ( 0 0 x x x )       ( 0 0 0 x x ) */

/*        where SIGMA = 3 and RO = 2. */
/*        Complex workspace: need maximum M+N-1. */

	irow = *nr;
	i__1 = sigma;
	for (icol = 1; icol <= i__1; ++icol) {
	    ++irow;
	    i__2 = ro + 1;
	    zlarfg_(&i__2, &abcd[irow + icol * abcd_dim1], &abcd[irow + 1 + 
		    icol * abcd_dim1], &c__1, &tc);
	    i__2 = ro + 1;
	    i__3 = mnr - icol;
	    d_cnjg(&z__1, &tc);
	    zlatzm_("L", &i__2, &i__3, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1, &z__1, &abcd[irow + (icol + 1) * abcd_dim1], &abcd[
		    irow + 1 + (icol + 1) * abcd_dim1], ldabcd, &zwork[1], (
		    ftnlen)1);
	    i__2 = *pr - icol;
	    zcopy_(&i__2, dum, &c__0, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1);
/* L20: */
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = mn - 1;
	wrkopt = max(i__1,i__2);

	if (*first) {

/*           Continue with Householder with column pivoting. */

/*              ( x x x x x )        ( x x x x x ) */
/*              ( 0 x x x x )        ( 0 x x x x ) */
/*              ( 0 0 x x x )  - >   ( 0 0 x x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 0 0 ) */

/*           Real workspace:    need maximum 2*M; */
/*           Complex workspace: need maximum min(P,M)+3*M-1; */
/*           Integer workspace: need maximum M. */

/* Computing MIN */
	    i__1 = *nr + sigma + 1;
	    irow = min(i__1,pn);
/* Computing MIN */
	    i__1 = sigma + 1;
	    icol = min(i__1,*m);
	    i__1 = *m - sigma;
	    mb3oyz_(&ro1, &i__1, &abcd[irow + icol * abcd_dim1], ldabcd, &
		    rcond, svlmax, &rank, sval, &iwork[1], &zwork[itau], &
		    dwork[1], &zwork[jwork1], info);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = jwork1 + *m * 3 - 2;
	    wrkopt = max(i__1,i__2);

/*           Apply the column permutations to B and part of D. */

	    i__1 = *nr + sigma;
	    i__2 = *m - sigma;
	    zlapmt_(&c_true, &i__1, &i__2, &abcd[icol * abcd_dim1 + 1], 
		    ldabcd, &iwork[1]);

	    if (rank > 0) {

/*              Apply the Householder transformations to the submatrix C. */
/*              Complex workspace: need   maximum min(P,M) + N; */
/*                                 prefer maximum min(P,M) + N*NB. */

		i__1 = *lzwork - jwork1 + 1;
		zunmqr_("Left", "ConjTranspose", &ro1, nr, &rank, &abcd[irow 
			+ icol * abcd_dim1], ldabcd, &zwork[itau], &abcd[irow 
			+ mp1 * abcd_dim1], ldabcd, &zwork[jwork1], &i__1, 
			info, (ftnlen)4, (ftnlen)13);
/* Computing MAX */
		i__3 = jwork1;
		i__1 = wrkopt, i__2 = jwork1 + (integer) zwork[i__3].r - 1;
		wrkopt = max(i__1,i__2);
		i__1 = ro1 - 1;
/* Computing MIN */
		i__3 = ro1 - 1;
		i__2 = min(i__3,rank);
/* Computing MIN */
		i__4 = irow + 1;
		zlaset_("Lower", &i__1, &i__2, &c_b2, &c_b2, &abcd[min(i__4,
			pn) + icol * abcd_dim1], ldabcd, (ftnlen)5);
		ro1 -= rank;
	    }
	}

/*        Terminate if Dr has maximal row rank. */

	if (ro1 == 0) {
	    goto L90;
	}

    }

/*     Update SIGMA. */

    sigma = *pr - ro1;

    ++nblcks;
    taui = ro1;

/*     Compress the columns of current C to separate a TAUI-by-MUI */
/*     full column rank block. */

    if (*nr == 0) {

/*        Finish for zero state dimension. */

	*pr = sigma;
	rank = 0;
    } else {

/*        Perform RQ-decomposition with row pivoting on the current C */
/*        while keeping E upper triangular. */
/*        The current C is the TAUI-by-NR matrix delimited by rows */
/*        IRC+1 to IRC+TAUI and columns M+1 to M+NR of ABCD. */
/*        The rank of current C is computed in MUI. */
/*        Real workspace:    need maximum 2*P; */
/*        Complex workspace: need maximum 3*P. */

	irc = *nr + sigma;
	n1 = *nr;
	if (taui > 1) {

/*           Compute norms. */

	    i__1 = taui;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dwork[i__] = dznrm2_(nr, &abcd[irc + i__ + mp1 * abcd_dim1], 
			ldabcd);
		dwork[*p + i__] = dwork[i__];
/* L30: */
	    }
	}

	rank = 0;
	mntau = min(taui,*nr);

/*        ICOL and IROW will point to the current pivot position in C. */

	ilast = *nr + *pr;
	jlast = *m + *nr;
	irow = ilast;
	icol = jlast;
	i__ = taui;
L40:
	if (rank < mntau) {
	    mn1 = *m + n1;

/*           Pivot if necessary. */

	    if (i__ != 1) {
		j = idamax_(&i__, &dwork[1], &c__1);
		if (j != i__) {
		    dwork[j] = dwork[i__];
		    dwork[*p + j] = dwork[*p + i__];
		    zswap_(&n1, &abcd[irow + mp1 * abcd_dim1], ldabcd, &abcd[
			    irc + j + mp1 * abcd_dim1], ldabcd);
		}
	    }

/*           Zero elements left to ABCD(IROW,ICOL). */

	    i__1 = n1 - 1;
	    for (k = 1; k <= i__1; ++k) {
		j = *m + k;

/*              Rotate columns J, J+1 to zero ABCD(IROW,J). */

		i__2 = irow + (j + 1) * abcd_dim1;
		tc.r = abcd[i__2].r, tc.i = abcd[i__2].i;
		zlartg_(&tc, &abcd[irow + j * abcd_dim1], &c__, &s, &abcd[
			irow + (j + 1) * abcd_dim1]);
		i__2 = irow + j * abcd_dim1;
		abcd[i__2].r = 0., abcd[i__2].i = 0.;
		i__2 = irow - 1;
		zrot_(&i__2, &abcd[(j + 1) * abcd_dim1 + 1], &c__1, &abcd[j * 
			abcd_dim1 + 1], &c__1, &c__, &s);
		i__2 = k + 1;
		zrot_(&i__2, &e[(k + 1) * e_dim1 + 1], &c__1, &e[k * e_dim1 + 
			1], &c__1, &c__, &s);

/*              Rotate rows K, K+1 to zero E(K+1,K). */

		i__2 = k + k * e_dim1;
		tc.r = e[i__2].r, tc.i = e[i__2].i;
		zlartg_(&tc, &e[k + 1 + k * e_dim1], &c__, &s, &e[k + k * 
			e_dim1]);
		i__2 = k + 1 + k * e_dim1;
		e[i__2].r = 0., e[i__2].i = 0.;
		i__2 = n1 - k;
		zrot_(&i__2, &e[k + (k + 1) * e_dim1], lde, &e[k + 1 + (k + 1)
			 * e_dim1], lde, &c__, &s);
		zrot_(&mn1, &abcd[k + abcd_dim1], ldabcd, &abcd[k + 1 + 
			abcd_dim1], ldabcd, &c__, &s);
/* L50: */
	    }

	    if (rank == 0) {

/*              Initialize; exit if matrix is zero (RANK = 0). */

		smax = z_abs(&abcd[ilast + jlast * abcd_dim1]);
		if (smax == 0.) {
		    goto L80;
		}
		smin = smax;
		smaxpr = smax;
		sminpr = smin;
		c1.r = 1., c1.i = 0.;
		c2.r = 1., c2.i = 0.;
	    } else {

/*              One step of incremental condition estimation. */
/*              Complex workspace: need maximum 3*P. */

		zcopy_(&rank, &abcd[irow + (icol + 1) * abcd_dim1], ldabcd, &
			zwork[jwork2], &c__1);
		zlaic1_(&c__2, &rank, &zwork[ismin], &smin, &zwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &sminpr, &s1, &c1);
		zlaic1_(&c__1, &rank, &zwork[ismax], &smax, &zwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &smaxpr, &s2, &c2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = *p * 3;
		wrkopt = max(i__1,i__2);
	    }

/*           Check the rank; finish the loop if rank loss occurs. */

	    if (*svlmax * rcond <= smaxpr) {
		if (*svlmax * rcond <= sminpr) {
		    if (smaxpr * rcond <= sminpr) {

/*                    Finish the loop if last row. */

			if (n1 == 0) {
			    ++rank;
			    goto L80;
			}

			if (n1 > 1) {

/*                       Update norms. */

			    if (i__ - 1 > 1) {
				i__1 = i__ - 1;
				for (j = 1; j <= i__1; ++j) {
				    if (dwork[j] != 0.) {
/* Computing 2nd power */
					d__1 = z_abs(&abcd[irc + j + icol * 
						abcd_dim1]) / dwork[j];
					t = 1. - d__1 * d__1;
					t = max(t,0.);
/* Computing 2nd power */
					d__1 = dwork[j] / dwork[*p + j];
					tt = t * .05 * (d__1 * d__1) + 1.;
					if (tt != 1.) {
					    dwork[j] *= sqrt(t);
					} else {
					    i__2 = n1 - 1;
					    dwork[j] = dznrm2_(&i__2, &abcd[
						    irc + j + mp1 * abcd_dim1]
						    , ldabcd);
					    dwork[*p + j] = dwork[j];
					}
				    }
/* L60: */
				}
			    }
			}

			i__1 = rank;
			for (j = 1; j <= i__1; ++j) {
			    i__2 = ismin + j - 1;
			    i__3 = ismin + j - 1;
			    z__1.r = s1.r * zwork[i__3].r - s1.i * zwork[i__3]
				    .i, z__1.i = s1.r * zwork[i__3].i + s1.i *
				     zwork[i__3].r;
			    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
			    i__2 = ismax + j - 1;
			    i__3 = ismax + j - 1;
			    z__1.r = s2.r * zwork[i__3].r - s2.i * zwork[i__3]
				    .i, z__1.i = s2.r * zwork[i__3].i + s2.i *
				     zwork[i__3].r;
			    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L70: */
			}

			i__1 = ismin + rank;
			zwork[i__1].r = c1.r, zwork[i__1].i = c1.i;
			i__1 = ismax + rank;
			zwork[i__1].r = c2.r, zwork[i__1].i = c2.i;
			smin = sminpr;
			smax = smaxpr;
			++rank;
			--icol;
			--irow;
			--n1;
			--i__;
			goto L40;
		    }
		}
	    }
	}
    }

L80:
    mui = rank;
    *nr -= mui;
    *pr = sigma + mui;

/*     Set number of left Kronecker blocks of order (i-1)-by-i. */

    kronl[nblcks] = taui - mui;

/*     Set number of infinite divisors of order i-1. */

    if (*first && nblcks > 1) {
	infz[nblcks - 1] = muim1 - taui;
    }
    muim1 = mui;
    ro = mui;

/*     Continue reduction if rank of current C is positive. */

    if (mui > 0) {
	goto L10;
    }

/*     Determine the maximal degree of infinite zeros and */
/*     the number of infinite zeros. */

L90:
    if (*first) {
	if (mui == 0) {
/* Computing MAX */
	    i__1 = 0, i__2 = nblcks - 1;
	    *dinfz = max(i__1,i__2);
	} else {
	    *dinfz = nblcks;
	    infz[nblcks] = mui;
	}
	k = *dinfz;
	for (i__ = k; i__ >= 1; --i__) {
	    if (infz[i__] != 0) {
		goto L110;
	    }
	    --(*dinfz);
/* L100: */
	}
L110:
	i__1 = *dinfz;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    *ninfz += infz[i__] * i__;
/* L120: */
	}
    }

/*     Determine the maximal order of left elementary Kronecker blocks. */

    *nkronl = nblcks;
    for (i__ = nblcks; i__ >= 1; --i__) {
	if (kronl[i__] != 0) {
	    goto L140;
	}
	--(*nkronl);
/* L130: */
    }
L140:

    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
    return 0;
/* *** Last line of AG8BYZ *** */
} /* ag8byz_ */

