/* MB03XP.f -- translated by f2c (version 20100827).
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

static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__8 = 8;
static integer c__3 = 3;
static integer c__15 = 15;
static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b89 = -1.;

/* Subroutine */ int mb03xp_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *q, integer *ldq, doublereal *z__, 
	integer *ldz, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	job_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3[2], i__4, i__5;
    doublereal d__1, d__2;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal v[51];
    static integer i1, i2;
    static doublereal as[225]	/* was [15][15] */, bs[225]	/* was [15][
	    15] */;
    static integer kk, nh, nr, ns, nv, pv2, pv3, dum, itn, its;
    static doublereal ulp, tst;
    static integer maxb, ierr;
    static doublereal unfl, temp, ovfl, tauv, tauw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iseed[4];
    extern integer ue01md_(integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03ya_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), mb03yd_(logical *,
	     logical *, logical *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical initq;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical wantq;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical initz, wantt, wantz;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dlarfx_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen), 
	    dlarnv_(integer *, integer *, integer *, doublereal *);
    static doublereal smlnum;


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

/*     To compute the periodic Schur decomposition and the eigenvalues of */
/*     a product of matrices, H = A*B, with A upper Hessenberg and B */
/*     upper triangular without evaluating any part of the product. */
/*     Specifically, the matrices Q and Z are computed, so that */

/*          Q' * A * Z = S,    Z' * B * Q = T */

/*     where S is in real Schur form, and T is upper triangular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = 'E':  Compute the eigenvalues only; */
/*             = 'S':  compute the factors S and T of the full */
/*                     Schur form. */

/*     COMPQ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = 'N':  The matrix Q is not required; */
/*             = 'I':  Q is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Q is returned; */
/*             = 'V':  Q must contain an orthogonal matrix U on entry, */
/*                     and the product U*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = 'N':  The matrix Z is not required; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned; */
/*             = 'V':  Z must contain an orthogonal matrix U on entry, */
/*                     and the product U*Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B. N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that the matrices A and B are already upper */
/*             triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/*             The routine works primarily with the submatrices in rows */
/*             and columns ILO to IHI, but applies the transformations to */
/*             all the rows and columns of the matrices A and B, if */
/*             JOB = 'S'. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array A must */
/*             contain the upper Hessenberg matrix A. */
/*             On exit, if JOB = 'S', the leading N-by-N part of this */
/*             array is upper quasi-triangular with any 2-by-2 diagonal */
/*             blocks corresponding to a pair of complex conjugated */
/*             eigenvalues. */
/*             If JOB = 'E', the diagonal elements and 2-by-2 diagonal */
/*             blocks of A will be correct, but the remaining parts of A */
/*             are unspecified on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array B must */
/*             contain the upper triangular matrix B. */
/*             On exit, if JOB = 'S', the leading N-by-N part of this */
/*             array contains the transformed upper triangular matrix. */
/*             2-by-2 blocks in B corresponding to 2-by-2 blocks in A */
/*             will be reduced to positive diagonal form. (I.e., if */
/*             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) */
/*             and B(j+1,j+1) will be positive.) */
/*             If JOB = 'E', the elements corresponding to diagonal */
/*             elements and 2-by-2 diagonal blocks in A will be correct, */
/*             but the remaining parts of B are unspecified on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if COMPQ = 'V', then the leading N-by-N part of */
/*             this array must contain a matrix Q which is assumed to be */
/*             equal to the unit matrix except for the submatrix */
/*             Q(ILO:IHI,ILO:IHI). */
/*             If COMPQ = 'I', Q need not be set on entry. */
/*             On exit, if COMPQ = 'V' or COMPQ = 'I' the leading N-by-N */
/*             part of this array contains the transformation matrix */
/*             which produced the Schur form. */
/*             If COMPQ = 'N', Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= 1. */
/*             If COMPQ <> 'N', LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if COMPZ = 'V', then the leading N-by-N part of */
/*             this array must contain a matrix Z which is assumed to be */
/*             equal to the unit matrix except for the submatrix */
/*             Z(ILO:IHI,ILO:IHI). */
/*             If COMPZ = 'I', Z need not be set on entry. */
/*             On exit, if COMPZ = 'V' or COMPZ = 'I' the leading N-by-N */
/*             part of this array contains the transformation matrix */
/*             which produced the Schur form. */
/*             If COMPZ = 'N', Z is not referenced. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= 1. */
/*             If COMPZ <> 'N', LDZ >= MAX(1,N). */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             The i-th (1 <= i <= N) computed eigenvalue is given by */
/*             BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two */
/*             eigenvalues are computed as a complex conjugate pair, */
/*             they are stored in consecutive elements of ALPHAR, ALPHAI */
/*             and BETA. If JOB = 'S', the eigenvalues are stored in the */
/*             same order as on the diagonales of the Schur forms of A */
/*             and B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then MB03XP failed to compute the Schur */
/*                   form in a total of 30*(IHI-ILO+1) iterations; */
/*                   elements 1:ilo-1 and i+1:n of ALPHAR, ALPHAI and */
/*                   BETA contain successfully computed eigenvalues. */

/*     METHOD */

/*     The implemented algorithm is a multi-shift version of the periodic */
/*     QR algorithm described in [1,3] with some minor modifications */
/*     proposed in [2]. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P. */
/*         The periodic Schur decomposition: Algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Kressner, D. */
/*         An efficient and reliable implementation of the periodic QZ */
/*         algorithm. Proc. of the IFAC Workshop on Periodic Control */
/*         Systems, pp. 187-192, 2001. */

/*     [3] Van Loan, C. */
/*         Generalized Singular Values with Algorithms and Applications. */
/*         Ph. D. Thesis, University of Michigan, 1973. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**3) floating point operations and is */
/*     backward stable. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHGPQR). */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, Hessenberg form, orthogonal */
/*     transformation, (periodic) Schur form */

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

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --alphar;
    --alphai;
    --beta;
    --dwork;

    /* Function Body */
    wantt = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    *info = 0;
    if (! lsame_(job, "E", (ftnlen)1, (ftnlen)1) && ! wantt) {
	*info = -1;
    } else if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
	*info = -2;
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -5;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldq < 1 || wantq && *ldq < *n) {
	*info = -12;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
	*info = -14;
    } else if (*ldwork < max(1,*n)) {
	dwork[1] = (doublereal) max(1,*n);
	*info = -19;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03XP", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize Q and Z, if necessary. */

    if (initq) {
	dlaset_("All", n, n, &c_b12, &c_b13, &q[q_offset], ldq, (ftnlen)3);
    }
    if (initz) {
	dlaset_("All", n, n, &c_b12, &c_b13, &z__[z_offset], ldz, (ftnlen)3);
    }

/*     Store isolated eigenvalues and standardize B. */

/*     FOR I = [1:ILO-1, IHI+1:N] */
    i__ = 1;
L10:
    if (i__ == *ilo) {
	i__ = *ihi + 1;
    }
    if (i__ <= *n) {
	if (b[i__ + i__ * b_dim1] < 0.) {
	    if (wantt) {
		i__1 = i__;
		for (k = *ilo; k <= i__1; ++k) {
		    b[k + i__ * b_dim1] = -b[k + i__ * b_dim1];
/* L20: */
		}
		i__1 = *ihi;
		for (k = i__; k <= i__1; ++k) {
		    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1];
/* L30: */
		}
	    } else {
		b[i__ + i__ * b_dim1] = -b[i__ + i__ * b_dim1];
		a[i__ + i__ * a_dim1] = -a[i__ + i__ * a_dim1];
	    }
	    if (wantq) {
		i__1 = *ihi;
		for (k = *ilo; k <= i__1; ++k) {
		    q[k + i__ * q_dim1] = -q[k + i__ * q_dim1];
/* L40: */
		}
	    }
	}
	alphar[i__] = a[i__ + i__ * a_dim1];
	alphai[i__] = 0.;
	beta[i__] = b[i__ + i__ * b_dim1];
	++i__;
/*        END FOR */
	goto L10;
    }

/*     Quick return if possible. */

    if (*n == 0 || *ilo == *ihi + 1) {
	dwork[1] = 1.;
	return 0;
    }

/*     Set rows and coloms ILO to IHI of B (A) to zero below the first */
/*     (sub)diagonal. */

    i__1 = *ihi - 2;
    for (j = *ilo; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j + 2; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] = 0.;
/* L50: */
	}
/* L60: */
    }
    i__1 = *ihi - 1;
    for (j = *ilo; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    b[i__ + j * b_dim1] = 0.;
/* L70: */
	}
/* L80: */
    }
    nh = *ihi - *ilo + 1;

/*     Suboptimal choice of the number of shifts. */

    if (wantq) {
/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compq;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	ns = ue01md_(&c__4, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (ftnlen)
		2);
/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compq;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	maxb = ue01md_(&c__8, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (
		ftnlen)2);
    } else {
/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compz;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	ns = ue01md_(&c__4, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (ftnlen)
		2);
/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compz;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	maxb = ue01md_(&c__8, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (
		ftnlen)2);
    }

    if (ns <= 2 || ns > nh || maxb >= nh) {

/*        Standard double-shift product QR. */

	mb03yd_(&wantt, &wantq, &wantz, n, ilo, ihi, ilo, ihi, &a[a_offset], 
		lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], 
		ldz, &alphar[1], &alphai[1], &beta[1], &dwork[1], ldwork, 
		info);
	return 0;
    }
    maxb = max(3,maxb);
/* Computing MIN */
    i__1 = min(ns,maxb);
    ns = min(i__1,15);

/*     Set machine-dependent constants for the stopping criterion. */
/*     If max(norm(A),norm(B)) <= sqrt(OVFL), then overflow should not */
/*     occur. */

    unfl = dlamch_("Safe minimum", (ftnlen)12);
    ovfl = 1. / unfl;
    dlabad_(&unfl, &ovfl);
    ulp = dlamch_("Precision", (ftnlen)9);
    smlnum = unfl * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first rows and last columns of */
/*     A and B to which transformations must be applied. */

    if (wantt) {
	i1 = 1;
	i2 = *n;
    }
    iseed[0] = 1;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 1;

/*     ITN is the maximal number of QR iterations. */

    itn = nh * 30;
    dum = 0;

/*     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO */
/*     or A(L,L-1) is negligible. */

    i__ = *ihi;
L90:
    l = *ilo;
    if (i__ < *ilo) {
	goto L210;
    }

    i__1 = itn;
    for (its = 0; its <= i__1; ++its) {
	dum += (*ihi - *ilo) * (*ihi - *ilo);

/*        Look for deflations in A. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    tst = (d__1 = a[k - 1 + (k - 1) * a_dim1], abs(d__1)) + (d__2 = a[
		    k + k * a_dim1], abs(d__2));
	    if (tst == 0.) {
		i__4 = i__ - l + 1;
		tst = dlanhs_("1", &i__4, &a[l + l * a_dim1], lda, &dwork[1], 
			(ftnlen)1);
	    }
/* Computing MAX */
	    d__2 = ulp * tst;
	    if ((d__1 = a[k + (k - 1) * a_dim1], abs(d__1)) <= max(d__2,
		    smlnum)) {
		goto L110;
	    }
/* L100: */
	}
L110:

/*        Look for deflation in B if problem size is greater than 1. */

	if (i__ - k >= 1) {
	    i__2 = k;
	    for (kk = i__; kk >= i__2; --kk) {
		if (kk == i__) {
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1));
		} else if (kk == k) {
		    tst = (d__1 = b[kk + (kk + 1) * b_dim1], abs(d__1));
		} else {
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1)) + (d__2 
			    = b[kk + (kk + 1) * b_dim1], abs(d__2));
		}
		if (tst == 0.) {
		    i__4 = i__ - k + 1;
		    tst = dlanhs_("1", &i__4, &b[k + k * b_dim1], ldb, &dwork[
			    1], (ftnlen)1);
		}
/* Computing MAX */
		d__2 = ulp * tst;
		if ((d__1 = b[kk + kk * b_dim1], abs(d__1)) <= max(d__2,
			smlnum)) {
		    goto L130;
		}
/* L120: */
	    }
	} else {
	    kk = k - 1;
	}
L130:
	if (kk >= k) {

/*           B has an element close to zero at position (KK,KK). */

	    b[kk + kk * b_dim1] = 0.;
	    mb03ya_(&wantt, &wantq, &wantz, n, &k, &i__, ilo, ihi, &kk, &a[
		    a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
		    z__[z_offset], ldz, info);
	    k = kk + 1;
	}
	l = k;
	if (l > *ilo) {

/*           A(L,L-1) is negligible. */

	    a[l + (l - 1) * a_dim1] = 0.;
	}

/*        Exit from loop if a submatrix of order <= MAXB has split off. */

	if (l >= i__ - maxb + 1) {
	    goto L200;
	}

/*        The active submatrices are now in rows and columns L:I. */

	if (! wantt) {
	    i1 = l;
	    i2 = i__;
	}
	if (its == 10 || its == 20) {

/*           Exceptional shift. The first column of the shift polynomial */
/*           is a pseudo-random vector. */

	    i__2 = ns + 1;
	    dlarnv_(&c__3, iseed, &i__2, v);
	} else {

/*           Use eigenvalues of trailing submatrix as shifts. */

	    dlacpy_("Full", &ns, &ns, &a[i__ - ns + 1 + (i__ - ns + 1) * 
		    a_dim1], lda, as, &c__15, (ftnlen)4);
	    dlacpy_("Full", &ns, &ns, &b[i__ - ns + 1 + (i__ - ns + 1) * 
		    b_dim1], ldb, bs, &c__15, (ftnlen)4);
	    mb03yd_(&c_false, &c_false, &c_false, &ns, &c__1, &ns, &c__1, &ns,
		     as, &c__15, bs, &c__15, &q[q_offset], ldq, &z__[z_offset]
		    , ldz, &alphar[i__ - ns + 1], &alphai[i__ - ns + 1], &
		    beta[i__ - ns + 1], &dwork[1], ldwork, &ierr);
	}

/*        Compute the nonzero elements of the first column of */
/*        (A*B-w(1)) (A*B-w(2)) .. (A*B-w(ns)). */

	v[0] = 1.;
	nv = 1;
/*        WHILE NV <= NS */
L140:
	if (nv <= ns) {
	    if (nv == ns || as[nv + 1 + nv * 15 - 16] == 0.) {

/*              Real shift. */

		v[nv] = 0.;
		pv2 = nv + 2;
		dcopy_(&nv, v, &c__1, &v[pv2 - 1], &c__1);
		dtrmv_("Upper", "No transpose", "No unit diagonal", &nv, &b[l 
			+ l * b_dim1], ldb, &v[pv2 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
		dscal_(&nv, &bs[nv + nv * 15 - 16], v, &c__1);
		i__2 = (nv << 1) + 1;
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
		temp = 1. / max(d__2,smlnum);
		i__2 = (nv << 1) + 1;
		dscal_(&i__2, &temp, v, &c__1);
		i__2 = nv + 1;
		d__1 = -as[nv + nv * 15 - 16];
		dgemv_("No transpose", &i__2, &nv, &c_b13, &a[l + l * a_dim1],
			 lda, &v[pv2 - 1], &c__1, &d__1, v, &c__1, (ftnlen)12)
			;
		++nv;
	    } else {

/*              Double shift using a product formulation of the shift */
/*              polynomial [2]. */

		v[nv] = 0.;
		v[nv + 1] = 0.;
		pv2 = nv + 3;
		pv3 = (nv << 1) + 5;
		i__2 = nv + 2;
		dcopy_(&i__2, v, &c__1, &v[pv2 - 1], &c__1);
		i__2 = nv + 1;
		dcopy_(&i__2, v, &c__1, &v[pv3 - 1], &c__1);
		dscal_(&nv, &bs[nv + 1 + (nv + 1) * 15 - 16], &v[pv2 - 1], &
			c__1);
		dtrmv_("Upper", "No transpose", "No unit diagonal", &nv, &b[l 
			+ l * b_dim1], ldb, &v[pv3 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
		i__2 = (nv << 1) + 3;
		itemp = idamax_(&i__2, &v[pv2 - 1], &c__1);
/* Computing MAX */
		d__2 = (d__1 = v[pv2 + itemp - 2], abs(d__1));
		temp = 1. / max(d__2,smlnum);
		i__2 = (nv << 1) + 3;
		dscal_(&i__2, &temp, &v[pv2 - 1], &c__1);

		dcopy_(&nv, &v[pv2 - 1], &c__1, v, &c__1);
		i__2 = nv + 1;
		dgemv_("No transpose", &i__2, &nv, &c_b89, &a[l + l * a_dim1],
			 lda, &v[pv3 - 1], &c__1, &as[nv + 1 + (nv + 1) * 15 
			- 16], &v[pv2 - 1], &c__1, (ftnlen)12);
		dscal_(&nv, &as[nv + (nv + 1) * 15 - 16], v, &c__1);
		i__2 = (nv << 1) + 3;
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
		temp = 1. / max(d__2,smlnum);
		i__2 = (nv << 1) + 3;
		dscal_(&i__2, &temp, v, &c__1);

		d__1 = -as[nv + 1 + nv * 15 - 16];
		dscal_(&nv, &d__1, v, &c__1);
		i__2 = nv + 1;
		daxpy_(&i__2, &as[nv + nv * 15 - 16], &v[pv2 - 1], &c__1, v, &
			c__1);
		i__2 = (nv << 1) + 3;
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
		temp = 1. / max(d__2,smlnum);
		i__2 = (nv << 1) + 3;
		dscal_(&i__2, &temp, v, &c__1);

		i__2 = nv + 1;
		dscal_(&i__2, &bs[nv + nv * 15 - 16], v, &c__1);
		i__2 = nv + 1;
		dtrmv_("Upper", "No transpose", "No unit diagonal", &i__2, &b[
			l + l * b_dim1], ldb, &v[pv2 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
		i__2 = (nv << 1) + 3;
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
		temp = 1. / max(d__2,smlnum);
		i__2 = (nv << 1) + 3;
		dscal_(&i__2, &temp, v, &c__1);

		i__2 = nv + 2;
		i__4 = nv + 1;
		dgemv_("No transpose", &i__2, &i__4, &c_b89, &a[l + l * 
			a_dim1], lda, &v[pv2 - 1], &c__1, &c_b13, v, &c__1, (
			ftnlen)12);
		nv += 2;
	    }
	    itemp = idamax_(&nv, v, &c__1);
	    temp = (d__1 = v[itemp - 1], abs(d__1));
	    if (temp == 0.) {
		v[0] = 1.;
		i__2 = nv;
		for (k = 2; k <= i__2; ++k) {
		    v[k - 1] = 0.;
/* L150: */
		}
	    } else {
		temp = max(temp,smlnum);
		d__1 = 1. / temp;
		dscal_(&nv, &d__1, v, &c__1);
	    }
	    goto L140;
/*        END WHILE */
	}

/*        Multi-shift product QR step. */

	pv2 = ns + 2;
	i__2 = i__ - 1;
	for (k = l; k <= i__2; ++k) {
/* Computing MIN */
	    i__4 = ns + 1, i__5 = i__ - k + 1;
	    nr = min(i__4,i__5);
	    if (k > l) {
		dcopy_(&nr, &a[k + (k - 1) * a_dim1], &c__1, v, &c__1);
	    }
	    dlarfg_(&nr, v, &v[1], &c__1, &tauv);
	    if (k > l) {
		a[k + (k - 1) * a_dim1] = v[0];
		i__4 = i__;
		for (kk = k + 1; kk <= i__4; ++kk) {
		    a[kk + (k - 1) * a_dim1] = 0.;
/* L160: */
		}
	    }

/*           Apply reflector V from the right to B in rows */
/*           I1:min(K+NS,I). */

	    v[0] = 1.;
/* Computing MIN */
	    i__5 = k + ns;
	    i__4 = min(i__5,i__) - i1 + 1;
	    dlarfx_("Right", &i__4, &nr, v, &tauv, &b[i1 + k * b_dim1], ldb, &
		    dwork[1], (ftnlen)5);

/*           Annihilate the introduced nonzeros in the K-th column. */

	    dcopy_(&nr, &b[k + k * b_dim1], &c__1, &v[pv2 - 1], &c__1);
	    dlarfg_(&nr, &v[pv2 - 1], &v[pv2], &c__1, &tauw);
	    b[k + k * b_dim1] = v[pv2 - 1];
	    i__4 = i__;
	    for (kk = k + 1; kk <= i__4; ++kk) {
		b[kk + k * b_dim1] = 0.;
/* L170: */
	    }
	    v[pv2 - 1] = 1.;

/*           Apply reflector W from the left to transform the rows of the */
/*           matrix B in columns K+1:I2. */

	    i__4 = i2 - k;
	    dlarfx_("Left", &nr, &i__4, &v[pv2 - 1], &tauw, &b[k + (k + 1) * 
		    b_dim1], ldb, &dwork[1], (ftnlen)4);

/*           Apply reflector V from the left to transform the rows of the */
/*           matrix A in columns K:I2. */

	    i__4 = i2 - k + 1;
	    dlarfx_("Left", &nr, &i__4, v, &tauv, &a[k + k * a_dim1], lda, &
		    dwork[1], (ftnlen)4);

/*           Apply reflector W from the right to transform the columns of */
/*           the matrix A in rows I1:min(K+NS,I). */

/* Computing MIN */
	    i__5 = k + ns + 1;
	    i__4 = min(i__5,i__) - i1 + 1;
	    dlarfx_("Right", &i__4, &nr, &v[pv2 - 1], &tauw, &a[i1 + k * 
		    a_dim1], lda, &dwork[1], (ftnlen)5);

/*           Accumulate transformations in the matrices Q and Z. */

	    if (wantq) {
		dlarfx_("Right", &nh, &nr, v, &tauv, &q[*ilo + k * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
	    }
	    if (wantz) {
		dlarfx_("Right", &nh, &nr, &v[pv2 - 1], &tauw, &z__[*ilo + k *
			 z_dim1], ldz, &dwork[1], (ftnlen)5);
	    }
/* L180: */
	}
/* L190: */
    }

/*     Failure to converge. */

    *info = i__;
    return 0;
L200:

/*     Submatrix of order <= MAXB has split off. Use double-shift */
/*     periodic QR algorithm. */

    mb03yd_(&wantt, &wantq, &wantz, n, &l, &i__, ilo, ihi, &a[a_offset], lda, 
	    &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
	    alphar[1], &alphai[1], &beta[1], &dwork[1], ldwork, info);
    if (*info > 0) {
	return 0;
    }
    itn -= its;
    i__ = l - 1;
    goto L90;

L210:
    dwork[1] = (doublereal) max(1,*n);
    return 0;
/* *** Last line of MB03XP *** */
} /* mb03xp_ */

