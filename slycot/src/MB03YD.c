/* MB03YD.f -- translated by f2c (version 20100827).
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

static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int mb03yd_(logical *wantt, logical *wantq, logical *wantz, 
	integer *n, integer *ilo, integer *ihi, integer *iloq, integer *ihiq, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *dwork, integer *
	ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal v[3], w[3];
    static integer i1, i2, kk, nh, nq, nr;
    static doublereal cs1, cs2, cs3, sn1, sn2, sn3;
    static integer itn, its;
    static doublereal ulp, tst, temp, ovfl, unfl;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal tauv, tauw, gamma, alpha, delta;
    static integer iseed[4];
    extern /* Subroutine */ int mb03ya_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal betax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), mb03yt_(doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen), dlarfx_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dlarnv_(integer *, integer *, integer *, doublereal *);
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

/*     To deal with small subtasks of the product eigenvalue problem. */

/*     MB03YD is an auxiliary routine called by SLICOT Library routine */
/*     MB03XP. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WANTT   LOGICAL */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = .TRUE. :  Compute the full Schur form; */
/*             = .FALSE.:  compute the eigenvalues only. */

/*     WANTQ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = .TRUE. :  The matrix Q is updated; */
/*             = .FALSE.:  the matrix Q is not required. */

/*     WANTZ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = .TRUE. :  The matrix Z is updated; */
/*             = .FALSE.:  the matrix Z is not required. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B. N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that the matrices A and B are already */
/*             (quasi) upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N. The routine works primarily with the submatrices */
/*             in rows and columns ILO to IHI, but applies the */
/*             transformations to all the rows and columns of the */
/*             matrices A and B, if WANTT = .TRUE.. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     ILOQ    (input) INTEGER */
/*     IHIQ    (input) INTEGER */
/*             Specify the rows of Q and Z to which transformations */
/*             must be applied if WANTQ = .TRUE. and WANTZ = .TRUE., */
/*             respectively. */
/*             1 <= ILOQ <= ILO; IHI <= IHIQ <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper Hessenberg matrix A. */
/*             On exit, if WANTT = .TRUE., the leading N-by-N part of */
/*             this array is upper quasi-triangular in rows and columns */
/*             ILO:IHI. */
/*             If WANTT = .FALSE., the diagonal elements and 2-by-2 */
/*             diagonal blocks of A will be correct, but the remaining */
/*             parts of A are unspecified on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix B. */
/*             On exit, if WANTT = .TRUE., the leading N-by-N part of */
/*             this array contains the transformed upper triangular */
/*             matrix. 2-by-2 blocks in B corresponding to 2-by-2 blocks */
/*             in A will be reduced to positive diagonal form. (I.e., if */
/*             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) */
/*             and B(j+1,j+1) will be positive.) */
/*             If WANTT = .FALSE., the elements corresponding to diagonal */
/*             elements and 2-by-2 diagonal blocks in A will be correct, */
/*             but the remaining parts of B are unspecified on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Q of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Q updated in the */
/*             submatrix Q(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTQ = .FALSE., Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= 1. */
/*             If WANTQ = .TRUE., LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Z of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Z updated in the */
/*             submatrix Z(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTZ = .FALSE., Z is not referenced. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= 1. */
/*             If WANTZ = .TRUE., LDZ >= MAX(1,N). */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             The i-th (ILO <= i <= IHI) computed eigenvalue is given */
/*             by BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two */
/*             eigenvalues are computed as a complex conjugate pair, */
/*             they are stored in consecutive elements of ALPHAR, ALPHAI */
/*             and BETA. If WANTT = .TRUE., the eigenvalues are stored in */
/*             the same order as on the diagonals of the Schur forms of */
/*             A and B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then MB03YD failed to compute the Schur */
/*                   form in a total of 30*(IHI-ILO+1) iterations; */
/*                   elements i+1:n of ALPHAR, ALPHAI and BETA contain */
/*                   successfully computed eigenvalues. */

/*     METHOD */

/*     The implemented algorithm is a double-shift version of the */
/*     periodic QR algorithm described in [1,3] with some minor */
/*     modifications [2]. The eigenvalues are computed via an implicit */
/*     complex single shift algorithm. */

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

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPQR). */

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

/*     Check the scalar input parameters. */

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
    *info = 0;
    nh = *ihi - *ilo + 1;
    nq = *ihiq - *iloq + 1;
    if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -5;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -6;
    } else if (*iloq < 1 || *iloq > *ilo) {
	*info = -7;
    } else if (*ihiq < *ihi || *ihiq > *n) {
	*info = -8;
    } else if (*lda < max(1,*n)) {
	*info = -10;
    } else if (*ldb < max(1,*n)) {
	*info = -12;
    } else if (*ldq < 1 || *wantq && *ldq < *n) {
	*info = -14;
    } else if (*ldz < 1 || *wantz && *ldz < *n) {
	*info = -16;
    } else if (*ldwork < max(1,*n)) {
	dwork[1] = (doublereal) max(1,*n);
	*info = -21;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03YD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     Set machine-dependent constants for the stopping criterion. */

    unfl = dlamch_("Safe minimum", (ftnlen)12);
    ovfl = 1. / unfl;
    dlabad_(&unfl, &ovfl);
    ulp = dlamch_("Precision", (ftnlen)9);
    smlnum = unfl * (nh / ulp);

/*     I1 and I2 are the indices of the first rows and last columns of */
/*     A and B to which transformations must be applied. */

    i1 = 1;
    i2 = *n;
    iseed[0] = 1;
    iseed[1] = 0;
    iseed[2] = 0;
    iseed[3] = 1;

/*     ITN is the maximal number of QR iterations. */

    itn = nh * 30;

/*     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO */
/*     or A(L,L-1) is negligible. */

    i__ = *ihi;
L10:
    l = *ilo;
    if (i__ < *ilo) {
	goto L120;
    }

/*     Perform periodic QR iteration on rows and columns ILO to I of A */
/*     and B until a submatrix of order 1 or 2 splits off at the bottom. */

    i__1 = itn;
    for (its = 0; its <= i__1; ++its) {

/*        Look for deflations in A. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {
	    tst = (d__1 = a[k - 1 + (k - 1) * a_dim1], abs(d__1)) + (d__2 = a[
		    k + k * a_dim1], abs(d__2));
	    if (tst == 0.) {
		i__3 = i__ - l + 1;
		tst = dlanhs_("1", &i__3, &a[l + l * a_dim1], lda, &dwork[1], 
			(ftnlen)1);
	    }
/* Computing MAX */
	    d__2 = ulp * tst;
	    if ((d__1 = a[k + (k - 1) * a_dim1], abs(d__1)) <= max(d__2,
		    smlnum)) {
		goto L30;
	    }
/* L20: */
	}
L30:

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
		    i__3 = i__ - k + 1;
		    tst = dlanhs_("1", &i__3, &b[k + k * b_dim1], ldb, &dwork[
			    1], (ftnlen)1);
		}
/* Computing MAX */
		d__2 = ulp * tst;
		if ((d__1 = b[kk + kk * b_dim1], abs(d__1)) <= max(d__2,
			smlnum)) {
		    goto L50;
		}
/* L40: */
	    }
	} else {
	    kk = k - 1;
	}
L50:
	if (kk >= k) {

/*           B has an element close to zero at position (KK,KK). */

	    b[kk + kk * b_dim1] = 0.;
	    mb03ya_(wantt, wantq, wantz, n, &k, &i__, iloq, ihiq, &kk, &a[
		    a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
		    z__[z_offset], ldz, info);
	    k = kk + 1;
	}
	l = k;
	if (l > *ilo) {

/*           A(L,L-1) is negligible. */

	    a[l + (l - 1) * a_dim1] = 0.;
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

	if (l >= i__ - 1) {
	    goto L80;
	}

/*        The active submatrices are now in rows and columns L:I. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i__;
	}
	if (its == 10 || its == 20) {

/*           Exceptional shift. The first column of the shift polynomial */
/*           is a pseudo-random vector. */

	    dlarnv_(&c__3, iseed, &c__3, v);
	} else {

/*           The implicit double shift is constructed via a partial */
/*           product QR factorization [2]. */

	    dlartg_(&b[l + l * b_dim1], &b[i__ + i__ * b_dim1], &cs2, &sn2, &
		    temp);
	    dlartg_(&temp, &b[i__ - 1 + i__ * b_dim1], &cs1, &sn1, &alpha);

	    alpha = a[l + l * a_dim1] * cs2 - a[i__ + i__ * a_dim1] * sn2;
	    betax = cs1 * (cs2 * a[l + 1 + l * a_dim1]);
	    gamma = cs1 * (sn2 * a[i__ - 1 + i__ * a_dim1]) + sn1 * a[i__ - 1 
		    + (i__ - 1) * a_dim1];
	    alpha = alpha * cs1 - a[i__ + (i__ - 1) * a_dim1] * sn1;
	    dlartg_(&alpha, &betax, &cs1, &sn1, &temp);

	    dlartg_(&temp, &gamma, &cs2, &sn2, &alpha);
	    alpha = cs2;
	    gamma = a[i__ - 1 + (i__ - 1) * a_dim1] * cs1 * cs2 + a[i__ + (
		    i__ - 1) * a_dim1] * sn2;
	    delta = a[i__ - 1 + (i__ - 1) * a_dim1] * sn1 * cs2;
	    dlartg_(&gamma, &delta, &cs3, &sn3, &temp);
	    dlartg_(&alpha, &temp, &cs2, &sn2, &alpha);

	    alpha = (b[l + l * b_dim1] * cs1 + b[l + (l + 1) * b_dim1] * sn1) 
		    * cs2;
	    betax = b[l + 1 + (l + 1) * b_dim1] * sn1 * cs2;
	    gamma = b[i__ - 1 + (i__ - 1) * b_dim1] * sn2;
	    dlartg_(&alpha, &betax, &cs1, &sn1, &temp);
	    dlartg_(&temp, &gamma, &cs2, &sn2, &alpha);

	    alpha = cs1 * a[l + l * a_dim1] + sn1 * a[l + (l + 1) * a_dim1];
	    betax = cs1 * a[l + 1 + l * a_dim1] + sn1 * a[l + 1 + (l + 1) * 
		    a_dim1];
	    gamma = sn1 * a[l + 2 + (l + 1) * a_dim1];

	    v[0] = cs2 * alpha - sn2 * cs3;
	    v[1] = cs2 * betax - sn2 * sn3;
	    v[2] = gamma * cs2;
	}

/*        Double-shift QR step */

	i__2 = i__ - 1;
	for (k = l; k <= i__2; ++k) {

/* Computing MIN */
	    i__3 = 3, i__4 = i__ - k + 1;
	    nr = min(i__3,i__4);
	    if (k > l) {
		dcopy_(&nr, &a[k + (k - 1) * a_dim1], &c__1, v, &c__1);
	    }
	    dlarfg_(&nr, v, &v[1], &c__1, &tauv);
	    if (k > l) {
		a[k + (k - 1) * a_dim1] = v[0];
		a[k + 1 + (k - 1) * a_dim1] = 0.;
		if (k < i__ - 1) {
		    a[k + 2 + (k - 1) * a_dim1] = 0.;
		}
	    }

/*           Apply reflector V from the right to B in rows I1:min(K+2,I). */

	    v[0] = 1.;
/* Computing MIN */
	    i__4 = k + 2;
	    i__3 = min(i__4,i__) - i1 + 1;
	    dlarfx_("Right", &i__3, &nr, v, &tauv, &b[i1 + k * b_dim1], ldb, &
		    dwork[1], (ftnlen)5);

/*           Annihilate the introduced nonzeros in the K-th column. */

	    dcopy_(&nr, &b[k + k * b_dim1], &c__1, w, &c__1);
	    dlarfg_(&nr, w, &w[1], &c__1, &tauw);
	    b[k + k * b_dim1] = w[0];
	    b[k + 1 + k * b_dim1] = 0.;
	    if (k < i__ - 1) {
		b[k + 2 + k * b_dim1] = 0.;
	    }

/*           Apply reflector W from the left to transform the rows of the */
/*           matrix B in columns K+1:I2. */

	    w[0] = 1.;
	    i__3 = i2 - k;
	    dlarfx_("Left", &nr, &i__3, w, &tauw, &b[k + (k + 1) * b_dim1], 
		    ldb, &dwork[1], (ftnlen)4);

/*           Apply reflector V from the left to transform the rows of the */
/*           matrix A in columns K:I2. */

	    i__3 = i2 - k + 1;
	    dlarfx_("Left", &nr, &i__3, v, &tauv, &a[k + k * a_dim1], lda, &
		    dwork[1], (ftnlen)4);

/*           Apply reflector W from the right to transform the columns of */
/*           the matrix A in rows I1:min(K+3,I). */

/* Computing MIN */
	    i__4 = k + 3;
	    i__3 = min(i__4,i__) - i1 + 1;
	    dlarfx_("Right", &i__3, &nr, w, &tauw, &a[i1 + k * a_dim1], lda, &
		    dwork[1], (ftnlen)5);

/*           Accumulate transformations in the matrices Q and Z. */

	    if (*wantq) {
		dlarfx_("Right", &nq, &nr, v, &tauv, &q[*iloq + k * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
	    }
	    if (*wantz) {
		dlarfx_("Right", &nq, &nr, w, &tauw, &z__[*iloq + k * z_dim1],
			 ldz, &dwork[1], (ftnlen)5);
	    }
/* L60: */
	}
/* L70: */
    }

/*     Failure to converge. */

    *info = i__;
    return 0;

L80:

/*     Compute 1-by-1 or 2-by-2 subproblem. */

    if (l == i__) {

/*        Standardize B, set ALPHAR, ALPHAI and BETA. */

	if (b[i__ + i__ * b_dim1] < 0.) {
	    if (*wantt) {
		i__1 = i__;
		for (k = i1; k <= i__1; ++k) {
		    b[k + i__ * b_dim1] = -b[k + i__ * b_dim1];
/* L90: */
		}
		i__1 = i2;
		for (k = i__; k <= i__1; ++k) {
		    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1];
/* L100: */
		}
	    } else {
		b[i__ + i__ * b_dim1] = -b[i__ + i__ * b_dim1];
		a[i__ + i__ * a_dim1] = -a[i__ + i__ * a_dim1];
	    }
	    if (*wantq) {
		i__1 = *ihiq;
		for (k = *iloq; k <= i__1; ++k) {
		    q[k + i__ * q_dim1] = -q[k + i__ * q_dim1];
/* L110: */
		}
	    }
	}
	alphar[i__] = a[i__ + i__ * a_dim1];
	alphai[i__] = 0.;
	beta[i__] = b[i__ + i__ * b_dim1];
    } else if (l == i__ - 1) {

/*        A double block has converged. */
/*        Compute eigenvalues and standardize double block. */

	mb03yt_(&a[i__ - 1 + (i__ - 1) * a_dim1], lda, &b[i__ - 1 + (i__ - 1) 
		* b_dim1], ldb, &alphar[i__ - 1], &alphai[i__ - 1], &beta[i__ 
		- 1], &cs1, &sn1, &cs2, &sn2);

/*        Apply transformation to rest of A and B. */

	if (i2 > i__) {
	    i__1 = i2 - i__;
	    drot_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], lda, &a[i__ + (i__ 
		    + 1) * a_dim1], lda, &cs1, &sn1);
	}
	i__1 = i__ - i1 - 1;
	drot_(&i__1, &a[i1 + (i__ - 1) * a_dim1], &c__1, &a[i1 + i__ * a_dim1]
		, &c__1, &cs2, &sn2);
	if (i2 > i__) {
	    i__1 = i2 - i__;
	    drot_(&i__1, &b[i__ - 1 + (i__ + 1) * b_dim1], ldb, &b[i__ + (i__ 
		    + 1) * b_dim1], ldb, &cs2, &sn2);
	}
	i__1 = i__ - i1 - 1;
	drot_(&i__1, &b[i1 + (i__ - 1) * b_dim1], &c__1, &b[i1 + i__ * b_dim1]
		, &c__1, &cs1, &sn1);

/*        Apply transformation to rest of Q and Z if desired. */

	if (*wantq) {
	    drot_(&nq, &q[*iloq + (i__ - 1) * q_dim1], &c__1, &q[*iloq + i__ *
		     q_dim1], &c__1, &cs1, &sn1);
	}
	if (*wantz) {
	    drot_(&nq, &z__[*iloq + (i__ - 1) * z_dim1], &c__1, &z__[*iloq + 
		    i__ * z_dim1], &c__1, &cs2, &sn2);
	}
    }

/*     Decrement number of remaining iterations, and return to start of */
/*     the main loop with new value of I. */

    itn -= its;
    i__ = l - 1;
    goto L10;

L120:
    dwork[1] = (doublereal) max(1,*n);
    return 0;
/* *** Last line of MB03YD *** */
} /* mb03yd_ */

