/* AB13MD.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b15 = 0.;
static doublereal c_b136 = 2.;

/* Subroutine */ int ab13md_(char *fact, integer *n, doublecomplex *z__, 
	integer *ldz, integer *m, integer *nblock, integer *itype, doublereal 
	*x, doublereal *bound, doublereal *d__, doublereal *g, integer *iwork,
	 doublereal *dwork, integer *ldwork, doublecomplex *zwork, integer *
	lzwork, integer *info, ftnlen fact_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double pow_di(doublereal *, integer *), log(doublereal);

    /* Local variables */
    static doublereal c__, e;
    static integer i__, j, k, l;
    static doublereal t1, t2, t3, hn;
    static integer mr;
    static doublereal pp;
    static integer mt, iw2, iw3, iw4, iw5, iw6, iw7, iw8, iw9, iz2, iz3, iz4, 
	    iz5, iz6, iz7, iz8, iz9, iw10, iw11, iw12, iw13, iw14, iw15, iw16,
	     iw17, iw18, iw19, iw20, iw21, iw22, iw23, iw24, iw25, iw26, iw27,
	     iw28, iw29, iw30, iw31, iw32, iw33, iz10, iz11, iz12, iz13, iz14,
	     iz15, iz16, iz17, iz18, iz19, iz20, iz21, iz22, iz23, iz24, lwa, 
	    lza;
    static doublereal eps, phi, rat, tau, tol;
    static logical pos;
    static doublereal tol2, tol3, tol4, tol5;
    static doublecomplex detf;
    static doublereal emin, emax;
    static integer sdim;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer iter;
    static doublereal prod, temp;
    static integer iwrk, isum, nsum, info2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale, delta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    static logical xfact;
    extern /* Subroutine */ int zgees_(char *, char *, L_fp, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, logical *, integer *, ftnlen, ftnlen), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    zgemm_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static doublereal hnorm, svlam;
    static logical bwork[1], gtest;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal snorm;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal znorm;
    static integer izwrk;
    extern /* Subroutine */ int dsysv_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static doublereal hnorm1, dlambd, ynorm1, ynorm2, znorm2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static doublereal regpar;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublecomplex tempij;
    static integer lwamax;
    static doublecomplex tempji;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer lzamax;
    extern /* Subroutine */ int dsycon_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen), zgetrf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *);
    static doublereal colsum;
    extern /* Subroutine */ int zgesvd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen), zgetri_(integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *), zlacpy_(char *, integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, ftnlen);
    static integer minwrk, minzrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static doublereal stsize;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal rowsum;


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

/*     To compute an upper bound on the structured singular value for a */
/*     given square complex matrix and a given block structure of the */
/*     uncertainty. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not an information from the */
/*             previous call is supplied in the vector X. */
/*             = 'F':  On entry, X contains information from the */
/*                     previous call. */
/*             = 'N':  On entry, X does not contain an information from */
/*                     the previous call. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix Z.  N >= 0. */

/*     Z       (input) COMPLEX*16 array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array must contain the */
/*             complex matrix Z for which the upper bound on the */
/*             structured singular value is to be computed. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= max(1,N). */

/*     M       (input) INTEGER */
/*             The number of diagonal blocks in the block structure of */
/*             the uncertainty.  M >= 1. */

/*     NBLOCK  (input) INTEGER array, dimension (M) */
/*             The vector of length M containing the block structure */
/*             of the uncertainty. NBLOCK(I), I = 1:M, is the size of */
/*             each block. */

/*     ITYPE   (input) INTEGER array, dimension (M) */
/*             The vector of length M indicating the type of each block. */
/*             For I = 1:M, */
/*             ITYPE(I) = 1 indicates that the corresponding block is a */
/*                          real block, and */
/*             ITYPE(I) = 2 indicates that the corresponding block is a */
/*                          complex block. */
/*             NBLOCK(I) must be equal to 1 if ITYPE(I) is equal to 1. */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             ( M + MR - 1 ), where MR is the number of the real blocks. */
/*             On entry, if FACT = 'F' and NBLOCK(1) < N, this array */
/*             must contain information from the previous call to AB13MD. */
/*             If NBLOCK(1) = N, this array is not used. */
/*             On exit, if NBLOCK(1) < N, this array contains information */
/*             that can be used in the next call to AB13MD for a matrix */
/*             close to Z. */

/*     BOUND   (output) DOUBLE PRECISION */
/*             The upper bound on the structured singular value. */

/*     D, G    (output) DOUBLE PRECISION arrays, dimension (N) */
/*             The vectors of length N containing the diagonal entries */
/*             of the diagonal N-by-N matrices D and G, respectively, */
/*             such that the matrix */
/*             Z'*D^2*Z + sqrt(-1)*(G*Z-Z'*G) - BOUND^2*D^2 */
/*             is negative semidefinite. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(4*M-2,N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 2*N*N*M - N*N + 9*M*M + N*M + 11*N + 33*M - 11. */
/*             For best performance */
/*             LDWORK >= 2*N*N*M - N*N + 9*M*M + N*M + 6*N + 33*M - 11 + */
/*                       MAX( 5*N,2*N*NB ) */
/*             where NB is the optimal blocksize returned by ILAENV. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) contains the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The dimension of the array ZWORK. */
/*             LZWORK >= 6*N*N*M + 12*N*N + 6*M + 6*N - 3. */
/*             For best performance */
/*             LZWORK >= 6*N*N*M + 12*N*N + 6*M + 3*N - 3 + */
/*                       MAX( 3*N,N*NB ) */
/*             where NB is the optimal blocksize returned by ILAENV. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the block sizes must be positive integers; */
/*             = 2:  the sum of block sizes must be equal to N; */
/*             = 3:  the size of a real block must be equal to 1; */
/*             = 4:  the block type must be either 1 or 2; */
/*             = 5:  errors in solving linear equations or in matrix */
/*                   inversion; */
/*             = 6:  errors in computing eigenvalues or singular values. */

/*     METHOD */

/*     The routine computes the upper bound proposed in [1]. */

/*     REFERENCES */

/*     [1] Fan, M.K.H., Tits, A.L., and Doyle, J.C. */
/*         Robustness in the presence of mixed parametric uncertainty */
/*         and unmodeled dynamics. */
/*         IEEE Trans. Automatic Control, vol. AC-36, 1991, pp. 25-38. */

/*     NUMERICAL ASPECTS */

/*     The accuracy and speed of computation depend on the value of */
/*     the internal threshold TOL. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, F. Delebecque, D.W. Gu, M.M. Konstantinov and */
/*     S. Steer with the assistance of V. Sima, September 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Universiteit Leuven, February 2001. */

/*     KEYWORDS */

/*     H-infinity optimal control, Robust control, Structured singular */
/*     value. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute workspace. */

    /* Parameter adjustments */
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --nblock;
    --itype;
    --x;
    --d__;
    --g;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    minwrk = (*n << 1) * *n * *m - *n * *n + *m * 9 * *m + *n * *m + *n * 11 
	    + *m * 33 - 11;
    minzrk = *n * 6 * *n * *m + *n * 12 * *n + *m * 6 + *n * 6 - 3;

/*     Decode and Test input parameters. */

    *info = 0;
    xfact = lsame_(fact, "F", (ftnlen)1, (ftnlen)1);
    if (! xfact && ! lsame_(fact, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < max(1,*n)) {
	*info = -4;
    } else if (*m < 1) {
	*info = -5;
    } else if (*ldwork < minwrk) {
	*info = -14;
    } else if (*lzwork < minzrk) {
	*info = -16;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("AB13MD", &i__1, (ftnlen)6);
	return 0;
    }

    nsum = 0;
    isum = 0;
    mr = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nblock[i__] < 1) {
	    *info = 1;
	    return 0;
	}
	if (itype[i__] == 1 && nblock[i__] > 1) {
	    *info = 3;
	    return 0;
	}
	nsum += nblock[i__];
	if (itype[i__] == 1) {
	    ++mr;
	}
	if (itype[i__] == 1 || itype[i__] == 2) {
	    ++isum;
	}
/* L10: */
    }
    if (nsum != *n) {
	*info = 2;
	return 0;
    }
    if (isum != *m) {
	*info = 4;
	return 0;
    }
    mt = *m + mr - 1;

    lwamax = 0;
    lzamax = 0;

/*     Set D = In, G = 0. */

    dlaset_("Full", n, &c__1, &c_b11, &c_b11, &d__[1], n, (ftnlen)4);
    dlaset_("Full", n, &c__1, &c_b15, &c_b15, &g[1], n, (ftnlen)4);

/*     Quick return if possible. */

    znorm = zlange_("F", n, n, &z__[z_offset], ldz, &dwork[1], (ftnlen)1);
    if (znorm == 0.) {
	*bound = 0.;
	dwork[1] = 1.;
	zwork[1].r = 1., zwork[1].i = 0.;
	return 0;
    }

/*     Copy Z into ZWORK. */

    zlacpy_("Full", n, n, &z__[z_offset], ldz, &zwork[1], n, (ftnlen)4);

/*     Exact bound for the case NBLOCK( 1 ) = N. */

    if (nblock[1] == *n) {
	if (itype[1] == 1) {

/*           1-by-1 real block. */

	    *bound = 0.;
	    dwork[1] = 1.;
	    zwork[1].r = 1., zwork[1].i = 0.;
	} else {

/*           N-by-N complex block. */

	    zgesvd_("N", "N", n, n, &zwork[1], n, &dwork[1], &zwork[1], &c__1,
		     &zwork[1], &c__1, &zwork[*n * *n + 1], lzwork, &dwork[*n 
		    + 1], &info2, (ftnlen)1, (ftnlen)1);
	    if (info2 > 0) {
		*info = 6;
		return 0;
	    }
	    *bound = dwork[1];
	    i__1 = *n * *n + 1;
	    lza = *n * *n + (integer) zwork[i__1].r;
	    dwork[1] = (doublereal) (*n * 5);
	    z__1.r = (doublereal) lza, z__1.i = 0.;
	    zwork[1].r = z__1.r, zwork[1].i = z__1.i;
	}
	return 0;
    }

/*     Get machine precision. */

    eps = dlamch_("P", (ftnlen)1);

/*     Set tolerances. */

    tol = sqrt(eps) * 100.;
    tol2 = eps * 1e4;
    tol3 = eps * 10.;
    tol4 = .001;
    tol5 = .001;
    regpar = eps * 1e3;

/*     Real workspace usage. */

    iw2 = *m * *m;
    iw3 = iw2 + *m;
    iw4 = iw3 + *n;
    iw5 = iw4 + *m;
    iw6 = iw5 + *m;
    iw7 = iw6 + *n;
    iw8 = iw7 + *n;
    iw9 = iw8 + *n * (*m - 1);
    iw10 = iw9 + *n * *n * mt;
    iw11 = iw10 + mt;
    iw12 = iw11 + mt * mt;
    iw13 = iw12 + *n;
    iw14 = iw13 + mt + 1;
    iw15 = iw14 + mt + 1;
    iw16 = iw15 + mt + 1;
    iw17 = iw16 + mt + 1;
    iw18 = iw17 + mt + 1;
    iw19 = iw18 + mt;
    iw20 = iw19 + mt;
    iw21 = iw20 + mt;
    iw22 = iw21 + *n;
    iw23 = iw22 + *m - 1;
    iw24 = iw23 + mr;
    iw25 = iw24 + *n;
    iw26 = iw25 + (mt << 1);
    iw27 = iw26 + mt;
    iw28 = iw27 + mt;
    iw29 = iw28 + *m - 1;
    iw30 = iw29 + mr;
    iw31 = iw30 + *n + (mt << 1);
    iw32 = iw31 + mt * mt;
    iw33 = iw32 + mt;
    iwrk = iw33 + mt + 1;

/*     Double complex workspace usage. */

    iz2 = *n * *n;
    iz3 = iz2 + *n * *n;
    iz4 = iz3 + *n * *n;
    iz5 = iz4 + *n * *n;
    iz6 = iz5 + *n * *n;
    iz7 = iz6 + *n * *n * mt;
    iz8 = iz7 + *n * *n;
    iz9 = iz8 + *n * *n;
    iz10 = iz9 + *n * *n;
    iz11 = iz10 + mt;
    iz12 = iz11 + *n * *n;
    iz13 = iz12 + *n;
    iz14 = iz13 + *n * *n;
    iz15 = iz14 + *n;
    iz16 = iz15 + *n * *n;
    iz17 = iz16 + *n;
    iz18 = iz17 + *n * *n;
    iz19 = iz18 + *n * *n * mt;
    iz20 = iz19 + mt;
    iz21 = iz20 + *n * *n * mt;
    iz22 = iz21 + *n * *n;
    iz23 = iz22 + *n * *n;
    iz24 = iz23 + *n * *n;
    izwrk = iz24 + mt;

/*     Compute the cumulative sums of blocks dimensions. */

    iwork[1] = 0;
    i__1 = *m + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	iwork[i__] = iwork[i__ - 1] + nblock[i__ - 1];
/* L20: */
    }

/*     Find Osborne scaling if initial scaling is not given. */

    if (! xfact) {
	dlaset_("Full", m, m, &c_b15, &c_b15, &dwork[1], m, (ftnlen)4);
	dlaset_("Full", m, &c__1, &c_b11, &c_b11, &dwork[iw2 + 1], m, (ftnlen)
		4);
	znorm = zlange_("F", n, n, &zwork[1], n, &dwork[1], (ftnlen)1);
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ != j) {
		    i__3 = iwork[i__ + 1] - iwork[i__];
		    i__4 = iwork[j + 1] - iwork[j];
		    zlacpy_("Full", &i__3, &i__4, &z__[iwork[i__] + 1 + (
			    iwork[j] + 1) * z_dim1], ldz, &zwork[iz2 + 1], n, 
			    (ftnlen)4);
		    i__3 = iwork[i__ + 1] - iwork[i__];
		    i__4 = iwork[j + 1] - iwork[j];
		    i__5 = *lzwork - izwrk;
		    zgesvd_("N", "N", &i__3, &i__4, &zwork[iz2 + 1], n, &
			    dwork[iw3 + 1], &zwork[1], &c__1, &zwork[1], &
			    c__1, &zwork[izwrk + 1], &i__5, &dwork[iwrk + 1], 
			    &info2, (ftnlen)1, (ftnlen)1);
		    if (info2 > 0) {
			*info = 6;
			return 0;
		    }
		    i__3 = izwrk + 1;
		    lza = (integer) zwork[i__3].r;
		    lzamax = max(lza,lzamax);
		    znorm2 = dwork[iw3 + 1];
		    dwork[i__ + (j - 1) * *m] = znorm2 + znorm * tol2;
		}
/* L30: */
	    }
/* L40: */
	}
	dlaset_("Full", m, &c__1, &c_b15, &c_b15, &dwork[iw4 + 1], m, (ftnlen)
		4);
L50:
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw5 + i__] = dwork[iw4 + i__] - 1.;
/* L60: */
	}
	hnorm = dlange_("F", m, &c__1, &dwork[iw5 + 1], m, &dwork[1], (ftnlen)
		1);
	if (hnorm <= tol2) {
	    goto L120;
	}
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    colsum = 0.;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		colsum += dwork[i__ + (k - 1) * *m];
/* L70: */
	    }
	    rowsum = 0.;
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		rowsum += dwork[k + (j - 1) * *m];
/* L80: */
	    }
	    rat = sqrt(colsum / rowsum);
	    dwork[iw4 + k] = rat;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[i__ + (k - 1) * *m] /= rat;
/* L90: */
	    }
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		dwork[k + (j - 1) * *m] *= rat;
/* L100: */
	    }
	    dwork[iw2 + k] *= rat;
/* L110: */
	}
	goto L50;
L120:
	scale = 1. / dwork[iw2 + 1];
	dscal_(m, &scale, &dwork[iw2 + 1], &c__1);
    } else {
	dwork[iw2 + 1] = 1.;
	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    dwork[iw2 + i__] = sqrt(x[i__ - 1]);
/* L130: */
	}
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ != j) {
		i__3 = iwork[i__ + 1] - iwork[i__];
		i__4 = iwork[j + 1] - iwork[j];
		zlascl_("G", m, m, &dwork[iw2 + j], &dwork[iw2 + i__], &i__3, 
			&i__4, &zwork[iwork[i__] + 1 + iwork[j] * *n], n, &
			info2, (ftnlen)1);
	    }
/* L140: */
	}
/* L150: */
    }

/*     Scale Z by its 2-norm. */

    zlacpy_("Full", n, n, &zwork[1], n, &zwork[iz2 + 1], n, (ftnlen)4);
    i__1 = *lzwork - izwrk;
    zgesvd_("N", "N", n, n, &zwork[iz2 + 1], n, &dwork[iw3 + 1], &zwork[1], &
	    c__1, &zwork[1], &c__1, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 1]
	    , &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    znorm = dwork[iw3 + 1];
    zlascl_("G", m, m, &znorm, &c_b11, n, n, &zwork[1], n, &info2, (ftnlen)1);

/*     Set BB. */

    i__1 = *n * *n;
    i__2 = *n * *n;
    dlaset_("Full", &i__1, &mt, &c_b15, &c_b15, &dwork[iw9 + 1], &i__2, (
	    ftnlen)4);

/*     Set P. */

    i__1 = nblock[1];
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw6 + i__] = 1.;
/* L160: */
    }
    i__1 = *n;
    for (i__ = nblock[1] + 1; i__ <= i__1; ++i__) {
	dwork[iw6 + i__] = 0.;
/* L170: */
    }

/*     Compute P*Z. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz3 + i__ + (j - 1) * *n;
	    i__4 = iw6 + i__;
	    z__2.r = dwork[i__4], z__2.i = 0.;
	    i__5 = i__ + (j - 1) * *n;
	    z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, z__1.i =
		     z__2.r * zwork[i__5].i + z__2.i * zwork[i__5].r;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L180: */
	}
/* L190: */
    }

/*     Compute Z'*P*Z. */

    zgemm_("C", "N", n, n, n, &c_b2, &zwork[1], n, &zwork[iz3 + 1], n, &c_b1, 
	    &zwork[iz4 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Copy Z'*P*Z into A0. */

    zlacpy_("Full", n, n, &zwork[iz4 + 1], n, &zwork[iz5 + 1], n, (ftnlen)4);

/*     Copy diag(P) into B0d. */

    dcopy_(n, &dwork[iw6 + 1], &c__1, &dwork[iw7 + 1], &c__1);

    i__1 = *m;
    for (k = 2; k <= i__1; ++k) {

/*        Set P. */

	i__2 = iwork[k];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw6 + i__] = 0.;
/* L200: */
	}
	i__2 = iwork[k] + nblock[k];
	for (i__ = iwork[k] + 1; i__ <= i__2; ++i__) {
	    dwork[iw6 + i__] = 1.;
/* L210: */
	}
	if (k < *m) {
	    i__2 = *n;
	    for (i__ = iwork[k + 1] + 1; i__ <= i__2; ++i__) {
		dwork[iw6 + i__] = 0.;
/* L220: */
	    }
	}

/*        Compute P*Z. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = iz3 + i__ + (j - 1) * *n;
		i__5 = iw6 + i__;
		z__2.r = dwork[i__5], z__2.i = 0.;
		i__6 = i__ + (j - 1) * *n;
		z__1.r = z__2.r * zwork[i__6].r - z__2.i * zwork[i__6].i, 
			z__1.i = z__2.r * zwork[i__6].i + z__2.i * zwork[i__6]
			.r;
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L230: */
	    }
/* L240: */
	}

/*        Compute t = Z'*P*Z. */

	zgemm_("C", "N", n, n, n, &c_b2, &zwork[1], n, &zwork[iz3 + 1], n, &
		c_b1, &zwork[iz4 + 1], n, (ftnlen)1, (ftnlen)1);

/*        Copy t(:) into the (k-1)-th column of AA. */

	i__2 = *n * *n;
	zcopy_(&i__2, &zwork[iz4 + 1], &c__1, &zwork[iz6 + 1 + (k - 2) * *n * 
		*n], &c__1);

/*        Copy diag(P) into the (k-1)-th column of BBd. */

	dcopy_(n, &dwork[iw6 + 1], &c__1, &dwork[iw8 + 1 + (k - 2) * *n], &
		c__1);

/*        Copy P(:) into the (k-1)-th column of BB. */

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw9 + i__ + (i__ - 1) * *n + (k - 2) * *n * *n] = dwork[iw6 
		    + i__];
/* L260: */
	}
/* L270: */
    }

    l = 0;

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	if (itype[k] == 1) {
	    ++l;

/*           Set P. */

	    i__2 = iwork[k];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[iw6 + i__] = 0.;
/* L280: */
	    }
	    i__2 = iwork[k] + nblock[k];
	    for (i__ = iwork[k] + 1; i__ <= i__2; ++i__) {
		dwork[iw6 + i__] = 1.;
/* L290: */
	    }
	    if (k < *m) {
		i__2 = *n;
		for (i__ = iwork[k + 1] + 1; i__ <= i__2; ++i__) {
		    dwork[iw6 + i__] = 0.;
/* L300: */
		}
	    }

/*           Compute P*Z. */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = iz3 + i__ + (j - 1) * *n;
		    i__5 = iw6 + i__;
		    z__2.r = dwork[i__5], z__2.i = 0.;
		    i__6 = i__ + (j - 1) * *n;
		    z__1.r = z__2.r * zwork[i__6].r - z__2.i * zwork[i__6].i, 
			    z__1.i = z__2.r * zwork[i__6].i + z__2.i * zwork[
			    i__6].r;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L310: */
		}
/* L320: */
	    }

/*           Compute t = sqrt(-1)*( P*Z - Z'*P ). */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = j;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = iz3 + i__ + (j - 1) * *n;
		    tempij.r = zwork[i__4].r, tempij.i = zwork[i__4].i;
		    i__4 = iz3 + j + (i__ - 1) * *n;
		    tempji.r = zwork[i__4].r, tempji.i = zwork[i__4].i;
		    i__4 = iz4 + i__ + (j - 1) * *n;
		    d_cnjg(&z__3, &tempji);
		    z__2.r = tempij.r - z__3.r, z__2.i = tempij.i - z__3.i;
		    z__1.r = z__2.r * 0. - z__2.i * 1., z__1.i = z__2.i * 0. 
			    + z__2.r * 1.;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		    i__4 = iz4 + j + (i__ - 1) * *n;
		    d_cnjg(&z__3, &tempij);
		    z__2.r = tempji.r - z__3.r, z__2.i = tempji.i - z__3.i;
		    z__1.r = z__2.r * 0. - z__2.i * 1., z__1.i = z__2.i * 0. 
			    + z__2.r * 1.;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L330: */
		}
/* L340: */
	    }

/*           Copy t(:) into the (m-1+l)-th column of AA. */

	    i__2 = *n * *n;
	    zcopy_(&i__2, &zwork[iz4 + 1], &c__1, &zwork[iz6 + 1 + (*m - 2 + 
		    l) * *n * *n], &c__1);
	}
/* L350: */
    }

/*     Set initial X. */

    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 1.;
/* L360: */
    }
    if (mr > 0) {
	if (! xfact) {
	    i__1 = mr;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[*m - 1 + i__] = 0.;
/* L370: */
	    }
	} else {
	    l = 0;
	    i__1 = *m;
	    for (k = 1; k <= i__1; ++k) {
		if (itype[k] == 1) {
		    ++l;
/* Computing 2nd power */
		    d__1 = dwork[iw2 + k];
		    x[*m - 1 + l] /= d__1 * d__1;
		}
/* L380: */
	    }
	}
    }

/*     Set constants. */

    svlam = 1. / eps;
    c__ = 1.;

/*     Set H. */

    dlaset_("Full", &mt, &mt, &c_b15, &c_b11, &dwork[iw11 + 1], &mt, (ftnlen)
	    4);

    iter = -1;

/*     Main iteration loop. */

L390:
    ++iter;

/*        Compute A(:) = A0 + AA*x. */

    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz10 + i__;
	i__3 = i__;
	z__1.r = x[i__3], z__1.i = 0.;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L400: */
    }
    i__1 = *n * *n;
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute diag( Binv ). */

    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw12 + 1], &c__1);
    i__1 = *m - 1;
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &x[1], &c__1, &c_b11, &
	    dwork[iw12 + 1], &c__1, (ftnlen)1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw12 + i__] = 1. / dwork[iw12 + i__];
/* L410: */
    }

/*        Compute Binv*A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz11 + i__ + (j - 1) * *n;
	    i__4 = iw12 + i__;
	    z__2.r = dwork[i__4], z__2.i = 0.;
	    i__5 = iz7 + i__ + (j - 1) * *n;
	    z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, z__1.i =
		     z__2.r * zwork[i__5].i + z__2.i * zwork[i__5].r;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L420: */
	}
/* L430: */
    }

/*        Compute eig( Binv*A ). */

    i__1 = *lzwork - izwrk;
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz11 + 1], n, &sdim, &zwork[
	    iz12 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    i__1 = iz12 + 1;
    e = zwork[i__1].r;
    if (*n > 1) {
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = iz12 + i__;
	    if (zwork[i__2].r > e) {
		i__3 = iz12 + i__;
		e = zwork[i__3].r;
	    }
/* L440: */
	}
    }

/*        Set tau. */

    if (mr > 0) {
	snorm = (d__1 = x[*m], abs(d__1));
	if (mr > 1) {
	    i__1 = mt;
	    for (i__ = *m + 1; i__ <= i__1; ++i__) {
		if ((d__1 = x[i__], abs(d__1)) > snorm) {
		    snorm = (d__2 = x[i__], abs(d__2));
		}
/* L450: */
	    }
	}
	if (snorm > 40.) {
	    tau = 100.;
	} else if (snorm > 8.) {
	    tau = 50.;
	} else if (snorm > 4.) {
	    tau = 10.;
	} else if (snorm > 1.) {
	    tau = 5.;
	} else {
	    tau = 2.;
	}
    }
    if (iter == 0) {
	dlambd = e + .001;
    } else {
	dwork[iw13 + 1] = e;
	dcopy_(&mt, &x[1], &c__1, &dwork[iw13 + 2], &c__1);
	dlambd = dwork[iw13 + 1] * .98999999999999999 + dwork[iw14 + 1] * .01;
	dcopy_(&mt, &dwork[iw13 + 2], &c__1, &dwork[iw18 + 1], &c__1);
	dcopy_(&mt, &dwork[iw14 + 2], &c__1, &dwork[iw19 + 1], &c__1);
	l = 0;
L460:
	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = (1. - .01 / pow_di(&c_b136, &l)) * dwork[iw18 + i__] + 
		    .01 / pow_di(&c_b136, &l) * dwork[iw19 + i__];
/* L470: */
	}

/*           Compute At(:) = A0 + AA*x. */

	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = iz10 + i__;
	    i__3 = i__;
	    z__1.r = x[i__3], z__1.i = 0.;
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L480: */
	}
	i__1 = *n * *n;
	zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz9 + 1], &c__1);
	i__1 = *n * *n;
	i__2 = *n * *n;
	zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz9 + 1], &c__1, (ftnlen)1);

/*           Compute diag(Bt). */

	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw21 + 1], &c__1);
	i__1 = *m - 1;
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &x[1], &c__1, &
		c_b11, &dwork[iw21 + 1], &c__1, (ftnlen)1);

/*           Compute W. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ == j) {
		    i__3 = iz13 + i__ + (i__ - 1) * *n;
		    d__1 = (dwork[iw14 + 1] - dwork[iw13 + 1]) * 1e-4 / 2. - 
			    dlambd * dwork[iw21 + i__];
		    z__2.r = d__1, z__2.i = 0.;
		    i__4 = iz9 + i__ + (i__ - 1) * *n;
		    z__1.r = z__2.r + zwork[i__4].r, z__1.i = z__2.i + zwork[
			    i__4].i;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
		} else {
		    i__3 = iz13 + i__ + (j - 1) * *n;
		    i__4 = iz9 + i__ + (j - 1) * *n;
		    zwork[i__3].r = zwork[i__4].r, zwork[i__3].i = zwork[i__4]
			    .i;
		}
/* L490: */
	    }
/* L500: */
	}

/*           Compute eig( W ). */

	i__1 = *lzwork - izwrk;
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz13 + 1], n, &sdim, &zwork[
		iz14 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
	if (info2 > 0) {
	    *info = 6;
	    return 0;
	}
	i__1 = izwrk + 1;
	lza = (integer) zwork[i__1].r;
	lzamax = max(lza,lzamax);
	i__1 = iz14 + 1;
	emax = zwork[i__1].r;
	if (*n > 1) {
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		i__2 = iz14 + i__;
		if (zwork[i__2].r > emax) {
		    i__3 = iz14 + i__;
		    emax = zwork[i__3].r;
		}
/* L510: */
	    }
	}
	if (emax <= 0.) {
	    goto L515;
	} else {
	    ++l;
	    goto L460;
	}
    }

/*        Set y. */

L515:
    dwork[iw13 + 1] = dlambd;
    dcopy_(&mt, &x[1], &c__1, &dwork[iw13 + 2], &c__1);

    if (svlam - dlambd < tol) {
	*bound = sqrt((max(e,0.))) * znorm;
	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = dwork[iw2 + i__ + 1];
	    x[i__] *= d__1 * d__1;
/* L520: */
	}

/*           Compute sqrt( x ). */

	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw20 + i__] = sqrt(x[i__]);
/* L530: */
	}

/*           Compute diag( D ). */

	dcopy_(n, &dwork[iw7 + 1], &c__1, &d__[1], &c__1);
	i__1 = *m - 1;
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw20 + 1], &
		c__1, &c_b11, &d__[1], &c__1, (ftnlen)1);

/*           Compute diag( G ). */

	j = 0;
	l = 0;
	i__1 = *m;
	for (k = 1; k <= i__1; ++k) {
	    j += nblock[k];
	    if (itype[k] == 1) {
		++l;
/* Computing 2nd power */
		d__1 = dwork[iw2 + k];
		x[*m - 1 + l] *= d__1 * d__1;
		g[j] = x[*m - 1 + l];
	    }
/* L540: */
	}
	dscal_(n, &znorm, &g[1], &c__1);
	dwork[1] = (doublereal) (minwrk - *n * 5 + lwamax);
	i__1 = minzrk - *n * 3 + lzamax;
	z__1.r = (doublereal) i__1, z__1.i = 0.;
	zwork[1].r = z__1.r, zwork[1].i = z__1.i;
	return 0;
    }
    svlam = dlambd;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {

/*           Store xD. */

	i__2 = *m - 1;
	dcopy_(&i__2, &x[1], &c__1, &dwork[iw22 + 1], &c__1);
	if (mr > 0) {

/*              Store xG. */

	    dcopy_(&mr, &x[*m], &c__1, &dwork[iw23 + 1], &c__1);
	}

/*           Compute A(:) = A0 + AA*x. */

	i__2 = mt;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz10 + i__;
	    i__4 = i__;
	    z__1.r = x[i__4], z__1.i = 0.;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L550: */
	}
	i__2 = *n * *n;
	zcopy_(&i__2, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
	i__2 = *n * *n;
	i__3 = *n * *n;
	zgemv_("N", &i__2, &mt, &c_b2, &zwork[iz6 + 1], &i__3, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute B = B0d + BBd*xD. */

	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
	i__2 = *m - 1;
	dgemv_("N", n, &i__2, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &
		c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute F. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (i__ == j) {
		    i__4 = iz15 + i__ + (i__ - 1) * *n;
		    d__1 = dlambd * dwork[iw24 + i__];
		    z__2.r = d__1, z__2.i = 0.;
		    i__5 = iz7 + i__ + (i__ - 1) * *n;
		    z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - zwork[
			    i__5].i;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		} else {
		    i__4 = iz15 + i__ + (j - 1) * *n;
		    i__5 = iz7 + i__ + (j - 1) * *n;
		    z__1.r = -zwork[i__5].r, z__1.i = -zwork[i__5].i;
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		}
/* L555: */
	    }
/* L556: */
	}
	zlacpy_("Full", n, n, &zwork[iz15 + 1], n, &zwork[iz17 + 1], n, (
		ftnlen)4);

/*           Compute det( F ). */

	i__2 = *lzwork - izwrk;
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
		iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__2, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
	if (info2 > 0) {
	    *info = 6;
	    return 0;
	}
	i__2 = izwrk + 1;
	lza = (integer) zwork[i__2].r;
	lzamax = max(lza,lzamax);
	detf.r = 1., detf.i = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz16 + i__;
	    z__1.r = detf.r * zwork[i__3].r - detf.i * zwork[i__3].i, z__1.i =
		     detf.r * zwork[i__3].i + detf.i * zwork[i__3].r;
	    detf.r = z__1.r, detf.i = z__1.i;
/* L560: */
	}

/*           Compute Finv. */

	zgetrf_(n, n, &zwork[iz17 + 1], n, &iwork[1], &info2);
	if (info2 > 0) {
	    *info = 5;
	    return 0;
	}
	i__2 = *ldwork - iwrk;
	zgetri_(n, &zwork[iz17 + 1], n, &iwork[1], &zwork[izwrk + 1], &i__2, &
		info2);
	i__2 = izwrk + 1;
	lza = (integer) zwork[i__2].r;
	lzamax = max(lza,lzamax);

/*           Compute phi. */

	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
	    dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L570: */
	}
	if (mr > 0) {
	    i__2 = mr;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
		dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
			i__];
/* L580: */
	    }
	}
	prod = 1.;
	i__2 = mt << 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    prod *= dwork[iw25 + i__];
/* L590: */
	}
	temp = detf.r;
	if (temp < eps) {
	    temp = eps;
	}
	phi = -log(temp) - log(prod);

/*           Compute g. */

	i__2 = mt;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n * *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = iz18 + i__ + (j - 1) * *n * *n;
		d__1 = dlambd * dwork[iw9 + i__ + (j - 1) * *n * *n];
		z__2.r = d__1, z__2.i = 0.;
		i__5 = iz6 + i__ + (j - 1) * *n * *n;
		z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - zwork[i__5]
			.i;
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L600: */
	    }
/* L610: */
	}
	i__2 = *n * *n;
	i__3 = *n * *n;
	zgemv_("C", &i__2, &mt, &c_b2, &zwork[iz18 + 1], &i__3, &zwork[iz17 + 
		1], &c__1, &c_b1, &zwork[iz19 + 1], &c__1, (ftnlen)1);
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[iw26 + i__] = 1. / (dwork[iw22 + i__] - .01) - 1. / (100. - 
		    dwork[iw22 + i__]);
/* L620: */
	}
	if (mr > 0) {
	    i__2 = mr;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[iw26 + *m - 1 + i__] = 1. / (dwork[iw23 + i__] + tau) - 
			1. / (tau - dwork[iw23 + i__]);
/* L630: */
	    }
	}
	i__2 = mt;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz19 + i__;
	    dwork[iw26 + i__] = -zwork[i__3].r - dwork[iw26 + i__];
/* L640: */
	}

/*           Compute h. */

	dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &
		mt, (ftnlen)4);
	dcopy_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
	i__2 = *ldwork - iwrk;
	dsysv_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw27 
		+ 1], &mt, &dwork[iwrk + 1], &i__2, &info2, (ftnlen)1);
	if (info2 > 0) {
	    *info = 5;
	    return 0;
	}
	lwa = (integer) dwork[iwrk + 1];
	lwamax = max(lwa,lwamax);
	stsize = 1.;

/*           Store hD. */

	i__2 = *m - 1;
	dcopy_(&i__2, &dwork[iw27 + 1], &c__1, &dwork[iw28 + 1], &c__1);

/*           Determine stepsize. */

	l = 0;
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (dwork[iw28 + i__] > 0.) {
		++l;
		if (l == 1) {
		    temp = (dwork[iw22 + i__] - .01) / dwork[iw28 + i__];
		} else {
/* Computing MIN */
		    d__1 = temp, d__2 = (dwork[iw22 + i__] - .01) / dwork[
			    iw28 + i__];
		    temp = min(d__1,d__2);
		}
	    }
/* L650: */
	}
	if (l > 0) {
	    stsize = min(stsize,temp);
	}
	l = 0;
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (dwork[iw28 + i__] < 0.) {
		++l;
		if (l == 1) {
		    temp = (100. - dwork[iw22 + i__]) / (-dwork[iw28 + i__]);
		} else {
/* Computing MIN */
		    d__1 = temp, d__2 = (100. - dwork[iw22 + i__]) / (-dwork[
			    iw28 + i__]);
		    temp = min(d__1,d__2);
		}
	    }
/* L660: */
	}
	if (l > 0) {
	    stsize = min(stsize,temp);
	}
	if (mr > 0) {

/*              Store hG. */

	    dcopy_(&mr, &dwork[iw27 + *m], &c__1, &dwork[iw29 + 1], &c__1);

/*              Determine stepsize. */

	    l = 0;
	    i__2 = mr;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (dwork[iw29 + i__] > 0.) {
		    ++l;
		    if (l == 1) {
			temp = (dwork[iw23 + i__] + tau) / dwork[iw29 + i__];
		    } else {
/* Computing MIN */
			d__1 = temp, d__2 = (dwork[iw23 + i__] + tau) / dwork[
				iw29 + i__];
			temp = min(d__1,d__2);
		    }
		}
/* L670: */
	    }
	    if (l > 0) {
		stsize = min(stsize,temp);
	    }
	    l = 0;
	    i__2 = mr;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (dwork[iw29 + i__] < 0.) {
		    ++l;
		    if (l == 1) {
			temp = (tau - dwork[iw23 + i__]) / (-dwork[iw29 + i__]
				);
		    } else {
/* Computing MIN */
			d__1 = temp, d__2 = (tau - dwork[iw23 + i__]) / (
				-dwork[iw29 + i__]);
			temp = min(d__1,d__2);
		    }
		}
/* L680: */
	    }
	}
	if (l > 0) {
	    stsize = min(stsize,temp);
	}
	stsize *= .9;
	if (stsize >= tol4) {

/*              Compute x_new. */

	    i__2 = mt;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[iw20 + i__] = x[i__] - stsize * dwork[iw27 + i__];
/* L700: */
	    }

/*              Store xD. */

	    i__2 = *m - 1;
	    dcopy_(&i__2, &dwork[iw20 + 1], &c__1, &dwork[iw22 + 1], &c__1);
	    if (mr > 0) {

/*                 Store xG. */

		dcopy_(&mr, &dwork[iw20 + *m], &c__1, &dwork[iw23 + 1], &c__1)
			;
	    }

/*              Compute A(:) = A0 + AA*x_new. */

	    i__2 = mt;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = iz10 + i__;
		i__4 = iw20 + i__;
		z__1.r = dwork[i__4], z__1.i = 0.;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L710: */
	    }
	    i__2 = *n * *n;
	    zcopy_(&i__2, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
	    i__2 = *n * *n;
	    i__3 = *n * *n;
	    zgemv_("N", &i__2, &mt, &c_b2, &zwork[iz6 + 1], &i__3, &zwork[
		    iz10 + 1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)
		    1);

/*              Compute B = B0d + BBd*xD. */

	    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
	    i__2 = *m - 1;
	    dgemv_("N", n, &i__2, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1]
		    , &c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*              Compute lambda*diag(B) - A. */

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (i__ == j) {
			i__4 = iz15 + i__ + (i__ - 1) * *n;
			d__1 = dlambd * dwork[iw24 + i__];
			z__2.r = d__1, z__2.i = 0.;
			i__5 = iz7 + i__ + (i__ - 1) * *n;
			z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - 
				zwork[i__5].i;
			zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		    } else {
			i__4 = iz15 + i__ + (j - 1) * *n;
			i__5 = iz7 + i__ + (j - 1) * *n;
			z__1.r = -zwork[i__5].r, z__1.i = -zwork[i__5].i;
			zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		    }
/* L720: */
		}
/* L730: */
	    }

/*              Compute eig( lambda*diag(B)-A ). */

	    i__2 = *lzwork - izwrk;
	    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &
		    zwork[iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__2, &
		    dwork[iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
	    if (info2 > 0) {
		*info = 6;
		return 0;
	    }
	    i__2 = izwrk + 1;
	    lza = (integer) zwork[i__2].r;
	    lzamax = max(lza,lzamax);
	    i__2 = iz16 + 1;
	    emin = zwork[i__2].r;
	    if (*n > 1) {
		i__2 = *n;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__3 = iz16 + i__;
		    if (zwork[i__3].r < emin) {
			i__4 = iz16 + i__;
			emin = zwork[i__4].r;
		    }
/* L740: */
		}
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = iz16 + i__;
		dwork[iw30 + i__] = zwork[i__3].r;
/* L750: */
	    }
	    i__2 = *m - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[iw30 + *n + i__] = dwork[iw22 + i__] - .01;
		dwork[iw30 + *n + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L760: */
	    }
	    if (mr > 0) {
		i__2 = mr;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[iw30 + *n + (*m - 1 << 1) + i__] = dwork[iw23 + i__]
			     + tau;
		    dwork[iw30 + *n + (*m - 1 << 1) + mr + i__] = tau - dwork[
			    iw23 + i__];
/* L770: */
		}
	    }
	    prod = 1.;
	    i__2 = *n + (mt << 1);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		prod *= dwork[iw30 + i__];
/* L780: */
	    }
	    if (emin <= 0. || -log(prod) >= phi) {
		stsize /= 10.;
	    } else {
		dcopy_(&mt, &dwork[iw20 + 1], &c__1, &x[1], &c__1);
	    }
	}
	if (stsize < tol4) {
	    goto L810;
	}
/* L800: */
    }

L810:

/*           Store xD. */

    i__1 = *m - 1;
    dcopy_(&i__1, &x[1], &c__1, &dwork[iw22 + 1], &c__1);
    if (mr > 0) {

/*              Store xG. */

	dcopy_(&mr, &x[*m], &c__1, &dwork[iw23 + 1], &c__1);
    }

/*           Compute A(:) = A0 + AA*x. */

    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz10 + i__;
	i__3 = i__;
	z__1.r = x[i__3], z__1.i = 0.;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L820: */
    }
    i__1 = *n * *n;
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
    i__1 = *m - 1;
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute F. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		i__3 = iz15 + i__ + (i__ - 1) * *n;
		d__1 = dlambd * dwork[iw24 + i__];
		z__2.r = d__1, z__2.i = 0.;
		i__4 = iz7 + i__ + (i__ - 1) * *n;
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    } else {
		i__3 = iz15 + i__ + (j - 1) * *n;
		i__4 = iz7 + i__ + (j - 1) * *n;
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    }
/* L830: */
	}
/* L840: */
    }
    zlacpy_("Full", n, n, &zwork[iz15 + 1], n, &zwork[iz17 + 1], n, (ftnlen)4)
	    ;

/*           Compute det( F ). */

    i__1 = *lzwork - izwrk;
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    detf.r = 1., detf.i = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz16 + i__;
	z__1.r = detf.r * zwork[i__2].r - detf.i * zwork[i__2].i, z__1.i = 
		detf.r * zwork[i__2].i + detf.i * zwork[i__2].r;
	detf.r = z__1.r, detf.i = z__1.i;
/* L850: */
    }

/*           Compute Finv. */

    zgetrf_(n, n, &zwork[iz17 + 1], n, &iwork[1], &info2);
    if (info2 > 0) {
	*info = 5;
	return 0;
    }
    i__1 = *ldwork - iwrk;
    zgetri_(n, &zwork[iz17 + 1], n, &iwork[1], &zwork[izwrk + 1], &i__1, &
	    info2);
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);

/*           Compute the barrier function. */

    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L860: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
/* L870: */
	}
    }
    prod = 1.;
    i__1 = mt << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prod *= dwork[iw25 + i__];
/* L880: */
    }
    temp = detf.r;
    if (temp < eps) {
	temp = eps;
    }
    phi = -log(temp) - log(prod);

/*           Compute the gradient of the barrier function. */

    i__1 = mt;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n * *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz18 + i__ + (j - 1) * *n * *n;
	    d__1 = dlambd * dwork[iw9 + i__ + (j - 1) * *n * *n];
	    z__2.r = d__1, z__2.i = 0.;
	    i__4 = iz6 + i__ + (j - 1) * *n * *n;
	    z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4].i;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L890: */
	}
/* L900: */
    }
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("C", &i__1, &mt, &c_b2, &zwork[iz18 + 1], &i__2, &zwork[iz17 + 1], 
	    &c__1, &c_b1, &zwork[iz19 + 1], &c__1, (ftnlen)1);
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw26 + i__] = 1. / (dwork[iw22 + i__] - .01) - 1. / (100. - 
		dwork[iw22 + i__]);
/* L910: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw26 + *m - 1 + i__] = 1. / (dwork[iw23 + i__] + tau) - 1. /
		     (tau - dwork[iw23 + i__]);
/* L920: */
	}
    }
    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz19 + i__;
	dwork[iw26 + i__] = -zwork[i__2].r - dwork[iw26 + i__];
/* L925: */
    }

/*           Compute the Hessian of the barrier function. */

    i__1 = *n * mt;
    zgemm_("N", "N", n, &i__1, n, &c_b2, &zwork[iz17 + 1], n, &zwork[iz18 + 1]
	    , n, &c_b1, &zwork[iz20 + 1], n, (ftnlen)1, (ftnlen)1);
    dlaset_("Full", &mt, &mt, &c_b15, &c_b15, &dwork[iw11 + 1], &mt, (ftnlen)
	    4);
    i__1 = mt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n * *n;
	zcopy_(&i__2, &zwork[iz20 + 1 + (k - 1) * *n * *n], &c__1, &zwork[
		iz22 + 1], &c__1);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = iz23 + i__ + (j - 1) * *n;
		d_cnjg(&z__1, &zwork[iz22 + j + (i__ - 1) * *n]);
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L930: */
	    }
/* L940: */
	}
	i__2 = *n * *n;
	i__3 = *n * *n;
	zgemv_("C", &i__2, &k, &c_b2, &zwork[iz20 + 1], &i__3, &zwork[iz23 + 
		1], &c__1, &c_b1, &zwork[iz24 + 1], &c__1, (ftnlen)1);
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    d_cnjg(&z__1, &zwork[iz24 + j]);
	    dwork[iw11 + k + (j - 1) * mt] = z__1.r;
/* L950: */
	}
/* L960: */
    }
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = dwork[iw22 + i__] - .01;
/* Computing 2nd power */
	d__2 = 100. - dwork[iw22 + i__];
	dwork[iw10 + i__] = 1. / (d__1 * d__1) + 1. / (d__2 * d__2);
/* L970: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = dwork[iw23 + i__] + tau;
/* Computing 2nd power */
	    d__2 = tau - dwork[iw23 + i__];
	    dwork[iw10 + *m - 1 + i__] = 1. / (d__1 * d__1) + 1. / (d__2 * 
		    d__2);
/* L980: */
	}
    }
    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw11 + i__ + (i__ - 1) * mt] += dwork[iw10 + i__];
/* L990: */
    }
    i__1 = mt;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ != j) {
		t1 = dwork[iw11 + i__ + (j - 1) * mt];
		t2 = dwork[iw11 + j + (i__ - 1) * mt];
		dwork[iw11 + i__ + (j - 1) * mt] = t1 + t2;
		dwork[iw11 + j + (i__ - 1) * mt] = t1 + t2;
	    }
/* L1000: */
	}
/* L1100: */
    }

/*           Compute norm( H ). */

L1110:
    hnorm = dlange_("F", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[1], (ftnlen)
	    1);

/*           Compute rcond( H ). */

    dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &mt, (
	    ftnlen)4);
    hnorm1 = dlange_("1", &mt, &mt, &dwork[iw31 + 1], &mt, &dwork[1], (ftnlen)
	    1);
    i__1 = *ldwork - iwrk;
    dsytrf_("U", &mt, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iwrk + 1], &
	    i__1, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 5;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1];
    lwamax = max(lwa,lwamax);
    dsycon_("U", &mt, &dwork[iw31 + 1], &mt, &iwork[1], &hnorm1, &rcond, &
	    dwork[iwrk + 1], &iwork[mt + 1], &info2, (ftnlen)1);
    if (rcond < tol3) {
	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw11 + i__ + (i__ - 1) * mt] += hnorm * regpar;
/* L1120: */
	}
	goto L1110;
    }

/*           Compute the tangent line to path of center. */

    dcopy_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
    dsytrs_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw27 + 
	    1], &mt, &info2, (ftnlen)1);

/*           Check if x-h satisfies the Goldstein test. */

    gtest = FALSE_;
    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw20 + i__] = x[i__] - dwork[iw27 + i__];
/* L1130: */
    }

/*           Store xD. */

    i__1 = *m - 1;
    dcopy_(&i__1, &dwork[iw20 + 1], &c__1, &dwork[iw22 + 1], &c__1);
    if (mr > 0) {

/*              Store xG. */

	dcopy_(&mr, &dwork[iw20 + *m], &c__1, &dwork[iw23 + 1], &c__1);
    }

/*           Compute A(:) = A0 + AA*x_new. */

    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz10 + i__;
	i__3 = iw20 + i__;
	z__1.r = dwork[i__3], z__1.i = 0.;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L1140: */
    }
    i__1 = *n * *n;
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
    i__1 = *m - 1;
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute lambda*diag(B) - A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		i__3 = iz15 + i__ + (i__ - 1) * *n;
		d__1 = dlambd * dwork[iw24 + i__];
		z__2.r = d__1, z__2.i = 0.;
		i__4 = iz7 + i__ + (i__ - 1) * *n;
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    } else {
		i__3 = iz15 + i__ + (j - 1) * *n;
		i__4 = iz7 + i__ + (j - 1) * *n;
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    }
/* L1150: */
	}
/* L1160: */
    }

/*           Compute eig( lambda*diag(B)-A ). */

    i__1 = *lzwork - izwrk;
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz16 + i__;
	dwork[iw30 + i__] = zwork[i__2].r;
/* L1190: */
    }
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw30 + *n + i__] = dwork[iw22 + i__] - .01;
	dwork[iw30 + *n + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L1200: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw30 + *n + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
	    dwork[iw30 + *n + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
		    i__];
/* L1210: */
	}
    }
    emin = dwork[iw30 + 1];
    i__1 = *n + (mt << 1);
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dwork[iw30 + i__] < emin) {
	    emin = dwork[iw30 + i__];
	}
/* L1220: */
    }
    if (emin <= 0.) {
	gtest = FALSE_;
    } else {
	pp = ddot_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
	prod = 1.;
	i__1 = *n + (mt << 1);
	for (i__ = 1; i__ <= i__1; ++i__) {
	    prod *= dwork[iw30 + i__];
/* L1230: */
	}
	t1 = -log(prod);
	t2 = phi - pp * .01;
	t3 = phi - pp * .9;
	if (t1 >= t3 && t1 < t2) {
	    gtest = TRUE_;
	}
    }

/*           Use x-h if Goldstein test is satisfied. Otherwise use */
/*           Nesterov-Nemirovsky's stepsize length. */

    pp = ddot_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
    delta = sqrt(pp);
    if (gtest || delta <= .25) {
	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] -= dwork[iw27 + i__];
/* L1240: */
	}
    } else {
	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] -= dwork[iw27 + i__] / (delta + 1.);
/* L1250: */
	}
    }

/*           Analytic center is found if delta is sufficiently small. */

    if (delta < tol5) {
	goto L1260;
    }
    goto L810;

/*        Set yf. */

L1260:
    dwork[iw14 + 1] = dlambd;
    dcopy_(&mt, &x[1], &c__1, &dwork[iw14 + 2], &c__1);

/*        Set yw. */

    i__1 = mt + 1;
    dcopy_(&i__1, &dwork[iw14 + 1], &c__1, &dwork[iw15 + 1], &c__1);

/*        Compute Fb. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = iz21 + i__ + (j - 1) * *n;
	    i__4 = iw24 + i__;
	    z__2.r = dwork[i__4], z__2.i = 0.;
	    d_cnjg(&z__3, &zwork[iz17 + j + (i__ - 1) * *n]);
	    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
		    z__3.i + z__2.i * z__3.r;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L1270: */
	}
/* L1280: */
    }
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("C", &i__1, &mt, &c_b2, &zwork[iz20 + 1], &i__2, &zwork[iz21 + 1], 
	    &c__1, &c_b1, &zwork[iz24 + 1], &c__1, (ftnlen)1);
    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz24 + i__;
	dwork[iw32 + i__] = zwork[i__2].r;
/* L1300: */
    }

/*        Compute h1. */

    dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &mt, (
	    ftnlen)4);
    i__1 = *ldwork - iwrk;
    dsysv_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw32 + 1]
	    , &mt, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1);
    if (info2 > 0) {
	*info = 5;
	return 0;
    }
    lwa = (integer) dwork[iwrk + 1];
    lwamax = max(lwa,lwamax);

/*        Compute hn. */

    hn = dlange_("F", &mt, &c__1, &dwork[iw32 + 1], &mt, &dwork[1], (ftnlen)1)
	    ;

/*        Compute y. */

    dwork[iw13 + 1] = dlambd - c__ / hn;
    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw13 + 1 + i__] = x[i__] + c__ * dwork[iw32 + i__] / hn;
/* L1310: */
    }

/*        Store xD. */

    i__1 = *m - 1;
    dcopy_(&i__1, &dwork[iw13 + 2], &c__1, &dwork[iw22 + 1], &c__1);
    if (mr > 0) {

/*           Store xG. */

	dcopy_(&mr, &dwork[iw13 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1);
    }

/*        Compute A(:) = A0 + AA*y(2:mt+1). */

    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz10 + i__;
	i__3 = iw13 + 1 + i__;
	z__1.r = dwork[i__3], z__1.i = 0.;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L1320: */
    }
    i__1 = *n * *n;
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute B = B0d + BBd*xD. */

    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
    i__1 = *m - 1;
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*        Compute y(1)*diag(B) - A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		i__3 = iz15 + i__ + (i__ - 1) * *n;
		d__1 = dwork[iw13 + 1] * dwork[iw24 + i__];
		z__2.r = d__1, z__2.i = 0.;
		i__4 = iz7 + i__ + (i__ - 1) * *n;
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    } else {
		i__3 = iz15 + i__ + (j - 1) * *n;
		i__4 = iz7 + i__ + (j - 1) * *n;
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    }
/* L1330: */
	}
/* L1340: */
    }

/*        Compute eig( y(1)*diag(B)-A ). */

    i__1 = *lzwork - izwrk;
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    i__1 = iz16 + 1;
    emin = zwork[i__1].r;
    if (*n > 1) {
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = iz16 + i__;
	    if (zwork[i__2].r < emin) {
		i__3 = iz16 + i__;
		emin = zwork[i__3].r;
	    }
/* L1350: */
	}
    }
    pos = TRUE_;
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L1360: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
/* L1370: */
	}
    }
    temp = dwork[iw25 + 1];
    i__1 = mt << 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dwork[iw25 + i__] < temp) {
	    temp = dwork[iw25 + i__];
	}
/* L1380: */
    }
    if (temp <= 0. || emin <= 0.) {
	pos = FALSE_;
    }
L1390:
    if (pos) {

/*           Set y2 = y. */

	i__1 = mt + 1;
	dcopy_(&i__1, &dwork[iw13 + 1], &c__1, &dwork[iw17 + 1], &c__1);

/*           Compute y = y + 1.5*( y - yw ). */

	i__1 = mt + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw13 + i__] += (dwork[iw13 + i__] - dwork[iw15 + i__]) * 
		    1.5;
/* L1400: */
	}

/*           Store xD. */

	i__1 = *m - 1;
	dcopy_(&i__1, &dwork[iw13 + 2], &c__1, &dwork[iw22 + 1], &c__1);
	if (mr > 0) {

/*              Store xG. */

	    dcopy_(&mr, &dwork[iw13 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1)
		    ;
	}

/*           Compute A(:) = A0 + AA*y(2:mt+1). */

	i__1 = mt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = iz10 + i__;
	    i__3 = iw13 + 1 + i__;
	    z__1.r = dwork[i__3], z__1.i = 0.;
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L1420: */
	}
	i__1 = *n * *n;
	zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
	i__1 = *n * *n;
	i__2 = *n * *n;
	zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
	i__1 = *m - 1;
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &
		c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Set yw = y2. */

	i__1 = mt + 1;
	dcopy_(&i__1, &dwork[iw17 + 1], &c__1, &dwork[iw15 + 1], &c__1);

/*           Compute y(1)*diag(B) - A. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ == j) {
		    i__3 = iz15 + i__ + (i__ - 1) * *n;
		    d__1 = dwork[iw13 + 1] * dwork[iw24 + i__];
		    z__2.r = d__1, z__2.i = 0.;
		    i__4 = iz7 + i__ + (i__ - 1) * *n;
		    z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[
			    i__4].i;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
		} else {
		    i__3 = iz15 + i__ + (j - 1) * *n;
		    i__4 = iz7 + i__ + (j - 1) * *n;
		    z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
		}
/* L1430: */
	    }
/* L1440: */
	}

/*           Compute eig( y(1)*diag(B)-A ). */

	i__1 = *lzwork - izwrk;
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
		iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
	if (info2 > 0) {
	    *info = 6;
	    return 0;
	}
	i__1 = izwrk + 1;
	lza = (integer) zwork[i__1].r;
	lzamax = max(lza,lzamax);
	i__1 = iz16 + 1;
	emin = zwork[i__1].r;
	if (*n > 1) {
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		i__2 = iz16 + i__;
		if (zwork[i__2].r < emin) {
		    i__3 = iz16 + i__;
		    emin = zwork[i__3].r;
		}
/* L1450: */
	    }
	}
	pos = TRUE_;
	i__1 = *m - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
	    dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L1460: */
	}
	if (mr > 0) {
	    i__1 = mr;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
		dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
			i__];
/* L1470: */
	    }
	}
	temp = dwork[iw25 + 1];
	i__1 = mt << 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if (dwork[iw25 + i__] < temp) {
		temp = dwork[iw25 + i__];
	    }
/* L1480: */
	}
	if (temp <= 0. || emin <= 0.) {
	    pos = FALSE_;
	}
	goto L1390;
    }
L1490:

/*        Set y1 = ( y + yw ) / 2. */

    i__1 = mt + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw16 + i__] = (dwork[iw13 + i__] + dwork[iw15 + i__]) / 2.;
/* L1500: */
    }

/*        Store xD. */

    i__1 = *m - 1;
    dcopy_(&i__1, &dwork[iw16 + 2], &c__1, &dwork[iw22 + 1], &c__1);
    if (mr > 0) {

/*           Store xG. */

	dcopy_(&mr, &dwork[iw16 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1);
    }

/*        Compute A(:) = A0 + AA*y1(2:mt+1). */

    i__1 = mt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iz10 + i__;
	i__3 = iw16 + 1 + i__;
	z__1.r = dwork[i__3], z__1.i = 0.;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L1510: */
    }
    i__1 = *n * *n;
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
    i__1 = *n * *n;
    i__2 = *n * *n;
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute diag( B ) = B0d + BBd*xD. */

    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
    i__1 = *m - 1;
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*        Compute y1(1)*diag(B) - A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		i__3 = iz15 + i__ + (i__ - 1) * *n;
		d__1 = dwork[iw16 + 1] * dwork[iw24 + i__];
		z__2.r = d__1, z__2.i = 0.;
		i__4 = iz7 + i__ + (i__ - 1) * *n;
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    } else {
		i__3 = iz15 + i__ + (j - 1) * *n;
		i__4 = iz7 + i__ + (j - 1) * *n;
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    }
/* L1520: */
	}
/* L1530: */
    }

/*        Compute eig( y1(1)*diag(B)-A ). */

    i__1 = *lzwork - izwrk;
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
    if (info2 > 0) {
	*info = 6;
	return 0;
    }
    i__1 = izwrk + 1;
    lza = (integer) zwork[i__1].r;
    lzamax = max(lza,lzamax);
    i__1 = iz16 + 1;
    emin = zwork[i__1].r;
    if (*n > 1) {
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = iz16 + i__;
	    if (zwork[i__2].r < emin) {
		i__3 = iz16 + i__;
		emin = zwork[i__3].r;
	    }
/* L1540: */
	}
    }
    pos = TRUE_;
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
/* L1550: */
    }
    if (mr > 0) {
	i__1 = mr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
/* L1560: */
	}
    }
    temp = dwork[iw25 + 1];
    i__1 = mt << 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dwork[iw25 + i__] < temp) {
	    temp = dwork[iw25 + i__];
	}
/* L1570: */
    }
    if (temp <= 0. || emin <= 0.) {
	pos = FALSE_;
    }
    if (pos) {

/*           Set yw = y1. */

	i__1 = mt + 1;
	dcopy_(&i__1, &dwork[iw16 + 1], &c__1, &dwork[iw15 + 1], &c__1);
    } else {

/*           Set y = y1. */

	i__1 = mt + 1;
	dcopy_(&i__1, &dwork[iw16 + 1], &c__1, &dwork[iw13 + 1], &c__1);
    }
    i__1 = mt + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw33 + i__] = dwork[iw13 + i__] - dwork[iw15 + i__];
/* L1580: */
    }
    i__1 = mt + 1;
    i__2 = mt + 1;
    ynorm1 = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);
    i__1 = mt + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw33 + i__] = dwork[iw13 + i__] - dwork[iw14 + i__];
/* L1590: */
    }
    i__1 = mt + 1;
    i__2 = mt + 1;
    ynorm2 = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);
    if (ynorm1 < ynorm2 * .01) {
	goto L1600;
    }
    goto L1490;

/*        Compute c. */

L1600:
    i__1 = mt + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iw33 + i__] = dwork[iw15 + i__] - dwork[iw14 + i__];
/* L1610: */
    }
    i__1 = mt + 1;
    i__2 = mt + 1;
    c__ = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);

/*        Set x = yw(2:mt+1). */

    dcopy_(&mt, &dwork[iw15 + 2], &c__1, &x[1], &c__1);
    goto L390;

/* *** Last line of AB13MD *** */
} /* ab13md_ */

