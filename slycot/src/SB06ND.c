/* SB06ND.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 0.;
static doublereal c_b22 = -1.;
static doublereal c_b25 = 1.;

/* Subroutine */ int sb06nd_(integer *n, integer *m, integer *kmax, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	kstair, doublereal *u, integer *ldu, doublereal *f, integer *ldf, 
	doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, u_dim1, 
	    u_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer j, j0, kk, kmin, jcur, kcur;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer jkcur;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer mkcur, ncont, kstep;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarfg_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer jmkcur;
    extern /* Subroutine */ int dlatzm_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, ftnlen);


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

/*     To construct the minimum norm feedback matrix F to perform */
/*     "deadbeat control" on a (A,B)-pair of a state-space model (which */
/*     must be preliminarily reduced to upper "staircase" form using */
/*     SLICOT Library routine AB01OD) such that the matrix R = A + BFU' */
/*     is nilpotent. */
/*     (The transformation matrix U reduces R to upper Schur form with */
/*     zero blocks on its diagonal (of dimension KSTAIR(i)) and */
/*     therefore contains bases for the i-th controllable subspaces, */
/*     where i = 1,...,KMAX). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the */
/*             matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The actual input dimension.  M >= 0. */

/*     KMAX    (input) INTEGER */
/*             The number of "stairs" in the staircase form as produced */
/*             by SLICOT Library routine AB01OD.  0 <= KMAX <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the transformed state-space matrix of the */
/*             (A,B)-pair with triangular stairs, as produced by SLICOT */
/*             Library routine AB01OD (with option STAGES = 'A'). */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix U'AU + U'BF. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the transformed triangular input matrix of the */
/*             (A,B)-pair as produced by SLICOT Library routine AB01OD */
/*             (with option STAGES = 'A'). */
/*             On exit, the leading N-by-M part of this array contains */
/*             the matrix U'B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     KSTAIR  (input) INTEGER array, dimension (KMAX) */
/*             The leading KMAX elements of this array must contain the */
/*             dimensions of each "stair" as produced by SLICOT Library */
/*             routine AB01OD. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain either a transformation matrix (e.g. from a */
/*             previous call to other SLICOT routine) or be initialised */
/*             as the identity matrix. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the product of the input matrix U and the state-space */
/*             transformation matrix which reduces A + BFU' to real */
/*             Schur form. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,N). */

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the */
/*             deadbeat feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Starting from the (A,B)-pair in "staircase form" with "triangular" */
/*     stairs, dimensions KSTAIR(i+1) x KSTAIR(i), (described by the */
/*     vector KSTAIR): */

/*                    | B | A      *  . . .  *  | */
/*                    |  1|  11       .      .  | */
/*                    |   | A     A     .    .  | */
/*                    |   |  21    22     .  .  | */
/*                    |   |    .      .     .   | */
/*      [ B | A ]  =  |   |      .      .    *  | */
/*                    |   |        .      .     | */
/*                    | 0 |   0                 | */
/*                    |   |          A      A   | */
/*                    |   |           r,r-1  rr | */

/*     where the i-th diagonal block of A has dimension KSTAIR(i), for */
/*     i = 1,2,...,r, the feedback matrix F is constructed recursively in */
/*     r steps (where the number of "stairs" r is given by KMAX). In each */
/*     step a unitary state-space transformation U and a part of F are */
/*     updated in order to achieve the final form: */

/*                       | 0   A      *   . . .  *   | */
/*                       |      12      .        .   | */
/*                       |                .      .   | */
/*                       |     0    A       .    .   | */
/*                       |           23       .  .   | */
/*                       |         .      .          | */
/*     [ U'AU + U'BF ] = |           .      .    *   | . */
/*                       |             .      .      | */
/*                       |                           | */
/*                       |                     A     | */
/*                       |                      r-1,r| */
/*                       |                           | */
/*                       |                       0   | */


/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         Deadbeat control: a special inverse eigenvalue problem. */
/*         BIT, 24, pp. 681-699, 1984. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O((N + M) * N**2) operations and is mixed */
/*     numerical stable (see [1]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB06BD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     1997, December 10; 2003, September 27. */

/*     KEYWORDS */

/*     Canonical form, deadbeat control, eigenvalue assignment, feedback */
/*     control, orthogonal transformation, real Schur form, staircase */
/*     form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --kstair;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*kmax < 0 || *kmax > *n) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldu < max(1,*n)) {
	*info = -10;
    } else if (*ldf < max(1,*m)) {
	*info = -12;
    } else {
	ncont = 0;

	i__1 = *kmax;
	for (kk = 1; kk <= i__1; ++kk) {
	    ncont += kstair[kk];
/* L10: */
	}

	if (ncont > *n) {
	    *info = -8;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB06ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0) {
	return 0;
    }

    i__1 = *kmax;
    for (kmin = 1; kmin <= i__1; ++kmin) {
	jcur = ncont;
	kstep = *kmax - kmin;

/*        Triangularize bottom part of A (if KSTEP > 0). */

	i__2 = *kmax - kstep + 1;
	for (kk = *kmax; kk >= i__2; --kk) {
	    kcur = kstair[kk];

/*           Construct Ukk and store in Fkk. */

	    i__3 = kcur;
	    for (j = 1; j <= i__3; ++j) {
		jmkcur = jcur - kcur;
		dcopy_(&kcur, &a[jcur + jmkcur * a_dim1], lda, &f[jcur * 
			f_dim1 + 1], &c__1);
		i__4 = kcur + 1;
		dlarfg_(&i__4, &a[jcur + jcur * a_dim1], &f[jcur * f_dim1 + 1]
			, &c__1, &dwork[jcur]);
		dlaset_("Full", &c__1, &kcur, &c_b11, &c_b11, &a[jcur + 
			jmkcur * a_dim1], lda, (ftnlen)4);

/*              Backmultiply A and U with Ukk. */

		i__4 = jcur - 1;
		i__5 = kcur + 1;
		dlatzm_("Right", &i__4, &i__5, &f[jcur * f_dim1 + 1], &c__1, &
			dwork[jcur], &a[jcur * a_dim1 + 1], &a[jmkcur * 
			a_dim1 + 1], lda, &dwork[1], (ftnlen)5);

		i__4 = kcur + 1;
		dlatzm_("Right", n, &i__4, &f[jcur * f_dim1 + 1], &c__1, &
			dwork[jcur], &u[jcur * u_dim1 + 1], &u[jmkcur * 
			u_dim1 + 1], ldu, &dwork[*n + 1], (ftnlen)5);
		--jcur;
/* L20: */
	    }

/* L40: */
	}

/*        Eliminate diagonal block Aii by feedback Fi. */

	kcur = kstair[kmin];
	j0 = jcur - kcur + 1;
	mkcur = *m - kcur + 1;

/*        Solve for Fi and add B x Fi to A. */

	dlacpy_("Full", &kcur, &kcur, &a[j0 + j0 * a_dim1], lda, &f[mkcur + 
		j0 * f_dim1], ldf, (ftnlen)4);
	dtrsm_("Left", "Upper", "No transpose", "Non-unit", &kcur, &kcur, &
		c_b22, &b[j0 + mkcur * b_dim1], ldb, &f[mkcur + j0 * f_dim1], 
		ldf, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	if (j0 > 1) {
	    i__2 = j0 - 1;
	    dgemm_("No transpose", "No transpose", &i__2, &kcur, &kcur, &
		    c_b25, &b[mkcur * b_dim1 + 1], ldb, &f[mkcur + j0 * 
		    f_dim1], ldf, &c_b25, &a[j0 * a_dim1 + 1], lda, (ftnlen)
		    12, (ftnlen)12);
	}
	dlaset_("Full", &kcur, &kcur, &c_b11, &c_b11, &a[j0 + j0 * a_dim1], 
		lda, (ftnlen)4);
	i__2 = *m - kcur;
	dlaset_("Full", &i__2, &kcur, &c_b11, &c_b11, &f[j0 * f_dim1 + 1], 
		ldf, (ftnlen)4);

	if (kstep != 0) {
	    jkcur = ncont;

/*           Premultiply A with Ukk. */

	    i__2 = *kmax - kstep + 1;
	    for (kk = *kmax; kk >= i__2; --kk) {
		kcur = kstair[kk];
		jcur = jkcur - kcur;

		i__3 = kcur;
		for (j = 1; j <= i__3; ++j) {
		    i__4 = kcur + 1;
		    i__5 = *n - jcur + 1;
		    dlatzm_("Left", &i__4, &i__5, &f[jkcur * f_dim1 + 1], &
			    c__1, &dwork[jkcur], &a[jkcur + jcur * a_dim1], &
			    a[jcur + jcur * a_dim1], lda, &dwork[*n + 1], (
			    ftnlen)4);
		    --jcur;
		    --jkcur;
/* L60: */
		}

/* L80: */
	    }

/*           Premultiply B with Ukk. */

	    jcur += kcur;
	    jkcur = jcur + kcur;

	    i__2 = *m - kcur + 1;
	    for (j = *m; j >= i__2; --j) {
		i__3 = kcur + 1;
		i__4 = *m - j + 1;
		dlatzm_("Left", &i__3, &i__4, &f[jkcur * f_dim1 + 1], &c__1, &
			dwork[jkcur], &b[jkcur + j * b_dim1], &b[jcur + j * 
			b_dim1], ldb, &dwork[*n + 1], (ftnlen)4);
		--jcur;
		--jkcur;
/* L100: */
	    }

	}
/* L120: */
    }

    if (ncont != *n) {
	i__1 = *n - ncont;
	dlaset_("Full", m, &i__1, &c_b11, &c_b11, &f[(ncont + 1) * f_dim1 + 1]
		, ldf, (ftnlen)4);
    }

    return 0;
/* *** Last line of SB06ND *** */
} /* sb06nd_ */

