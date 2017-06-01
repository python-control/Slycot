/* MB04XY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04xy_(char *jobu, char *jobv, integer *m, integer *n, 
	doublereal *x, integer *ldx, doublereal *taup, doublereal *tauq, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, logical *
	inul, integer *info, ftnlen jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, l, p, im, ioff, ncol;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwork[1], first;
    static logical wantu, wantv, ljobua, ljobva;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ljobus, ljobvs;


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

/*     To apply the Householder transformations Pj stored in factored */
/*     form into the columns of the array X, to the desired columns of */
/*     the matrix U by premultiplication, and/or the Householder */
/*     transformations Qj stored in factored form into the rows of the */
/*     array X, to the desired columns of the matrix V by */
/*     premultiplication. The Householder transformations Pj and Qj */
/*     are stored as produced by LAPACK Library routine DGEBRD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Specifies whether to transform the columns in U as */
/*             follows: */
/*             = 'N':  Do not transform the columns in U; */
/*             = 'A':  Transform the columns in U (U has M columns); */
/*             = 'S':  Transform the columns in U (U has min(M,N) */
/*                     columns). */

/*     JOBV    CHARACTER*1 */
/*             Specifies whether to transform the columns in V as */
/*             follows: */
/*             = 'N':  Do not transform the columns in V; */
/*             = 'A':  Transform the columns in V (V has N columns); */
/*             = 'S':  Transform the columns in V (V has min(M,N) */
/*                     columns). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix X.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix X.  N >= 0. */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading M-by-N part contains in the columns of its */
/*             lower triangle the Householder transformations Pj, and */
/*             in the rows of its upper triangle the Householder */
/*             transformations Qj in factored form. */
/*             X is modified by the routine but restored on exit. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.   LDX >= MAX(1,M). */

/*     TAUP    (input) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*             The scalar factors of the Householder transformations Pj. */

/*     TAUQ    (input) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*             The scalar factors of the Householder transformations Qj. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, U contains the M-by-M (if JOBU = 'A') or */
/*             M-by-min(M,N) (if JOBU = 'S') matrix U. */
/*             On exit, the Householder transformations Pj have been */
/*             applied to each column i of U corresponding to a parameter */
/*             INUL(i) = .TRUE. */
/*             NOTE that U is not referenced if JOBU = 'N'. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,M), if JOBU = 'A' or JOBU = 'S'; */
/*             LDU >= 1,        if JOBU = 'N'. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             On entry, V contains the N-by-N (if JOBV = 'A') or */
/*             N-by-min(M,N) (if JOBV = 'S') matrix V. */
/*             On exit, the Householder transformations Qj have been */
/*             applied to each column i of V corresponding to a parameter */
/*             INUL(i) = .TRUE. */
/*             NOTE that V is not referenced if JOBV = 'N'. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M), if JOBV = 'A' or JOBV = 'S'; */
/*             LDV >= 1,        if JOBV = 'N'. */

/*     INUL    (input) LOGICAL array, dimension (MAX(M,N)) */
/*             INUL(i) = .TRUE. if the i-th column of U and/or V is to be */
/*             transformed, and INUL(i) = .FALSE., otherwise. */
/*             (1 <= i <= MAX(M,N)). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Householder transformations Pj or Qj are applied to the */
/*     columns of U or V indexed by I for which INUL(I) = .TRUE.. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04PZ by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bidiagonalization, orthogonal transformation, singular subspace, */
/*     singular value decomposition. */

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

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --taup;
    --tauq;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --inul;

    /* Function Body */
    *info = 0;
    ljobua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
    ljobus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
    ljobva = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
    ljobvs = lsame_(jobv, "S", (ftnlen)1, (ftnlen)1);
    wantu = ljobua || ljobus;
    wantv = ljobva || ljobvs;

/*     Test the input scalar arguments. */

    if (! wantu && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! wantv && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldx < max(1,*m)) {
	*info = -6;
    } else if (wantu && *ldu < max(1,*m) || ! wantu && *ldu < 1) {
	*info = -10;
    } else if (wantv && *ldv < max(1,*n) || ! wantv && *ldv < 1) {
	*info = -12;
    }

    if (*info != 0) {

/*        Error return */

	i__1 = -(*info);
	xerbla_("MB04XY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    p = min(*m,*n);
    if (p == 0) {
	return 0;
    }

    if (*m < *n) {
	ioff = 1;
    } else {
	ioff = 0;
    }

/*     Apply the Householder transformations Pj onto the desired */
/*     columns of U. */

/* Computing MIN */
    i__1 = *m - 1;
    im = min(i__1,*n);
    if (wantu && im > 0) {
	if (ljobua) {
	    ncol = *m;
	} else {
	    ncol = p;
	}

	i__1 = ncol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inul[i__]) {

		for (l = im; l >= 1; --l) {
		    if (taup[l] != 0.) {
			first = x[l + ioff + l * x_dim1];
			x[l + ioff + l * x_dim1] = 1.;
			i__2 = *m - l + 1 - ioff;
			dlarf_("Left", &i__2, &c__1, &x[l + ioff + l * x_dim1]
				, &c__1, &taup[l], &u[l + ioff + i__ * u_dim1]
				, ldu, dwork, (ftnlen)4);
			x[l + ioff + l * x_dim1] = first;
		    }
/* L20: */
		}

	    }
/* L40: */
	}

    }

/*     Apply the Householder transformations Qj onto the desired columns */
/*     of V. */

/* Computing MIN */
    i__1 = *n - 1;
    im = min(i__1,*m);
    if (wantv && im > 0) {
	if (ljobva) {
	    ncol = *n;
	} else {
	    ncol = p;
	}

	i__1 = ncol;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inul[i__]) {

		for (l = im; l >= 1; --l) {
		    if (tauq[l] != 0.) {
			first = x[l + (l + 1 - ioff) * x_dim1];
			x[l + (l + 1 - ioff) * x_dim1] = 1.;
			i__2 = *n - l + ioff;
			dlarf_("Left", &i__2, &c__1, &x[l + (l + 1 - ioff) * 
				x_dim1], ldx, &tauq[l], &v[l + 1 - ioff + i__ 
				* v_dim1], ldv, dwork, (ftnlen)4);
			x[l + (l + 1 - ioff) * x_dim1] = first;
		    }
/* L60: */
		}

	    }
/* L80: */
	}

    }

    return 0;
/* *** Last line of MB04XY *** */
} /* mb04xy_ */

