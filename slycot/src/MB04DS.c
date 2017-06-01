/* MB04DS.f -- translated by f2c (version 20100827).
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
static doublereal c_b19 = -1.;

/* Subroutine */ int mb04ds_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, integer *ilo, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j;
    static doublereal r__, s;
    static integer ic;
    static doublereal maxc;
    static logical conv;
    static doublereal maxr;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical lperm;
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal sclfac;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer iloold;


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

/*     To balance a real skew-Hamiltonian matrix */

/*                   [  A   G  ] */
/*              S =  [       T ] , */
/*                   [  Q   A  ] */

/*     where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric */
/*     matrices. This involves, first, permuting S by a symplectic */
/*     similarity transformation to isolate eigenvalues in the first */
/*     1:ILO-1 elements on the diagonal of A; and second, applying a */
/*     diagonal similarity transformation to rows and columns */
/*     ILO:2*N-ILO+1 to make the rows and columns as close in 1-norm */
/*     as possible. Both steps are optional. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the operations to be performed on S: */
/*             = 'N':  none, set ILO = 1, SCALE(I) = 1.0, I = 1 .. N; */
/*             = 'P':  permute only; */
/*             = 'S':  scale only; */
/*             = 'B':  both permute and scale. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix A of the balanced skew-Hamiltonian. In */
/*             particular, the lower triangular part of the first ILO-1 */
/*             columns of A is zero. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the strictly lower triangular part */
/*             of the matrix Q and in columns 2:N+1 the strictly upper */
/*             triangular part of the matrix G. The parts containing the */
/*             diagonal and the first supdiagonal of this array are not */
/*             referenced. */
/*             On exit, the leading N-by-N+1 part of this array contains */
/*             the strictly lower and strictly upper triangular parts of */
/*             the matrices Q and G, respectively, of the balanced */
/*             skew-Hamiltonian. In particular, the strictly lower */
/*             triangular part of the first ILO-1 columns of QG is zero. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     ILO     (output) INTEGER */
/*             ILO-1 is the number of deflated eigenvalues in the */
/*             balanced skew-Hamiltonian matrix. */

/*     SCALE   (output) DOUBLE PRECISION array of dimension (N) */
/*             Details of the permutations and scaling factors applied to */
/*             S.  For j = 1,...,ILO-1 let P(j) = SCALE(j). If P(j) <= N, */
/*             then rows and columns P(j) and P(j)+N are interchanged */
/*             with rows and columns j and j+N, respectively. If */
/*             P(j) > N, then row and column P(j)-N are interchanged with */
/*             row and column j+N by a generalized symplectic */
/*             permutation. For j = ILO,...,N the j-th element of SCALE */
/*             contains the factor of the scaling applied to row and */
/*             column j. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DSHBAL). */

/*     KEYWORDS */

/*     Balancing, skew-Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    --scale;

    /* Function Body */
    *info = 0;
    lperm = lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
    lscal = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);

    if (! lperm && ! lscal && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldqg < max(1,*n)) {
	*info = -6;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04DS", &i__1, (ftnlen)6);
	return 0;
    }

    *ilo = 1;

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }
    if (! lperm && ! lscal) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scale[i__] = 1.;
/* L10: */
	}
	return 0;
    }

/*     Permutations to isolate eigenvalues if possible. */

    if (lperm) {
	iloold = 0;
/*        WHILE ( ILO.NE.ILOOLD ) */
L20:
	if (*ilo != iloold) {
	    iloold = *ilo;

/*           Scan columns ILO .. N. */

	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
L30:
	    if (i__ <= *n && *ilo == iloold) {
		i__1 = i__ - 1;
		for (j = *ilo; j <= i__1; ++j) {
		    if (a[j + i__ * a_dim1] != 0.) {
			++i__;
			goto L30;
		    }
/* L40: */
		}
		i__1 = *n;
		for (j = i__ + 1; j <= i__1; ++j) {
		    if (a[j + i__ * a_dim1] != 0.) {
			++i__;
			goto L30;
		    }
/* L50: */
		}
		i__1 = i__ - 1;
		for (j = *ilo; j <= i__1; ++j) {
		    if (qg[i__ + j * qg_dim1] != 0.) {
			++i__;
			goto L30;
		    }
/* L60: */
		}
		i__1 = *n;
		for (j = i__ + 1; j <= i__1; ++j) {
		    if (qg[j + i__ * qg_dim1] != 0.) {
			++i__;
			goto L30;
		    }
/* L70: */
		}

/*              Exchange columns/rows ILO <-> I. */

		scale[*ilo] = (doublereal) i__;
		if (*ilo != i__) {

		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
		    i__1 = *n - *ilo + 1;
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

		    if (i__ < *n) {
			i__1 = *n - i__;
			dswap_(&i__1, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
				qg[i__ + 1 + *ilo * qg_dim1], &c__1);
		    }
		    if (i__ > *ilo + 1) {
			i__1 = i__ - *ilo - 1;
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + *ilo * qg_dim1], 
				&c__1);
			i__1 = i__ - *ilo - 1;
			dswap_(&i__1, &qg[*ilo + 1 + *ilo * qg_dim1], &c__1, &
				qg[i__ + (*ilo + 1) * qg_dim1], ldqg);
		    }

		    i__1 = *ilo - 1;
		    dswap_(&i__1, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
		    if (*n > i__) {
			i__1 = *n - i__;
			dswap_(&i__1, &qg[i__ + (i__ + 2) * qg_dim1], ldqg, &
				qg[*ilo + (i__ + 2) * qg_dim1], ldqg);
		    }
		    if (i__ > *ilo + 1) {
			i__1 = i__ - *ilo - 1;
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + (i__ + 1) * 
				qg_dim1], &c__1);
			i__1 = i__ - *ilo - 1;
			dswap_(&i__1, &qg[*ilo + (*ilo + 2) * qg_dim1], ldqg, 
				&qg[*ilo + 1 + (i__ + 1) * qg_dim1], &c__1);
		    }
		    i__1 = i__ - *ilo;
		    dscal_(&i__1, &c_b19, &qg[*ilo + (i__ + 1) * qg_dim1], &
			    c__1);
		}
		++(*ilo);
	    }
/*           END WHILE 30 */

/*           Scan columns N+ILO .. 2*N. */

	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
L80:
	    if (i__ <= *n && *ilo == iloold) {
		i__1 = i__ - 1;
		for (j = *ilo; j <= i__1; ++j) {
		    if (a[i__ + j * a_dim1] != 0.) {
			++i__;
			goto L80;
		    }
/* L90: */
		}
		i__1 = *n;
		for (j = i__ + 1; j <= i__1; ++j) {
		    if (a[i__ + j * a_dim1] != 0.) {
			++i__;
			goto L80;
		    }
/* L100: */
		}
		i__1 = i__ - 1;
		for (j = *ilo; j <= i__1; ++j) {
		    if (qg[j + (i__ + 1) * qg_dim1] != 0.) {
			++i__;
			goto L80;
		    }
/* L110: */
		}
		i__1 = *n;
		for (j = i__ + 1; j <= i__1; ++j) {
		    if (qg[i__ + (j + 1) * qg_dim1] != 0.) {
			++i__;
			goto L80;
		    }
/* L120: */
		}
		scale[*ilo] = (doublereal) (*n + i__);

/*              Exchange columns/rows I <-> I+N with a symplectic */
/*              generalized permutation. */

		i__1 = i__ - *ilo;
		dswap_(&i__1, &a[i__ + *ilo * a_dim1], lda, &qg[i__ + *ilo * 
			qg_dim1], ldqg);
		i__1 = i__ - *ilo;
		dscal_(&i__1, &c_b19, &a[i__ + *ilo * a_dim1], lda);
		i__1 = *n - i__;
		dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &qg[i__ + 1 
			+ i__ * qg_dim1], &c__1);
		i__1 = *n - i__;
		dscal_(&i__1, &c_b19, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
		i__1 = i__ - 1;
		dswap_(&i__1, &a[i__ * a_dim1 + 1], &c__1, &qg[(i__ + 1) * 
			qg_dim1 + 1], &c__1);
		i__1 = i__ - 1;
		dscal_(&i__1, &c_b19, &a[i__ * a_dim1 + 1], &c__1);
		i__1 = *n - i__;
		dscal_(&i__1, &c_b19, &a[i__ + 1 + i__ * a_dim1], &c__1);
		i__1 = *n - i__;
		dswap_(&i__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &qg[i__ + (
			i__ + 2) * qg_dim1], ldqg);

/*              Exchange columns/rows ILO <-> I. */

		if (*ilo != i__) {

		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
		    i__1 = *n - *ilo + 1;
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

		    if (i__ < *n) {
			i__1 = *n - i__;
			dswap_(&i__1, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
				qg[i__ + 1 + *ilo * qg_dim1], &c__1);
		    }
		    if (i__ > *ilo + 1) {
			i__1 = i__ - *ilo - 1;
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + *ilo * qg_dim1], 
				&c__1);
			i__1 = i__ - *ilo - 1;
			dswap_(&i__1, &qg[*ilo + 1 + *ilo * qg_dim1], &c__1, &
				qg[i__ + (*ilo + 1) * qg_dim1], ldqg);
		    }

		    i__1 = *ilo - 1;
		    dswap_(&i__1, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
		    if (*n > i__) {
			i__1 = *n - i__;
			dswap_(&i__1, &qg[i__ + (i__ + 2) * qg_dim1], ldqg, &
				qg[*ilo + (i__ + 2) * qg_dim1], ldqg);
		    }
		    if (i__ > *ilo + 1) {
			i__1 = i__ - *ilo - 1;
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + (i__ + 1) * 
				qg_dim1], &c__1);
			i__1 = i__ - *ilo - 1;
			dswap_(&i__1, &qg[*ilo + (*ilo + 2) * qg_dim1], ldqg, 
				&qg[*ilo + 1 + (i__ + 1) * qg_dim1], &c__1);
		    }
		    i__1 = i__ - *ilo;
		    dscal_(&i__1, &c_b19, &qg[*ilo + (i__ + 1) * qg_dim1], &
			    c__1);
		}
		++(*ilo);
	    }
/*           END WHILE 80 */
	    goto L20;
	}
/*        END WHILE 20 */
    }

    i__1 = *n;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	scale[i__] = 1.;
/* L130: */
    }

/*     Scale to reduce the 1-norm of the remaining blocks. */

    if (lscal) {
	sclfac = dlamch_("B", (ftnlen)1);
	sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
	sfmax1 = 1. / sfmin1;
	sfmin2 = sfmin1 * sclfac;
	sfmax2 = 1. / sfmin2;

/*        Scale the rows and columns one at a time to minimize the */
/*        1-norm of the skew-Hamiltonian submatrix. */
/*        Stop when the 1-norm is very roughly minimal. */

L140:
	conv = TRUE_;
	i__1 = *n;
	for (i__ = *ilo; i__ <= i__1; ++i__) {

/*              Compute 1-norm of row and column I without diagonal */
/*              elements. */

	    i__2 = i__ - *ilo;
	    i__3 = *n - i__;
	    i__4 = i__ - *ilo;
	    i__5 = *n - i__;
	    r__ = dasum_(&i__2, &a[i__ + *ilo * a_dim1], lda) + dasum_(&i__3, 
		    &a[i__ + (i__ + 1) * a_dim1], lda) + dasum_(&i__4, &qg[*
		    ilo + (i__ + 1) * qg_dim1], &c__1) + dasum_(&i__5, &qg[
		    i__ + (i__ + 2) * qg_dim1], ldqg);
	    i__2 = i__ - *ilo;
	    i__3 = *n - i__;
	    i__4 = i__ - *ilo;
	    i__5 = *n - i__;
	    c__ = dasum_(&i__2, &a[*ilo + i__ * a_dim1], &c__1) + dasum_(&
		    i__3, &a[i__ + 1 + i__ * a_dim1], &c__1) + dasum_(&i__4, &
		    qg[i__ + *ilo * qg_dim1], ldqg) + dasum_(&i__5, &qg[i__ + 
		    1 + i__ * qg_dim1], &c__1);

/*              Compute inf-norms of row and column I. */

	    i__2 = *n - *ilo + 1;
	    ic = idamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
	    maxr = (d__1 = a[i__ + (ic + *ilo - 1) * a_dim1], abs(d__1));
	    if (i__ > 1) {
		i__2 = i__ - 1;
		ic = idamax_(&i__2, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
/* Computing MAX */
		d__2 = maxr, d__3 = (d__1 = qg[ic + (i__ + 1) * qg_dim1], abs(
			d__1));
		maxr = max(d__2,d__3);
	    }
	    if (*n > i__) {
		i__2 = *n - i__;
		ic = idamax_(&i__2, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
/* Computing MAX */
		d__2 = maxr, d__3 = (d__1 = qg[i__ + (ic + i__ + 1) * qg_dim1]
			, abs(d__1));
		maxr = max(d__2,d__3);
	    }
	    ic = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
	    maxc = (d__1 = a[ic + i__ * a_dim1], abs(d__1));
	    if (i__ > *ilo) {
		i__2 = i__ - *ilo;
		ic = idamax_(&i__2, &qg[i__ + *ilo * qg_dim1], ldqg);
/* Computing MAX */
		d__2 = maxc, d__3 = (d__1 = qg[i__ + (ic + *ilo - 1) * 
			qg_dim1], abs(d__1));
		maxc = max(d__2,d__3);
	    }
	    if (*n > i__) {
		i__2 = *n - i__;
		ic = idamax_(&i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
/* Computing MAX */
		d__2 = maxc, d__3 = (d__1 = qg[ic + i__ + i__ * qg_dim1], abs(
			d__1));
		maxc = max(d__2,d__3);
	    }

	    if (c__ == 0. || r__ == 0.) {
		goto L190;
	    }
	    g = r__ / sclfac;
	    f = 1.;
	    s = c__ + r__;
L150:
/* Computing MAX */
	    d__1 = max(f,c__);
/* Computing MIN */
	    d__2 = min(r__,g);
	    if (c__ >= g || max(d__1,maxc) >= sfmax2 || min(d__2,maxr) <= 
		    sfmin2) {
		goto L160;
	    }
	    f *= sclfac;
	    g /= sclfac;
	    c__ *= sclfac;
	    r__ /= sclfac;
	    maxc *= sclfac;
	    maxr /= sclfac;
	    goto L150;

L160:
	    g = c__ / sclfac;
L170:
/* Computing MIN */
	    d__1 = min(f,c__), d__1 = min(d__1,g);
	    if (g < r__ || max(r__,maxr) >= sfmax2 || min(d__1,maxc) <= 
		    sfmin2) {
		goto L180;
	    }
	    f /= sclfac;
	    g /= sclfac;
	    c__ /= sclfac;
	    r__ *= sclfac;
	    maxc /= sclfac;
	    maxr *= sclfac;
	    goto L170;

L180:

/*              Now balance if necessary. */

	    if (c__ + r__ >= s * .95) {
		goto L190;
	    }
	    if (f < 1. && scale[i__] < 1.) {
		if (f * scale[i__] <= sfmin1) {
		    goto L190;
		}
	    }
	    if (f > 1. && scale[i__] > 1.) {
		if (scale[i__] >= sfmax1 / f) {
		    goto L190;
		}
	    }
	    conv = FALSE_;
	    scale[i__] *= f;
	    i__2 = i__ - *ilo;
	    drscl_(&i__2, &f, &a[i__ + *ilo * a_dim1], lda);
	    i__2 = *n - i__;
	    drscl_(&i__2, &f, &a[i__ + (i__ + 1) * a_dim1], lda);
	    i__2 = i__ - 1;
	    dscal_(&i__2, &f, &a[i__ * a_dim1 + 1], &c__1);
	    i__2 = *n - i__;
	    dscal_(&i__2, &f, &a[i__ + 1 + i__ * a_dim1], &c__1);
	    i__2 = i__ - 1;
	    drscl_(&i__2, &f, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
	    i__2 = *n - i__;
	    drscl_(&i__2, &f, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
	    i__2 = i__ - *ilo;
	    dscal_(&i__2, &f, &qg[i__ + *ilo * qg_dim1], ldqg);
	    i__2 = *n - i__;
	    dscal_(&i__2, &f, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
L190:
	    ;
	}
	if (! conv) {
	    goto L140;
	}
    }
    return 0;
/* *** Last line of MB04DS *** */
} /* mb04ds_ */

