/* SB01MD.f -- translated by f2c (version 20100827).
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
static doublereal c_b23 = 1.;
static doublereal c_b25 = 0.;

/* Subroutine */ int sb01md_(integer *ncont, integer *n, doublereal *a, 
	integer *lda, doublereal *b, doublereal *wr, doublereal *wi, 
	doublereal *z__, integer *ldz, doublereal *g, doublereal *dwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal p, q, r__, s, t, b1;
    static integer ni, ll, nj, nl, im1, lp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dgemv_(char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical compl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer ncont2;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen);


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

/*     To determine the one-dimensional state feedback matrix G of the */
/*     linear time-invariant single-input system */

/*           dX/dt = A * X + B * U, */

/*     where A is an NCONT-by-NCONT matrix and B is an NCONT element */
/*     vector such that the closed-loop system */

/*           dX/dt = (A - B * G) * X */

/*     has desired poles. The system must be preliminarily reduced */
/*     to orthogonal canonical form using the SLICOT Library routine */
/*     AB01MD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NCONT   (input) INTEGER */
/*             The order of the matrix A as produced by SLICOT Library */
/*             routine AB01MD.  NCONT >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix Z.  N >= NCONT. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,NCONT) */
/*             On entry, the leading NCONT-by-NCONT part of this array */
/*             must contain the canonical form of the state dynamics */
/*             matrix A as produced by SLICOT Library routine AB01MD. */
/*             On exit, the leading NCONT-by-NCONT part of this array */
/*             contains the upper quasi-triangular form S of the closed- */
/*             loop system matrix (A - B * G), that is triangular except */
/*             for possible 2-by-2 diagonal blocks. */
/*             (To reconstruct the closed-loop system matrix see */
/*             FURTHER COMMENTS below.) */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,NCONT). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (NCONT) */
/*             On entry, this array must contain the canonical form of */
/*             the input/state vector B as produced by SLICOT Library */
/*             routine AB01MD. */
/*             On exit, this array contains the transformed vector Z * B */
/*             of the closed-loop system. */

/*     WR      (input) DOUBLE PRECISION array, dimension (NCONT) */
/*     WI      (input) DOUBLE PRECISION array, dimension (NCONT) */
/*             These arrays must contain the real and imaginary parts, */
/*             respectively, of the desired poles of the closed-loop */
/*             system. The poles can be unordered, except that complex */
/*             conjugate pairs of poles must appear consecutively. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the orthogonal transformation matrix as produced */
/*             by SLICOT Library routine AB01MD, which reduces the system */
/*             to canonical form. */
/*             On exit, the leading NCONT-by-NCONT part of this array */
/*             contains the orthogonal matrix Z which reduces the closed- */
/*             loop system matrix (A - B * G) to upper quasi-triangular */
/*             form. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,N). */

/*     G       (output) DOUBLE PRECISION array, dimension (NCONT) */
/*             This array contains the one-dimensional state feedback */
/*             matrix G of the original system. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*NCONT) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The method is based on the orthogonal reduction of the closed-loop */
/*     system matrix (A - B * G) to upper quasi-triangular form S whose */
/*     1-by-1 and 2-by-2 diagonal blocks correspond to the desired poles. */
/*     That is, S = Z'*(A - B * G)*Z, where Z is an orthogonal matrix. */

/*     REFERENCES */

/*     [1] Petkov, P. Hr. */
/*         A Computational Algorithm for Pole Assignment of Linear */
/*         Single Input Systems. */
/*         Internal Report 81/2, Control Systems Research Group, School */
/*         of Electronic Engineering and Computer Science, Kingston */
/*         Polytechnic, 1981. */

/*     NUMERICAL ASPECTS */
/*                                   3 */
/*     The algorithm requires 0(NCONT ) operations and is backward */
/*     stable. */

/*     FURTHER COMMENTS */

/*     If required, the closed-loop system matrix (A - B * G) can be */
/*     formed from the matrix product Z * S * Z' (where S and Z are the */
/*     matrices output in arrays A and Z respectively). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB01AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, May 1981. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Closed loop spectrum, closed loop systems, eigenvalue assignment, */
/*     orthogonal canonical form, orthogonal transformation, pole */
/*     placement, Schur form. */

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
    --b;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --g;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*ncont < 0) {
	*info = -1;
    } else if (*n < *ncont) {
	*info = -2;
    } else if (*lda < max(1,*ncont)) {
	*info = -4;
    } else if (*ldz < max(1,*n)) {
	*info = -9;
    }

    if (*info != 0) {

/*        Error return */

	i__1 = -(*info);
	xerbla_("SB01MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*ncont == 0 || *n == 0) {
	return 0;
    }

/*     Return if the system is not complete controllable. */

    if (b[1] == 0.) {
	return 0;
    }

    if (*ncont == 1) {

/*        1-by-1 case. */

	p = a[a_dim1 + 1] - wr[1];
	a[a_dim1 + 1] = wr[1];
	g[1] = p / b[1];
	z__[z_dim1 + 1] = 1.;
	return 0;
    }

/*     General case.  Save the contents of WI in DWORK. */

    ncont2 = *ncont << 1;
    dcopy_(ncont, &wi[1], &c__1, &dwork[ncont2 + 1], &c__1);

    b1 = b[1];
    b[1] = 1.;
    l = 0;
    ll = 0;
L20:
    ++l;
    ++ll;
    compl = dwork[ncont2 + l] != 0.;
    if (l != *ncont) {
	lp1 = l + 1;
	nl = *ncont - l;
	if (ll != 2) {
	    if (compl) {

/*              Compute complex eigenvector. */

		dwork[*ncont] = 1.;
		dwork[ncont2] = 1.;
		p = wr[l];
		t = dwork[ncont2 + l];
		q = t * dwork[ncont2 + lp1];
		dwork[ncont2 + l] = 1.;
		dwork[ncont2 + lp1] = q;

		i__1 = lp1;
		for (i__ = *ncont; i__ >= i__1; --i__) {
		    im1 = i__ - 1;
		    i__2 = *ncont - im1;
		    dwork[im1] = (p * dwork[i__] + q * dwork[*ncont + i__] - 
			    ddot_(&i__2, &a[i__ + i__ * a_dim1], lda, &dwork[
			    i__], &c__1)) / a[i__ + im1 * a_dim1];
		    i__2 = *ncont - im1;
		    dwork[*ncont + im1] = (p * dwork[*ncont + i__] + dwork[
			    i__] - ddot_(&i__2, &a[i__ + i__ * a_dim1], lda, &
			    dwork[*ncont + i__], &c__1)) / a[i__ + im1 * 
			    a_dim1];
/* L40: */
		}

	    } else {

/*              Compute real eigenvector. */

		dwork[*ncont] = 1.;
		p = wr[l];

		i__1 = lp1;
		for (i__ = *ncont; i__ >= i__1; --i__) {
		    im1 = i__ - 1;
		    i__2 = *ncont - im1;
		    dwork[im1] = (p * dwork[i__] - ddot_(&i__2, &a[i__ + i__ *
			     a_dim1], lda, &dwork[i__], &c__1)) / a[i__ + im1 
			    * a_dim1];
/* L60: */
		}

	    }
	}

/*        Transform eigenvector. */

	i__1 = l;
	for (k = *ncont - 1; k >= i__1; --k) {
	    if (ll != 2) {
		r__ = dwork[k];
		s = dwork[k + 1];
	    } else {
		r__ = dwork[*ncont + k];
		s = dwork[*ncont + k + 1];
	    }
	    dlartg_(&r__, &s, &p, &q, &t);
	    dwork[k] = t;
	    if (ll != 2) {
/* Computing MAX */
		i__2 = k - 1;
		nj = max(i__2,l);
	    } else {
		dwork[*ncont + k] = t;
		nj = l - 1;
	    }

/*           Transform  A. */

	    i__2 = *ncont - nj + 1;
	    drot_(&i__2, &a[k + nj * a_dim1], lda, &a[k + 1 + nj * a_dim1], 
		    lda, &p, &q);

	    if (compl && ll == 1) {
		ni = *ncont;
	    } else {
/* Computing MIN */
		i__2 = k + 2;
		ni = min(i__2,*ncont);
	    }
	    drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
		    c__1, &p, &q);

	    if (k == l) {

/*              Transform  B. */

		t = b[k];
		b[k] = p * t;
		b[k + 1] = -q * t;
	    }

/*           Accumulate transformations. */

	    drot_(ncont, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * z_dim1 + 
		    1], &c__1, &p, &q);

	    if (compl && ll != 2) {
		t = dwork[*ncont + k];
		dwork[*ncont + k] = p * t + q * dwork[*ncont + k + 1];
		dwork[*ncont + k + 1] = p * dwork[*ncont + k + 1] - q * t;
	    }
/* L80: */
	}

    }

    if (! compl) {

/*        Find one element of  G. */

	k = l;
	r__ = b[l];
	if (l != *ncont) {
	    if ((d__1 = b[lp1], abs(d__1)) > (d__2 = b[l], abs(d__2))) {
		k = lp1;
		r__ = b[lp1];
	    }
	}
	p = a[k + l * a_dim1];
	if (k == l) {
	    p -= wr[l];
	}
	p /= r__;

	d__1 = -p;
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[l * a_dim1 + 1], &c__1);

	g[l] = p / b1;
	if (l != *ncont) {
	    ll = 0;
	    goto L20;
	}
    } else if (ll == 1) {
	goto L20;
    } else {

/*        Find two elements of  G. */

	k = l;
	r__ = b[l];
	if (l != *ncont) {
	    if ((d__1 = b[lp1], abs(d__1)) > (d__2 = b[l], abs(d__2))) {
		k = lp1;
		r__ = b[lp1];
	    }
	}
	p = a[k + (l - 1) * a_dim1];
	q = a[k + l * a_dim1];
	if (k == l) {
	    p -= dwork[*ncont + l] / dwork[l - 1] * dwork[ncont2 + l];
	    q = q - wr[l] + dwork[*ncont + l - 1] / dwork[l - 1] * dwork[
		    ncont2 + l];
	}
	p /= r__;
	q /= r__;

	d__1 = -p;
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[(l - 1) * a_dim1 + 1], &c__1);
	d__1 = -q;
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[l * a_dim1 + 1], &c__1);

	g[l - 1] = p / b1;
	g[l] = q / b1;
	if (l != *ncont) {
	    ll = 0;
	    goto L20;
	}
    }

/*     Transform  G. */

    dgemv_("No transpose", ncont, ncont, &c_b23, &z__[z_offset], ldz, &g[1], &
	    c__1, &c_b25, &dwork[1], &c__1, (ftnlen)12);
    dcopy_(ncont, &dwork[1], &c__1, &g[1], &c__1);
    dscal_(ncont, &b1, &b[1], &c__1);

/*     Annihilate A after the first subdiagonal. */

    if (*ncont > 2) {
	i__1 = *ncont - 2;
	i__2 = *ncont - 2;
	dlaset_("Lower", &i__1, &i__2, &c_b25, &c_b25, &a[a_dim1 + 3], lda, (
		ftnlen)5);
    }

    return 0;
/* *** Last line of SB01MD *** */
} /* sb01md_ */

