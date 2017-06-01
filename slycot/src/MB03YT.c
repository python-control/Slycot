/* MB03YT.f -- translated by f2c (version 20100827).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03yt_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal r__, t, h1, h2, h3, wi, qq, rr, wr1, wr2, ulp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlag2_(
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal anorm, bnorm, scale1, scale2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


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

/*     To compute the periodic Schur factorization of a real 2-by-2 */
/*     matrix pair (A,B) where B is upper triangular. This routine */
/*     computes orthogonal (rotation) matrices given by CSL, SNL and CSR, */
/*     SNR such that */

/*     1) if the pair (A,B) has two real eigenvalues, then */

/*        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/*        [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */

/*        [ b11 b12 ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ] */
/*        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ], */

/*     2) if the pair (A,B) has a pair of complex conjugate eigenvalues, */
/*        then */

/*        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/*        [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */

/*        [ b11  0  ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ] */
/*        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ]. */

/*     This is a modified version of the LAPACK routine DLAGV2 for */
/*     computing the real, generalized Schur decomposition of a */
/*     two-by-two matrix pencil. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,2) */
/*             On entry, the leading 2-by-2 part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading 2-by-2 part of this array contains */
/*             the matrix A of the pair in periodic Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= 2. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,2) */
/*             On entry, the leading 2-by-2 part of this array must */
/*             contain the upper triangular matrix B. */
/*             On exit, the leading 2-by-2 part of this array contains */
/*             the matrix B of the pair in periodic Schur form. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= 2. */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (2) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (2) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (2) */
/*             (ALPHAR(k)+i*ALPHAI(k))*BETA(k) are the eigenvalues of the */
/*             pair (A,B), k=1,2, i = sqrt(-1). ALPHAI(1) >= 0. */

/*     CSL     (output) DOUBLE PRECISION */
/*             The cosine of the first rotation matrix. */

/*     SNL     (output) DOUBLE PRECISION */
/*             The sine of the first rotation matrix. */

/*     CSR     (output) DOUBLE PRECISION */
/*             The cosine of the second rotation matrix. */

/*     SNR     (output) DOUBLE PRECISION */
/*             The sine of the second rotation matrix. */

/*     REFERENCES */

/*     [1] Van Loan, C. */
/*         Generalized Singular Values with Algorithms and Applications. */
/*         Ph. D. Thesis, University of Michigan, 1973. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPV2). */
/*     V. Sima, July 2008, May 2009. */

/*     KEYWORDS */

/*     Eigenvalue, periodic Schur form */

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
    --alphar;
    --alphai;
    --beta;

    /* Function Body */
    safmin = dlamch_("S", (ftnlen)1);
    ulp = dlamch_("P", (ftnlen)1);

/*     Scale A. */

/* Computing MAX */
    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[a_dim1 + 2], abs(
	    d__2)), d__6 = (d__3 = a[(a_dim1 << 1) + 1], abs(d__3)) + (d__4 = 
	    a[(a_dim1 << 1) + 2], abs(d__4)), d__5 = max(d__5,d__6);
    anorm = max(d__5,safmin);
    a[a_dim1 + 1] /= anorm;
    a[(a_dim1 << 1) + 1] /= anorm;
    a[a_dim1 + 2] /= anorm;
    a[(a_dim1 << 1) + 2] /= anorm;

/*     Scale B. */

/* Computing MAX */
    d__4 = (d__3 = b[b_dim1 + 1], abs(d__3)), d__5 = (d__1 = b[(b_dim1 << 1) 
	    + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1) + 2], abs(d__2)), d__4 
	    = max(d__4,d__5);
    bnorm = max(d__4,safmin);
    b[b_dim1 + 1] /= bnorm;
    b[(b_dim1 << 1) + 1] /= bnorm;
    b[(b_dim1 << 1) + 2] /= bnorm;

/*     Check if A can be deflated. */

    if ((d__1 = a[a_dim1 + 2], abs(d__1)) <= ulp) {
	*csl = 1.;
	*snl = 0.;
	*csr = 1.;
	*snr = 0.;
	wi = 0.;
	a[a_dim1 + 2] = 0.;
	b[b_dim1 + 2] = 0.;

/*     Check if B is singular. */

    } else if ((d__1 = b[b_dim1 + 1], abs(d__1)) <= ulp) {
	dlartg_(&a[(a_dim1 << 1) + 2], &a[a_dim1 + 2], csr, snr, &t);
	*snr = -(*snr);
	drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr,
		 snr);
	drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csr, snr);
	*csl = 1.;
	*snl = 0.;
	wi = 0.;
	a[a_dim1 + 2] = 0.;
	b[b_dim1 + 1] = 0.;
	b[b_dim1 + 2] = 0.;
    } else if ((d__1 = b[(b_dim1 << 1) + 2], abs(d__1)) <= ulp) {
	dlartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);
	*csr = 1.;
	*snr = 0.;
	wi = 0.;
	drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
	drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csl,
		 snl);
	a[a_dim1 + 2] = 0.;
	b[b_dim1 + 2] = 0.;
	b[(b_dim1 << 1) + 2] = 0.;
    } else {

/*        B is nonsingular, first compute the eigenvalues of A / adj(B). */

	r__ = b[b_dim1 + 1];
	b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
	b[(b_dim1 << 1) + 2] = r__;
	b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];
	dlag2_(&a[a_offset], lda, &b[b_offset], ldb, &safmin, &scale1, &
		scale2, &wr1, &wr2, &wi);

	if (wi == 0.) {

/*           Two real eigenvalues, compute s*A-w*B. */

	    h1 = scale1 * a[a_dim1 + 1] - wr1 * b[b_dim1 + 1];
	    h2 = scale1 * a[(a_dim1 << 1) + 1] - wr1 * b[(b_dim1 << 1) + 1];
	    h3 = scale1 * a[(a_dim1 << 1) + 2] - wr1 * b[(b_dim1 << 1) + 2];

	    rr = dlapy2_(&h1, &h2);
	    d__1 = scale1 * a[a_dim1 + 2];
	    qq = dlapy2_(&d__1, &h3);

	    if (rr > qq) {

/*              Find right rotation matrix to zero 1,1 element of */
/*              (sA - wB). */

		dlartg_(&h2, &h1, csr, snr, &t);

	    } else {

/*              Find right rotation matrix to zero 2,1 element of */
/*              (sA - wB). */

		d__1 = scale1 * a[a_dim1 + 2];
		dlartg_(&h3, &d__1, csr, snr, &t);

	    }

	    *snr = -(*snr);
	    drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
	    drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csr, snr);

/*           Compute inf norms of A and B. */

/* Computing MAX */
	    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[(a_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = a[a_dim1 + 2], abs(d__3)
		    ) + (d__4 = a[(a_dim1 << 1) + 2], abs(d__4));
	    h1 = max(d__5,d__6);
/* Computing MAX */
	    d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = b[b_dim1 + 2], abs(d__3)
		    ) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
	    h2 = max(d__5,d__6);

	    if (scale1 * h1 >= abs(wr1) * h2) {

/*              Find left rotation matrix Q to zero out B(2,1). */

		dlartg_(&b[b_dim1 + 1], &b[b_dim1 + 2], csl, snl, &r__);

	    } else {

/*              Find left rotation matrix Q to zero out A(2,1). */

		dlartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);

	    }

	    drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
	    drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);

	    a[a_dim1 + 2] = 0.;
	    b[b_dim1 + 2] = 0.;

/*           Re-adjoint B. */

	    r__ = b[b_dim1 + 1];
	    b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
	    b[(b_dim1 << 1) + 2] = r__;
	    b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];

	} else {

/*           A pair of complex conjugate eigenvalues: */
/*           first compute the SVD of the matrix adj(B). */

	    r__ = b[b_dim1 + 1];
	    b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
	    b[(b_dim1 << 1) + 2] = r__;
	    b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];
	    dlasv2_(&b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], &b[(b_dim1 << 1) + 
		    2], &r__, &t, snl, csl, snr, csr);

/*           Form (A,B) := Q(A,adj(B))Z' where Q is left rotation matrix */
/*           and Z is right rotation matrix computed from DLASV2. */

	    drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
	    drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csr, snr);
	    drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
	    drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csl, snl);

	    b[b_dim1 + 2] = 0.;
	    b[(b_dim1 << 1) + 1] = 0.;
	}

    }

/*     Unscaling */

    r__ = b[b_dim1 + 1];
    t = b[(b_dim1 << 1) + 2];
    a[a_dim1 + 1] = anorm * a[a_dim1 + 1];
    a[a_dim1 + 2] = anorm * a[a_dim1 + 2];
    a[(a_dim1 << 1) + 1] = anorm * a[(a_dim1 << 1) + 1];
    a[(a_dim1 << 1) + 2] = anorm * a[(a_dim1 << 1) + 2];
    b[b_dim1 + 1] = bnorm * b[b_dim1 + 1];
    b[b_dim1 + 2] = bnorm * b[b_dim1 + 2];
    b[(b_dim1 << 1) + 1] = bnorm * b[(b_dim1 << 1) + 1];
    b[(b_dim1 << 1) + 2] = bnorm * b[(b_dim1 << 1) + 2];

    if (wi == 0.) {
	alphar[1] = a[a_dim1 + 1];
	alphar[2] = a[(a_dim1 << 1) + 2];
	alphai[1] = 0.;
	alphai[2] = 0.;
	beta[1] = b[b_dim1 + 1];
	beta[2] = b[(b_dim1 << 1) + 2];
    } else {
	wr1 = anorm * wr1;
	wi = anorm * wi;
	if (abs(wr1) > 1. || wi > 1.) {
	    wr1 *= r__;
	    wi *= r__;
	    r__ = 1.;
	}
	if (abs(wr1) > 1. || abs(wi) > 1.) {
	    wr1 *= t;
	    wi *= t;
	    t = 1.;
	}
	alphar[1] = wr1 / scale1 * r__ * t;
	alphai[1] = (d__1 = wi / scale1 * r__ * t, abs(d__1));
	alphar[2] = alphar[1];
	alphai[2] = -alphai[1];
	beta[1] = bnorm;
	beta[2] = bnorm;
    }
    return 0;
/* *** Last line of MB03YT *** */
} /* mb03yt_ */

