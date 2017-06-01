/* SB01BY.f -- translated by f2c (version 20100827).
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
static doublereal c_b4 = 0.;
static integer c__2 = 2;

/* Subroutine */ int sb01by_(integer *n, integer *m, doublereal *s, 
	doublereal *p, doublereal *a, doublereal *b, doublereal *f, 
	doublereal *tol, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublereal r__, x, y, z__, b1, b2, c0, c1, c3, c4, b21, c11, c12, 
	    c21, c22, cs, s12, cu, cv, s21;
    static integer ir;
    static doublereal rn, sn, wi, su, sv, wr, dc0, dc2, dc3, wi1, wr1, sig, 
	    tau1, tau2, absr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal diffr;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasv2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlatzm_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen);


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

/*     To solve an N-by-N pole placement problem for the simple cases */
/*     N = 1 or N = 2: given the N-by-N matrix A and N-by-M matrix B, */
/*     construct an M-by-N matrix F such that A + B*F has prescribed */
/*     eigenvalues. These eigenvalues are specified by their sum S and */
/*     product P (if N = 2). The resulting F has minimum Frobenius norm. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and also the number of rows of */
/*             the matrix B and the number of columns of the matrix F. */
/*             N is either 1, if a single real eigenvalue is prescribed */
/*             or 2, if a complex conjugate pair or a set of two real */
/*             eigenvalues are prescribed. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B and also the number */
/*             of rows of the matrix F.  M >= 1. */

/*     S       (input) DOUBLE PRECISION */
/*             The sum of the prescribed eigenvalues if N = 2 or the */
/*             value of prescribed eigenvalue if N = 1. */

/*     P       (input) DOUBLE PRECISION */
/*             The product of the prescribed eigenvalues if N = 2. */
/*             Not referenced if N = 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N,N) */
/*             On entry, this array must contain the N-by-N state */
/*             dynamics matrix whose eigenvalues have to be moved to */
/*             prescribed locations. */
/*             On exit, this array contains no useful information. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N,M) */
/*             On entry, this array must contain the N-by-M input/state */
/*             matrix B. */
/*             On exit, this array contains no useful information. */

/*     F       (output) DOUBLE PRECISION array, dimension (M,N) */
/*             The state feedback matrix F which assigns one pole or two */
/*             poles of the closed-loop matrix A + B*F. */
/*             If N = 2 and the pair (A,B) is not controllable */
/*             (INFO = 1), then F(1,1) and F(1,2) contain the elements of */
/*             an orthogonal rotation which can be used to remove the */
/*             uncontrollable part of the pair (A,B). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of A */
/*             and B are considered zero (used for controllability test). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if uncontrollability of the pair (A,B) is detected. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine SB01BY. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     May  2003, A. Varga, German Aerospace Center. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue assignment, feedback control, pole */
/*     placement, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

    /* Parameter adjustments */
    b_dim1 = *n;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *n;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    f_dim1 = *m;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n == 1) {

/*        The case N = 1. */

	if (*m > 1) {
	    dlarfg_(m, &b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], n, &tau1);
	}
	b1 = b[b_dim1 + 1];
	if (abs(b1) <= *tol) {

/*           The pair (A,B) is uncontrollable. */

	    *info = 1;
	    return 0;
	}

	f[f_dim1 + 1] = (*s - a[a_dim1 + 1]) / b1;
	if (*m > 1) {
	    i__1 = *m - 1;
	    dlaset_("Full", &i__1, &c__1, &c_b4, &c_b4, &f[f_dim1 + 2], m, (
		    ftnlen)4);
	    dlatzm_("Left", m, n, &b[(b_dim1 << 1) + 1], n, &tau1, &f[f_dim1 
		    + 1], &f[f_dim1 + 2], m, &dwork[1], (ftnlen)4);
	}
	return 0;
    }

/*     In the sequel N = 2. */

/*     Compute the singular value decomposition of B in the form */

/*                    ( V  0 )                ( B1 0  ) */
/*     B = U*( G1 0 )*(      )*H2*H1 ,   G1 = (       ), */
/*                    ( 0  I )                ( 0  B2 ) */

/*               ( CU   SU )          ( CV   SV ) */
/*     where U = (         )  and V = (         )  are orthogonal */
/*               (-SU   CU )          (-SV   CV ) */

/*     rotations and H1 and H2 are elementary Householder reflectors. */
/*     ABS(B1) and ABS(B2) are the singular values of matrix B, */
/*     with ABS(B1) >= ABS(B2). */

/*     Reduce first B to the lower bidiagonal form  ( B1  0  ... 0 ). */
/*                                                  ( B21 B2 ... 0 ) */
    if (*m == 1) {

/*        Initialization for the case M = 1; no reduction required. */

	b1 = b[b_dim1 + 1];
	b21 = b[b_dim1 + 2];
	b2 = 0.;
    } else {

/*        Postmultiply B with elementary Householder reflectors H1 */
/*        and H2. */

	dlarfg_(m, &b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], n, &tau1);
	i__1 = *n - 1;
	dlatzm_("Right", &i__1, m, &b[(b_dim1 << 1) + 1], n, &tau1, &b[b_dim1 
		+ 2], &b[(b_dim1 << 1) + 2], n, &dwork[1], (ftnlen)5);
	b1 = b[b_dim1 + 1];
	b21 = b[b_dim1 + 2];
	if (*m > 2) {
	    i__1 = *m - 1;
	    dlarfg_(&i__1, &b[(b_dim1 << 1) + 2], &b[b_dim1 * 3 + 2], n, &
		    tau2);
	}
	b2 = b[(b_dim1 << 1) + 2];
    }

/*     Reduce B to a diagonal form by premultiplying and postmultiplying */
/*     it with orthogonal rotations U and V, respectively, and order the */
/*     diagonal elements to have decreasing magnitudes. */
/*     Note: B2 has been set to zero if M = 1. Thus in the following */
/*     computations the case M = 1 need not to be distinguished. */
/*     Note also that LAPACK routine DLASV2 assumes an upper triangular */
/*     matrix, so the results should be adapted. */

    dlasv2_(&b1, &b21, &b2, &x, &y, &su, &cu, &sv, &cv);
    su = -su;
    b1 = y;
    b2 = x;

/*     Compute  A1 = U'*A*U. */

    drot_(&c__2, &a[a_dim1 + 2], &c__2, &a[a_dim1 + 1], &c__2, &cu, &su);
    drot_(&c__2, &a[(a_dim1 << 1) + 1], &c__1, &a[a_dim1 + 1], &c__1, &cu, &
	    su);

/*     Compute the rank of B and check the controllability of the */
/*     pair (A,B). */

    ir = 0;
    if (abs(b2) > *tol) {
	++ir;
    }
    if (abs(b1) > *tol) {
	++ir;
    }
    if (ir == 0 || ir == 1 && (d__1 = a[a_dim1 + 2], abs(d__1)) <= *tol) {
	f[f_dim1 + 1] = cu;
	f[(f_dim1 << 1) + 1] = -su;

/*        The pair (A,B) is uncontrollable. */

	*info = 1;
	return 0;
    }

/*     Compute F1 which assigns N poles for the reduced pair (A1,G1). */

    x = dlamc3_(&b1, &b2);
    if (x == b1) {

/*        Rank one G1. */

	f[f_dim1 + 1] = (*s - (a[a_dim1 + 1] + a[(a_dim1 << 1) + 2])) / b1;
	f[(f_dim1 << 1) + 1] = -(a[(a_dim1 << 1) + 2] * (a[(a_dim1 << 1) + 2] 
		- *s) + a[a_dim1 + 2] * a[(a_dim1 << 1) + 1] + *p) / a[a_dim1 
		+ 2] / b1;
	if (*m > 1) {
	    f[f_dim1 + 2] = 0.;
	    f[(f_dim1 << 1) + 2] = 0.;
	}
    } else {

/*        Rank two G1. */

	z__ = (*s - (a[a_dim1 + 1] + a[(a_dim1 << 1) + 2])) / (b1 * b1 + b2 * 
		b2);
	f[f_dim1 + 1] = b1 * z__;
	f[(f_dim1 << 1) + 2] = b2 * z__;

/*        Compute an approximation for the minimum norm parameter */
/*        selection. */

	x = a[a_dim1 + 1] + b1 * f[f_dim1 + 1];
	c__ = x * (*s - x) - *p;
	if (c__ >= 0.) {
	    sig = 1.;
	} else {
	    sig = -1.;
	}
	s12 = b1 / b2;
	s21 = b2 / b1;
	c11 = 0.;
	c12 = 1.;
	c21 = sig * s12 * c__;
	c22 = a[(a_dim1 << 1) + 1] - sig * s12 * a[a_dim1 + 2];
	dlanv2_(&c11, &c12, &c21, &c22, &wr, &wi, &wr1, &wi1, &cs, &sn);
	if ((d__1 = wr - a[(a_dim1 << 1) + 1], abs(d__1)) > (d__2 = wr1 - a[(
		a_dim1 << 1) + 1], abs(d__2))) {
	    r__ = wr1;
	} else {
	    r__ = wr;
	}

/*        Perform Newton iteration to solve the equation for minimum. */

	c0 = -c__ * c__;
	c1 = c__ * a[a_dim1 + 2];
	c4 = s21 * s21;
	c3 = -c4 * a[(a_dim1 << 1) + 1];
	dc0 = c1;
	dc2 = c3 * 3.;
	dc3 = c4 * 4.;

	for (j = 1; j <= 10; ++j) {
	    x = c0 + r__ * (c1 + r__ * r__ * (c3 + r__ * c4));
	    y = dc0 + r__ * r__ * (dc2 + r__ * dc3);
	    if (y == 0.) {
		goto L20;
	    }
	    rn = r__ - x / y;
	    absr = abs(r__);
	    diffr = (d__1 = r__ - rn, abs(d__1));
	    z__ = dlamc3_(&absr, &diffr);
	    if (z__ == absr) {
		goto L20;
	    }
	    r__ = rn;
/* L10: */
	}

L20:
	if (r__ == 0.) {
	    r__ = dlamch_("Epsilon", (ftnlen)7);
	}
	f[(f_dim1 << 1) + 1] = (r__ - a[(a_dim1 << 1) + 1]) / b1;
	f[f_dim1 + 2] = (c__ / r__ - a[a_dim1 + 2]) / b2;
    }

/*     Back-transform F1. Compute first F1*U'. */

    i__1 = min(*m,2);
    drot_(&i__1, &f[f_dim1 + 1], &c__1, &f[(f_dim1 << 1) + 1], &c__1, &cu, &
	    su);
    if (*m == 1) {
	return 0;
    }

/*     Compute V'*F1. */

    drot_(&c__2, &f[f_dim1 + 2], m, &f[f_dim1 + 1], m, &cv, &sv);

/*               ( F1 ) */
/*     Form  F = (    ) . */
/*               ( 0  ) */

    if (*m > *n) {
	i__1 = *m - *n;
	dlaset_("Full", &i__1, n, &c_b4, &c_b4, &f[*n + 1 + f_dim1], m, (
		ftnlen)4);
    }

/*     Compute H1*H2*F. */

    if (*m > 2) {
	i__1 = *m - 1;
	dlatzm_("Left", &i__1, n, &b[b_dim1 * 3 + 2], n, &tau2, &f[f_dim1 + 2]
		, &f[f_dim1 + 3], m, &dwork[1], (ftnlen)4);
    }
    dlatzm_("Left", m, n, &b[(b_dim1 << 1) + 1], n, &tau1, &f[f_dim1 + 1], &f[
	    f_dim1 + 2], m, &dwork[1], (ftnlen)4);

    return 0;
/* *** Last line of SB01BY *** */
} /* sb01by_ */

