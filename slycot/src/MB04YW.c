/* MB04YW.f -- translated by f2c (version 20100827).
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

static doublereal c_b16 = 1.;

/* Subroutine */ int mb04yw_(logical *qrit, logical *updatu, logical *updatv, 
	integer *m, integer *n, integer *l, integer *k, doublereal *shift, 
	doublereal *d__, doublereal *e, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *dwork)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__;
    static doublereal r__, cs, sn;
    static integer nm1, nm12, nm13, ncv;
    static doublereal cosl, sinl, cosr, sinr;
    static integer irot;
    static doublereal oldcs;
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal oldsn;
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

/*     To perform either one QR or QL iteration step onto the unreduced */
/*     bidiagonal submatrix Jk: */

/*              |D(l) E(l)    0  ...    0   | */
/*              | 0   D(l+1) E(l+1)     .   | */
/*         Jk = | .                     .   | */
/*              | .                     .   | */
/*              | .                   E(k-1)| */
/*              | 0   ...        ...   D(k) | */

/*     with k <= p and l >= 1, p = MIN(M,N), of the bidiagonal matrix J: */

/*              |D(1) E(1)  0    ...   0   | */
/*              | 0   D(2) E(2)        .   | */
/*          J = | .                    .   |. */
/*              | .                    .   | */
/*              | .                  E(p-1)| */
/*              | 0   ...        ...  D(p) | */

/*     Hereby, Jk is transformed to  S' Jk T with S and T products of */
/*     Givens rotations. These Givens rotations S (respectively, T) are */
/*     postmultiplied into U (respectively, V), if UPDATU (respectively, */
/*     UPDATV) is .TRUE.. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     QRIT    LOGICAL */
/*             Indicates whether a QR or QL iteration step is to be */
/*             taken (from larger end diagonal element towards smaller), */
/*             as follows: */
/*             = .TRUE. :  QR iteration step (chase bulge from top to */
/*                         bottom); */
/*             = .FALSE.:  QL iteration step (chase bulge from bottom to */
/*                         top). */

/*     UPDATU  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the left-hand Givens rotations S, as follows: */
/*             = .FALSE.:  Do not form U; */
/*             = .TRUE. :  The given matrix U is updated (postmultiplied) */
/*                         by the left-hand Givens rotations S. */

/*     UPDATV  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the right-hand Givens rotations S, as follows: */
/*             = .FALSE.:  Do not form V; */
/*             = .TRUE. :  The given matrix V is updated (postmultiplied) */
/*                         by the right-hand Givens rotations T. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix U.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix V.  N >= 0. */

/*     L       (input) INTEGER */
/*             The index of the first diagonal entry of the considered */
/*             unreduced bidiagonal submatrix Jk of J. */

/*     K       (input) INTEGER */
/*             The index of the last diagonal entry of the considered */
/*             unreduced bidiagonal submatrix Jk of J. */

/*     SHIFT   (input) DOUBLE PRECISION */
/*             Value of the shift used in the QR or QL iteration step. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (p) */
/*             where p = MIN(M,N) */
/*             On entry, D must contain the diagonal entries of the */
/*             bidiagonal matrix J. */
/*             On exit, D contains the diagonal entries of the */
/*             transformed bidiagonal matrix S' J T. */

/*     E       (input/output) DOUBLE PRECISION array, dimension (p-1) */
/*             On entry, E must contain the superdiagonal entries of J. */
/*             On exit, E contains the superdiagonal entries of the */
/*             transformed matrix S' J T. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,p) */
/*             On entry, if UPDATU = .TRUE., U must contain the M-by-p */
/*             left transformation matrix. */
/*             On exit, if UPDATU = .TRUE., the Givens rotations S on the */
/*             left have been postmultiplied into U, i.e., U * S is */
/*             returned. */
/*             U is not referenced if UPDATU = .FALSE.. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= max(1,M) if UPDATU = .TRUE.; */
/*             LDU >= 1        if UPDATU = .FALSE.. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,p) */
/*             On entry, if UPDATV = .TRUE., V must contain the N-by-p */
/*             right transformation matrix. */
/*             On exit, if UPDATV = .TRUE., the Givens rotations T on the */
/*             right have been postmultiplied into V, i.e., V * T is */
/*             returned. */
/*             V is not referenced if UPDATV = .FALSE.. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= max(1,N) if UPDATV = .TRUE.; */
/*             LDV >= 1        if UPDATV = .FALSE.. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(1,LDWORK)) */
/*             LDWORK >= 4*MIN(M,N)-4, if UPDATU = UPDATV = .TRUE.; */
/*             LDWORK >= 2*MIN(M,N)-2, if */
/*                             UPDATU = .TRUE. and UPDATV = .FALSE. or */
/*                             UPDATV = .TRUE. and UPDATU = .FALSE.; */
/*             LDWORK >= 1, if UPDATU = UPDATV = .FALSE.. */

/*     METHOD */

/*     QR iterations diagonalize the bidiagonal matrix by zeroing the */
/*     super-diagonal elements of Jk from bottom to top. */
/*     QL iterations diagonalize the bidiagonal matrix by zeroing the */
/*     super-diagonal elements of Jk from top to bottom. */
/*     The routine overwrites Jk with the bidiagonal matrix S' Jk T, */
/*     where S and T are products of Givens rotations. */
/*     T is essentially the orthogonal matrix that would be obtained by */
/*     applying one implicit symmetric shift QR (QL) step onto the matrix */
/*     Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a */
/*     product of an orthogonal matrix T and a upper (lower) triangular */
/*     matrix. See [1,Sec.8.2-8.3] and [2] for more details. */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         Matrix Computations. */
/*         The Johns Hopkins University Press, Baltimore, Maryland, 1983. */

/*     [2] Bowdler, H., Martin, R.S. and Wilkinson, J.H. */
/*         The QR and QL algorithms for symmetric matrices. */
/*         Numer. Math., 11, pp. 293-306, 1968. */

/*     [3] Demmel, J. and Kahan, W. */
/*         Computing small singular values of bidiagonal matrices with */
/*         guaranteed high relative accuracy. */
/*         SIAM J. Sci. Statist. Comput., 11, pp. 873-912, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routines MB04QY and MB04QZ by S. Van */
/*     Huffel, Katholieke University Leuven, Belgium. */
/*     This subroutine is based on the QR/QL step implemented in LAPACK */
/*     routine DBDSQR. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bidiagonal matrix, orthogonal transformation, singular values. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     For speed, no tests of the input scalar arguments are done. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --dwork;

    /* Function Body */
    ncv = min(*m,*n);
    if (ncv <= 1 || *l == *k) {
	return 0;
    }

    nm1 = ncv - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    if (! (*updatv)) {
	nm12 = 0;
	nm13 = nm1;
    }

/*     If SHIFT = 0, do simplified QR iteration. */

    if (*shift == 0.) {
	if (*qrit) {

/*           Chase bulge from top to bottom. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

	    cs = 1.;
	    oldcs = 1.;
	    d__1 = d__[*l] * cs;
	    dlartg_(&d__1, &e[*l], &cs, &sn, &r__);
	    d__1 = oldcs * r__;
	    d__2 = d__[*l + 1] * sn;
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[*l]);
	    if (*updatv) {
		dwork[1] = cs;
		dwork[nm1 + 1] = sn;
	    }
	    if (*updatu) {
		dwork[nm12 + 1] = oldcs;
		dwork[nm13 + 1] = oldsn;
	    }
	    irot = 1;

	    i__1 = *k - 1;
	    for (i__ = *l + 1; i__ <= i__1; ++i__) {
		d__1 = d__[i__] * cs;
		dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
		e[i__ - 1] = oldsn * r__;
		d__1 = oldcs * r__;
		d__2 = d__[i__ + 1] * sn;
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
		++irot;
		if (*updatv) {
		    dwork[irot] = cs;
		    dwork[irot + nm1] = sn;
		}
		if (*updatu) {
		    dwork[irot + nm12] = oldcs;
		    dwork[irot + nm13] = oldsn;
		}
/* L110: */
	    }

	    h__ = d__[*k] * cs;
	    d__[*k] = h__ * oldcs;
	    e[*k - 1] = h__ * oldsn;

/*           Update U and/or V. */

	    if (*updatv) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "F", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	    }
	    if (*updatu) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "F", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
	    }

	} else {

/*           Chase bulge from bottom to top. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

	    cs = 1.;
	    oldcs = 1.;
	    d__1 = d__[*k] * cs;
	    dlartg_(&d__1, &e[*k - 1], &cs, &sn, &r__);
	    d__1 = oldcs * r__;
	    d__2 = d__[*k - 1] * sn;
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[*k]);
	    if (*updatv) {
		dwork[*k - *l] = oldcs;
		dwork[*k - *l + nm1] = -oldsn;
	    }
	    if (*updatu) {
		dwork[*k - *l + nm12] = cs;
		dwork[*k - *l + nm13] = -sn;
	    }
	    irot = *k - *l;

	    i__1 = *l + 1;
	    for (i__ = *k - 1; i__ >= i__1; --i__) {
		d__1 = d__[i__] * cs;
		dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
		e[i__] = oldsn * r__;
		d__1 = oldcs * r__;
		d__2 = d__[i__ - 1] * sn;
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
		--irot;
		if (*updatv) {
		    dwork[irot] = oldcs;
		    dwork[irot + nm1] = -oldsn;
		}
		if (*updatu) {
		    dwork[irot + nm12] = cs;
		    dwork[irot + nm13] = -sn;
		}
/* L120: */
	    }

	    h__ = d__[*l] * cs;
	    d__[*l] = h__ * oldcs;
	    e[*l] = h__ * oldsn;

/*           Update U and/or V. */

	    if (*updatv) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "B", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	    }
	    if (*updatu) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "B", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
	    }
	}
    } else {

/*        Use nonzero shift. */

	if (*qrit) {

/*           Chase bulge from top to bottom. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

	    f = ((d__1 = d__[*l], abs(d__1)) - *shift) * (d_sign(&c_b16, &d__[
		    *l]) + *shift / d__[*l]);
	    g = e[*l];
	    dlartg_(&f, &g, &cosr, &sinr, &r__);
	    f = cosr * d__[*l] + sinr * e[*l];
	    e[*l] = cosr * e[*l] - sinr * d__[*l];
	    g = sinr * d__[*l + 1];
	    d__[*l + 1] = cosr * d__[*l + 1];
	    dlartg_(&f, &g, &cosl, &sinl, &r__);
	    d__[*l] = r__;
	    f = cosl * e[*l] + sinl * d__[*l + 1];
	    d__[*l + 1] = cosl * d__[*l + 1] - sinl * e[*l];
	    g = sinl * e[*l + 1];
	    e[*l + 1] = cosl * e[*l + 1];
	    if (*updatv) {
		dwork[1] = cosr;
		dwork[nm1 + 1] = sinr;
	    }
	    if (*updatu) {
		dwork[nm12 + 1] = cosl;
		dwork[nm13 + 1] = sinl;
	    }
	    irot = 1;

	    i__1 = *k - 2;
	    for (i__ = *l + 1; i__ <= i__1; ++i__) {
		dlartg_(&f, &g, &cosr, &sinr, &r__);
		e[i__ - 1] = r__;
		f = cosr * d__[i__] + sinr * e[i__];
		e[i__] = cosr * e[i__] - sinr * d__[i__];
		g = sinr * d__[i__ + 1];
		d__[i__ + 1] = cosr * d__[i__ + 1];
		dlartg_(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__] + sinl * d__[i__ + 1];
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
		g = sinl * e[i__ + 1];
		e[i__ + 1] = cosl * e[i__ + 1];
		++irot;
		if (*updatv) {
		    dwork[irot] = cosr;
		    dwork[irot + nm1] = sinr;
		}
		if (*updatu) {
		    dwork[irot + nm12] = cosl;
		    dwork[irot + nm13] = sinl;
		}
/* L130: */
	    }

	    if (*l < *k - 1) {
		dlartg_(&f, &g, &cosr, &sinr, &r__);
		e[*k - 2] = r__;
		f = cosr * d__[*k - 1] + sinr * e[*k - 1];
		e[*k - 1] = cosr * e[*k - 1] - sinr * d__[*k - 1];
		g = sinr * d__[*k];
		d__[*k] = cosr * d__[*k];
		dlartg_(&f, &g, &cosl, &sinl, &r__);
		d__[*k - 1] = r__;
		f = cosl * e[*k - 1] + sinl * d__[*k];
		d__[*k] = cosl * d__[*k] - sinl * e[*k - 1];
		++irot;
		if (*updatv) {
		    dwork[irot] = cosr;
		    dwork[irot + nm1] = sinr;
		}
		if (*updatu) {
		    dwork[irot + nm12] = cosl;
		    dwork[irot + nm13] = sinl;
		}
	    }
	    e[*k - 1] = f;

/*           Update U and/or V. */

	    if (*updatv) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "F", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	    }
	    if (*updatu) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "F", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
	    }

	} else {

/*           Chase bulge from bottom to top. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

	    f = ((d__1 = d__[*k], abs(d__1)) - *shift) * (d_sign(&c_b16, &d__[
		    *k]) + *shift / d__[*k]);
	    g = e[*k - 1];
	    if (*l < *k - 1) {
		dlartg_(&f, &g, &cosr, &sinr, &r__);
		f = cosr * d__[*k] + sinr * e[*k - 1];
		e[*k - 1] = cosr * e[*k - 1] - sinr * d__[*k];
		g = sinr * d__[*k - 1];
		d__[*k - 1] = cosr * d__[*k - 1];
		dlartg_(&f, &g, &cosl, &sinl, &r__);
		d__[*k] = r__;
		f = cosl * e[*k - 1] + sinl * d__[*k - 1];
		d__[*k - 1] = cosl * d__[*k - 1] - sinl * e[*k - 1];
		g = sinl * e[*k - 2];
		e[*k - 2] = cosl * e[*k - 2];
		if (*updatv) {
		    dwork[*k - *l] = cosl;
		    dwork[*k - *l + nm1] = -sinl;
		}
		if (*updatu) {
		    dwork[*k - *l + nm12] = cosr;
		    dwork[*k - *l + nm13] = -sinr;
		}
		irot = *k - *l;
	    } else {
		irot = *k - *l + 1;
	    }

	    i__1 = *l + 2;
	    for (i__ = *k - 1; i__ >= i__1; --i__) {
		dlartg_(&f, &g, &cosr, &sinr, &r__);
		e[i__] = r__;
		f = cosr * d__[i__] + sinr * e[i__ - 1];
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
		g = sinr * d__[i__ - 1];
		d__[i__ - 1] = cosr * d__[i__ - 1];
		dlartg_(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
		g = sinl * e[i__ - 2];
		e[i__ - 2] = cosl * e[i__ - 2];
		--irot;
		if (*updatv) {
		    dwork[irot] = cosl;
		    dwork[irot + nm1] = -sinl;
		}
		if (*updatu) {
		    dwork[irot + nm12] = cosr;
		    dwork[irot + nm13] = -sinr;
		}
/* L140: */
	    }

	    dlartg_(&f, &g, &cosr, &sinr, &r__);
	    e[*l + 1] = r__;
	    f = cosr * d__[*l + 1] + sinr * e[*l];
	    e[*l] = cosr * e[*l] - sinr * d__[*l + 1];
	    g = sinr * d__[*l];
	    d__[*l] = cosr * d__[*l];
	    dlartg_(&f, &g, &cosl, &sinl, &r__);
	    d__[*l + 1] = r__;
	    f = cosl * e[*l] + sinl * d__[*l];
	    d__[*l] = cosl * d__[*l] - sinl * e[*l];
	    --irot;
	    if (*updatv) {
		dwork[irot] = cosl;
		dwork[irot + nm1] = -sinl;
	    }
	    if (*updatu) {
		dwork[irot + nm12] = cosr;
		dwork[irot + nm13] = -sinr;
	    }
	    e[*l] = f;

/*           Update U and/or V if desired. */

	    if (*updatv) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "B", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	    }
	    if (*updatu) {
		i__1 = *k - *l + 1;
		dlasr_("R", "V", "B", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
	    }
	}
    }

    return 0;
/* *** Last line of MB04YW *** */
} /* mb04yw_ */

