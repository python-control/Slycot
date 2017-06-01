/* SB01FY.f -- translated by f2c (version 20100827).
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
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b15 = 0.;
static doublereal c_b16 = 1.;

/* Subroutine */ int sb01fy_(logical *discr, integer *n, integer *m, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	f, integer *ldf, doublereal *v, integer *ldv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, v_dim1, 
	    v_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal u[4]	/* was [2][2] */, r11, r12, cs, r22, at[4]	
	    /* was [2][2] */, sn, temp;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal scale;
    extern /* Subroutine */ int mb04ox_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), sb03oy_(logical *, logical *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal dummy[4]	/* was [2][2] */;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlatzm_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dtrtri_(char *, char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the inner denominator of a right-coprime factorization */
/*     of a system of order N, where N is either 1 or 2. Specifically, */
/*     given the N-by-N unstable system state matrix A and the N-by-M */
/*     system input matrix B, an M-by-N state-feedback matrix F and */
/*     an M-by-M matrix V are constructed, such that the system */
/*     (A + B*F, B*V, F, V) is inner. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the type of system as follows: */
/*             = .FALSE.:  continuous-time system; */
/*             = .TRUE. :  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and also the number of rows of */
/*             the matrix B and the number of columns of the matrix F. */
/*             N is either 1 or 2. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices B and V, and also */
/*             the number of rows of the matrix F.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A whose eigenvalues must have positive */
/*             real parts if DISCR = .FALSE. or moduli greater than unity */
/*             if DISCR = .TRUE.. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the state- */
/*             feedback matrix F which assigns one eigenvalue (if N = 1) */
/*             or two eigenvalues (if N = 2) of the matrix A + B*F in */
/*             symmetric positions with respect to the imaginary axis */
/*             (if DISCR = .FALSE.) or the unit circle (if */
/*             DISCR = .TRUE.). */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,M) */
/*             The leading M-by-M upper triangular part of this array */
/*             contains the input/output matrix V of the resulting inner */
/*             system in upper triangular form. */
/*             If DISCR = .FALSE., the resulting V is an identity matrix. */

/*     LDV     INTEGER */
/*             The leading dimension of array V.  LDF >= MAX(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if uncontrollability of the pair (A,B) is detected; */
/*             = 2:  if A is stable or at the stability limit; */
/*             = 3:  if N = 2 and A has a pair of real eigenvalues. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine RCFID2. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Feb. 1999, A. Varga, DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Coprime factorization, eigenvalue, eigenvalue assignment, */
/*     feedback control, pole placement, state-space model. */

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

/*     For efficiency reasons, the parameters are not checked. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;

    /* Function Body */
    *info = 0;

/*     Compute an N-by-N upper triangular R such that R'*R = B*B' and */
/*     find an upper triangular matrix U in the equation */

/*     A'*U'*U + U'*U*A = R'*R if DISCR = .FALSE. or */
/*     A'*U'*U*A - U'*U = R'*R if DISCR = .TRUE. . */

    ma02ad_("Full", n, m, &b[b_offset], ldb, &f[f_offset], ldf, (ftnlen)4);

    if (*n == 1) {

/*        The N = 1 case. */

	if (*m > 1) {
	    dlarfg_(m, &f[f_dim1 + 1], &f[f_dim1 + 2], &c__1, &temp);
	}
	r11 = (d__1 = f[f_dim1 + 1], abs(d__1));

/*        Make sure A is unstable or divergent and find U. */

	if (*discr) {
	    temp = (d__1 = a[a_dim1 + 1], abs(d__1));
	    if (temp <= 1.) {
		*info = 2;
		return 0;
	    } else {
		temp = r11 / sqrt((temp - 1.) * (temp + 1.));
	    }
	} else {
	    if (a[a_dim1 + 1] <= 0.) {
		*info = 2;
		return 0;
	    } else {
		temp = r11 / sqrt((d__1 = a[a_dim1 + 1] * 2., abs(d__1)));
	    }
	}
	u[0] = temp;
	scale = 1.;
    } else {

/*        The N = 2 case. */

	if (*m > 1) {
	    dlarfg_(m, &f[f_dim1 + 1], &f[f_dim1 + 2], &c__1, &temp);
	    i__1 = *n - 1;
	    dlatzm_("Left", m, &i__1, &f[f_dim1 + 2], &c__1, &temp, &f[(
		    f_dim1 << 1) + 1], &f[(f_dim1 << 1) + 2], ldf, &v[
		    v_offset], (ftnlen)4);
	}
	r11 = f[f_dim1 + 1];
	r12 = f[(f_dim1 << 1) + 1];
	if (*m > 2) {
	    i__1 = *m - 1;
	    dlarfg_(&i__1, &f[(f_dim1 << 1) + 2], &f[(f_dim1 << 1) + 3], &
		    c__1, &temp);
	}
	if (*m == 1) {
	    r22 = 0.;
	} else {
	    r22 = f[(f_dim1 << 1) + 2];
	}
	at[0] = a[a_dim1 + 1];
	at[2] = a[a_dim1 + 2];
	at[1] = a[(a_dim1 << 1) + 1];
	at[3] = a[(a_dim1 << 1) + 2];
	u[0] = r11;
	u[2] = r12;
	u[3] = r22;
	sb03oy_(discr, &c_false, &c_n1, at, &c__2, u, &c__2, dummy, &c__2, &
		scale, info);
	if (*info != 0) {
	    if (*info != 4) {
		*info = 2;
	    } else {
		*info = 3;
	    }
	    return 0;
	}
    }

/*     Check the controllability of the pair (A,B). */

/*     Warning. Only an exact controllability check is performed. */
/*              If the pair (A,B) is nearly uncontrollable, then */
/*              the computed results may be inaccurate. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (u[i__ + (i__ << 1) - 3] == 0.) {
	    *info = 1;
	    return 0;
	}
/* L10: */
    }

/*     Set V = I. */

    dlaset_("Upper", m, m, &c_b15, &c_b16, &v[v_offset], ldv, (ftnlen)5);

    if (*discr) {

/*        Compute an upper triangular matrix V such that */
/*                                 -1 */
/*        V*V' = (I+B'*inv(U'*U)*B)  . */

/*        First compute F = B'*inv(U) and the Cholesky factorization */
/*        of I + F*F'. */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    f[i__ + f_dim1] = b[i__ * b_dim1 + 1] / u[0] * scale;
/* L20: */
	}
	if (*n == 2) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + (f_dim1 << 1)] = (b[i__ * b_dim1 + 2] - f[i__ + 
			f_dim1] * u[2]) / u[3] * scale;
/* L30: */
	    }
	    mb04ox_(m, &v[v_offset], ldv, &f[(f_dim1 << 1) + 1], &c__1);
	}
	mb04ox_(m, &v[v_offset], ldv, &f[f_dim1 + 1], &c__1);
	dtrtri_("Upper", "NonUnit", m, &v[v_offset], ldv, info, (ftnlen)5, (
		ftnlen)7);
    }

/*     Compute the feedback matrix F as: */

/*     1)   If DISCR = .FALSE. */

/*             F = -B'*inv(U'*U); */

/*     2)   If DISCR = .TRUE. */
/*                                -1 */
/*             F = -B'*(U'*U+B*B')  *A. */

    if (*n == 1) {
	if (*discr) {
	    temp = -a[a_dim1 + 1];
	    r11 = dlapy2_(u, &r11);
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + f_dim1] = b[i__ * b_dim1 + 1] / r11 / r11 * temp;
/* L40: */
	    }
	} else {
	    r11 = u[0];
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		f[i__ + f_dim1] = -(b[i__ * b_dim1 + 1] / r11 / r11);
/* L50: */
	    }
	}
    } else {

/*        Set R = U  if DISCR = .FALSE. or compute the Cholesky */
/*        factorization of R'*R = U'*U+B*B' if DISCR = .TRUE.. */

	if (*discr) {
	    temp = u[0];
	    drotg_(&r11, &temp, &cs, &sn);
	    temp = -sn * r12 + cs * u[2];
	    r12 = cs * r12 + sn * u[2];
	    r22 = dlapy3_(&r22, &temp, &u[3]);
	} else {
	    r11 = u[0];
	    r12 = u[2];
	    r22 = u[3];
	}

/*        Compute F = -B'*inv(R'*R). */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    f[i__ + f_dim1] = -b[i__ * b_dim1 + 1] / r11;
	    f[i__ + (f_dim1 << 1)] = -(b[i__ * b_dim1 + 2] + f[i__ + f_dim1] *
		     r12) / r22;
	    f[i__ + (f_dim1 << 1)] /= r22;
	    f[i__ + f_dim1] = (f[i__ + f_dim1] - f[i__ + (f_dim1 << 1)] * r12)
		     / r11;
/* L60: */
	}
	if (*discr) {

/*           Compute F <-- F*A. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp = f[i__ + f_dim1] * a[a_dim1 + 1] + f[i__ + (f_dim1 << 1)
			] * a[a_dim1 + 2];
		f[i__ + (f_dim1 << 1)] = f[i__ + f_dim1] * a[(a_dim1 << 1) + 
			1] + f[i__ + (f_dim1 << 1)] * a[(a_dim1 << 1) + 2];
		f[i__ + f_dim1] = temp;
/* L70: */
	    }
	}
    }

    return 0;
/* *** Last line of SB01FY *** */
} /* sb01fy_ */

