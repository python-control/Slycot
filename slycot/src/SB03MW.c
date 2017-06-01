/* SB03MW.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb03mw_(logical *ltran, logical *lupper, doublereal *t, 
	integer *ldt, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, t_dim1, t_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t9[9]	/* was [3][3] */;
    static integer ip, jp;
    static doublereal eps, tmp[3], btmp[3], temp, smin;
    static integer jpiv[3];
    static doublereal xmax;
    static integer ipsv, jpsv;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
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

/*     To solve for the 2-by-2 symmetric matrix X in */

/*            op(T)'*X + X*op(T) = SCALE*B, */

/*     where T is 2-by-2, B is symmetric 2-by-2, and op(T) = T or T', */
/*     where T' denotes the transpose of T. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRAN   LOGICAL */
/*             Specifies the form of op(T) to be used, as follows: */
/*             = .FALSE.:  op(T) = T, */
/*             = .TRUE. :  op(T) = T'. */

/*     LUPPER  LOGICAL */
/*             Specifies which triangle of the matrix B is used, and */
/*             which triangle of the matrix X is computed, as follows: */
/*             = .TRUE. :  The upper triangular part; */
/*             = .FALSE.:  The lower triangular part. */

/*     Input/Output Parameters */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,2) */
/*             The leading 2-by-2 part of this array must contain the */
/*             matrix T. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= 2. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
/*             On entry with LUPPER = .TRUE., the leading 2-by-2 upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix B and the strictly */
/*             lower triangular part of B is not referenced. */
/*             On entry with LUPPER = .FALSE., the leading 2-by-2 lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix B and the strictly */
/*             upper triangular part of B is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor. SCALE is chosen less than or equal to 1 */
/*             to prevent the solution overflowing. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,2) */
/*             On exit with LUPPER = .TRUE., the leading 2-by-2 upper */
/*             triangular part of this array contains the upper */
/*             triangular part of the symmetric solution matrix X and the */
/*             strictly lower triangular part of X is not referenced. */
/*             On exit with LUPPER = .FALSE., the leading 2-by-2 lower */
/*             triangular part of this array contains the lower */
/*             triangular part of the symmetric solution matrix X and the */
/*             strictly upper triangular part of X is not referenced. */
/*             Note that X may be identified with B in the calling */
/*             statement. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 2. */

/*     XNORM   (output) DOUBLE PRECISION */
/*             The infinity-norm of the solution. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if T and -T have too close eigenvalues, so T */
/*                   is perturbed to get a nonsingular equation. */

/*             NOTE: In the interests of speed, this routine does not */
/*                   check the inputs for errors. */

/*     METHOD */

/*     The equivalent linear algebraic system of equations is formed and */
/*     solved using Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is stable and reliable, since Gaussian elimination */
/*     with complete pivoting is used. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Based on DLALY2 by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Continuous-time system, Lyapunov equation, matrix algebra. */

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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    *info = 0;

/*     Set constants to control overflow */

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1) / eps;

/*     Solve equivalent 3-by-3 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

/* Computing MAX */
/* Computing MAX */
    d__6 = (d__1 = t[t_dim1 + 1], abs(d__1)), d__7 = (d__2 = t[(t_dim1 << 1) 
	    + 1], abs(d__2)), d__6 = max(d__6,d__7), d__7 = (d__3 = t[t_dim1 
	    + 2], abs(d__3)), d__6 = max(d__6,d__7), d__7 = (d__4 = t[(t_dim1 
	    << 1) + 2], abs(d__4));
    d__5 = max(d__6,d__7) * eps;
    smin = max(d__5,smlnum);
    t9[6] = 0.;
    t9[2] = 0.;
    t9[0] = t[t_dim1 + 1];
    t9[4] = t[t_dim1 + 1] + t[(t_dim1 << 1) + 2];
    t9[8] = t[(t_dim1 << 1) + 2];
    if (*ltran) {
	t9[3] = t[(t_dim1 << 1) + 1];
	t9[1] = t[t_dim1 + 2];
	t9[7] = t[(t_dim1 << 1) + 1];
	t9[5] = t[t_dim1 + 2];
    } else {
	t9[3] = t[t_dim1 + 2];
	t9[1] = t[(t_dim1 << 1) + 1];
	t9[7] = t[t_dim1 + 2];
	t9[5] = t[(t_dim1 << 1) + 1];
    }
    btmp[0] = b[b_dim1 + 1] / 2.;
    if (*lupper) {
	btmp[1] = b[(b_dim1 << 1) + 1];
    } else {
	btmp[1] = b[b_dim1 + 2];
    }
    btmp[2] = b[(b_dim1 << 1) + 2] / 2.;

/*     Perform elimination */

    for (i__ = 1; i__ <= 2; ++i__) {
	xmax = 0.;

	for (ip = i__; ip <= 3; ++ip) {

	    for (jp = i__; jp <= 3; ++jp) {
		if ((d__1 = t9[ip + jp * 3 - 4], abs(d__1)) >= xmax) {
		    xmax = (d__1 = t9[ip + jp * 3 - 4], abs(d__1));
		    ipsv = ip;
		    jpsv = jp;
		}
/* L10: */
	    }

/* L20: */
	}

	if (ipsv != i__) {
	    dswap_(&c__3, &t9[ipsv - 1], &c__3, &t9[i__ - 1], &c__3);
	    temp = btmp[i__ - 1];
	    btmp[i__ - 1] = btmp[ipsv - 1];
	    btmp[ipsv - 1] = temp;
	}
	if (jpsv != i__) {
	    dswap_(&c__3, &t9[jpsv * 3 - 3], &c__1, &t9[i__ * 3 - 3], &c__1);
	}
	jpiv[i__ - 1] = jpsv;
	if ((d__1 = t9[i__ + i__ * 3 - 4], abs(d__1)) < smin) {
	    *info = 1;
	    t9[i__ + i__ * 3 - 4] = smin;
	}

	for (j = i__ + 1; j <= 3; ++j) {
	    t9[j + i__ * 3 - 4] /= t9[i__ + i__ * 3 - 4];
	    btmp[j - 1] -= t9[j + i__ * 3 - 4] * btmp[i__ - 1];

	    for (k = i__ + 1; k <= 3; ++k) {
		t9[j + k * 3 - 4] -= t9[j + i__ * 3 - 4] * t9[i__ + k * 3 - 4]
			;
/* L30: */
	    }

/* L40: */
	}

/* L50: */
    }

    if (abs(t9[8]) < smin) {
	t9[8] = smin;
    }
    *scale = 1.;
    if (smlnum * 4. * abs(btmp[0]) > abs(t9[0]) || smlnum * 4. * abs(btmp[1]) 
	    > abs(t9[4]) || smlnum * 4. * abs(btmp[2]) > abs(t9[8])) {
/* Computing MAX */
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]);
	*scale = .25 / max(d__1,d__2);
	btmp[0] *= *scale;
	btmp[1] *= *scale;
	btmp[2] *= *scale;
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	k = 4 - i__;
	temp = 1. / t9[k + k * 3 - 4];
	tmp[k - 1] = btmp[k - 1] * temp;

	for (j = k + 1; j <= 3; ++j) {
	    tmp[k - 1] -= temp * t9[k + j * 3 - 4] * tmp[j - 1];
/* L60: */
	}

/* L70: */
    }

    for (i__ = 1; i__ <= 2; ++i__) {
	if (jpiv[3 - i__ - 1] != 3 - i__) {
	    temp = tmp[3 - i__ - 1];
	    tmp[3 - i__ - 1] = tmp[jpiv[3 - i__ - 1] - 1];
	    tmp[jpiv[3 - i__ - 1] - 1] = temp;
	}
/* L80: */
    }

    x[x_dim1 + 1] = tmp[0];
    if (*lupper) {
	x[(x_dim1 << 1) + 1] = tmp[1];
    } else {
	x[x_dim1 + 2] = tmp[1];
    }
    x[(x_dim1 << 1) + 2] = tmp[2];
/* Computing MAX */
    d__1 = abs(tmp[0]) + abs(tmp[1]), d__2 = abs(tmp[1]) + abs(tmp[2]);
    *xnorm = max(d__1,d__2);

    return 0;
/* *** Last line of SB03MW *** */
} /* sb03mw_ */

