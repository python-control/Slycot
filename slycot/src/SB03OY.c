/* SB03OY.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 1.;

/* Subroutine */ int sb03oy_(logical *discr, logical *ltrans, integer *isgn, 
	doublereal *s, integer *lds, doublereal *r__, integer *ldr, 
	doublereal *a, integer *lda, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, s_dim1, s_offset;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal g[2], t[2], y[2], e1, e2, p1, p2[2], p3, v1, v2[2], v3, 
	    dp[2], s11, s12, s21, s22, dt[2], x11[2], x12[2], x21[2], x22[2], 
	    p3i, p3r, eta, csp[2], csq[2], eps, sgn, cst[2], snp, snq, snt, 
	    absb, absg, abst, temp[2], smin, gamma[2], alpha, delta[2];
    extern /* Subroutine */ int sb03ov_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tempi, tempr;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scaloc, bignum, smlnum;


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

/*     To solve for the Cholesky factor  U  of  X, */

/*        op(U)'*op(U) = X, */

/*     where  U  is a two-by-two upper triangular matrix, either the */
/*     continuous-time two-by-two Lyapunov equation */
/*                                         2 */
/*         op(S)'*X + X*op(S) = -ISGN*scale *op(R)'*op(R), */

/*     when DISCR = .FALSE., or the discrete-time two-by-two Lyapunov */
/*     equation */
/*                                         2 */
/*         op(S)'*X*op(S) - X = -ISGN*scale *op(R)'*op(R), */

/*     when DISCR = .TRUE., where op(K) = K or K' (i.e., the transpose of */
/*     the matrix K),  S  is a two-by-two matrix with complex conjugate */
/*     eigenvalues,  R  is a two-by-two upper triangular matrix, */
/*     ISGN = -1 or 1,  and  scale  is an output scale factor, set less */
/*     than or equal to 1 to avoid overflow in  X.  The routine also */
/*     computes two matrices, B and A, so that */
/*                                   2 */
/*        B*U = U*S  and  A*U = scale *R,  if  LTRANS = .FALSE.,  or */
/*                                   2 */
/*        U*B = S*U  and  U*A = scale *R,  if  LTRANS = .TRUE., */
/*     which are used by the general Lyapunov solver. */
/*     In the continuous-time case  ISGN*S  must be stable, so that its */
/*     eigenvalues must have strictly negative real parts. */
/*     In the discrete-time case  S  must be convergent if ISGN = 1, that */
/*     is, its eigenvalues must have moduli less than unity, or  S  must */
/*     be completely divergent if ISGN = -1, that is, its eigenvalues */
/*     must have moduli greater than unity. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the equation to be solved:       2 */
/*             = .FALSE.: op(S)'*X + X*op(S) = -ISGN*scale *op(R)'*op(R); */
/*                                                        2 */
/*             = .TRUE. : op(S)'*X*op(S) - X = -ISGN*scale *op(R)'*op(R). */

/*     LTRANS  LOGICAL */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = .FALSE.:  op(K) = K    (No transpose); */
/*             = .TRUE. :  op(K) = K**T (Transpose). */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     S       (input/output) DOUBLE PRECISION array, dimension (LDS,2) */
/*             On entry, S must contain a 2-by-2 matrix. */
/*             On exit, S contains a 2-by-2 matrix B such that B*U = U*S, */
/*             if LTRANS = .FALSE., or U*B = S*U, if LTRANS = .TRUE.. */
/*             Notice that if U is nonsingular then */
/*               B = U*S*inv( U ),  if LTRANS = .FALSE. */
/*               B = inv( U )*S*U,  if LTRANS = .TRUE.. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= 2. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,2) */
/*             On entry, R must contain a 2-by-2 upper triangular matrix. */
/*             The element R( 2, 1 ) is not referenced. */
/*             On exit, R contains U, the 2-by-2 upper triangular */
/*             Cholesky factor of the solution X, X = op(U)'*op(U). */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= 2. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,2) */
/*             A contains a 2-by-2 upper triangular matrix A satisfying */
/*             A*U/scale = scale*R, if LTRANS = .FALSE., or */
/*             U*A/scale = scale*R, if LTRANS = .TRUE.. */
/*             Notice that if U is nonsingular then */
/*               A = scale*scale*R*inv( U ),  if LTRANS = .FALSE. */
/*               A = scale*scale*inv( U )*R,  if LTRANS = .TRUE.. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if the Lyapunov equation is (nearly) singular */
/*                   (warning indicator); */
/*                   if DISCR = .FALSE., this means that while the */
/*                   matrix S has computed eigenvalues with negative real */
/*                   parts, it is only just stable in the sense that */
/*                   small perturbations in S can make one or more of the */
/*                   eigenvalues have a non-negative real part; */
/*                   if DISCR = .TRUE., this means that while the */
/*                   matrix S has computed eigenvalues inside the unit */
/*                   circle, it is nevertheless only just convergent, in */
/*                   the sense that small perturbations in S can make one */
/*                   or more of the eigenvalues lie outside the unit */
/*                   circle; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrix S is unchanged); */
/*             = 2:  if DISCR = .FALSE., and ISGN*S is not stable or */
/*                   if DISCR = .TRUE., ISGN = 1 and S is not convergent */
/*                   or if DISCR = .TRUE., ISGN = -1 and S is not */
/*                   completely divergent; */
/*             = 4:  if S has real eigenvalues. */

/*     NOTE: In the interests of speed, this routine does not check all */
/*           inputs for errors. */

/*     METHOD */

/*     The LAPACK scheme for solving 2-by-2 Sylvester equations is */
/*     adapted for 2-by-2 Lyapunov equations, but directly computing the */
/*     Cholesky factor of the solution. */

/*     REFERENCES */

/*     [1] Hammarling S. J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB03CY by Sven Hammarling, */
/*     NAG Ltd., United Kingdom, November 1986. */
/*     Partly based on SB03CY and PLYAP2 by A. Varga, University of */
/*     Bochum, May 1992. */

/*     REVISIONS */

/*     Dec. 1997, April 1998. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ***************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     The comments in this routine refer to notation and equation */
/*     numbers in sections 6 and 10 of [1]. */

/*     Find the eigenvalue  lambda = E1 - i*E2  of s11. */

    /* Parameter adjustments */
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    sgn = (doublereal) (*isgn);
    s11 = s[s_dim1 + 1];
    s12 = s[(s_dim1 << 1) + 1];
    s21 = s[s_dim1 + 2];
    s22 = s[(s_dim1 << 1) + 2];

/*     Set constants to control overflow. */

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1);
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = smlnum * 4. / eps;
    bignum = 1. / smlnum;

/* Computing MAX */
/* Computing MAX */
    d__3 = abs(s11), d__4 = abs(s12), d__3 = max(d__3,d__4), d__4 = abs(s21), 
	    d__3 = max(d__3,d__4), d__4 = abs(s22);
    d__1 = smlnum, d__2 = eps * max(d__3,d__4);
    smin = max(d__1,d__2);
    *scale = 1.;

    dlanv2_(&s11, &s12, &s21, &s22, &tempr, &tempi, &e1, &e2, csp, csq);
    if (tempi == 0.) {
	*info = 4;
	return 0;
    }
    absb = dlapy2_(&e1, &e2);
    if (*discr) {
	if (sgn * (absb - 1.) >= 0.) {
	    *info = 2;
	    return 0;
	}
    } else {
	if (sgn * e1 >= 0.) {
	    *info = 2;
	    return 0;
	}
    }

/*     Compute the cos and sine that define  Qhat.  The sine is real. */

    temp[0] = s[s_dim1 + 1] - e1;
    temp[1] = e2;
    if (*ltrans) {
	temp[1] = -e2;
    }
    sb03ov_(temp, &s[s_dim1 + 2], csq, &snq);

/*     beta in (6.9) is given by  beta = E1 + i*E2,  compute  t. */

    temp[0] = csq[0] * s[(s_dim1 << 1) + 1] - snq * s[s_dim1 + 1];
    temp[1] = csq[1] * s[(s_dim1 << 1) + 1];
    tempr = csq[0] * s[(s_dim1 << 1) + 2] - snq * s[s_dim1 + 2];
    tempi = csq[1] * s[(s_dim1 << 1) + 2];
    t[0] = csq[0] * temp[0] - csq[1] * temp[1] + snq * tempr;
    t[1] = csq[0] * temp[1] + csq[1] * temp[0] + snq * tempi;

    if (*ltrans) {
/*                                                         (     -- ) */
/*        Case op(M) = M'.  Note that the modified  R  is  ( p3  p2 ). */
/*                                                         ( 0   p1 ) */

/*        Compute the cos and sine that define  Phat. */

	temp[0] = csq[0] * r__[(r_dim1 << 1) + 2] - snq * r__[(r_dim1 << 1) + 
		1];
	temp[1] = -csq[1] * r__[(r_dim1 << 1) + 2];
	d__1 = -snq * r__[r_dim1 + 1];
	sb03ov_(temp, &d__1, csp, &snp);

/*        Compute p1, p2 and p3 of the relation corresponding to (6.11). */

	p1 = temp[0];
	temp[0] = csq[0] * r__[(r_dim1 << 1) + 1] + snq * r__[(r_dim1 << 1) + 
		2];
	temp[1] = -csq[1] * r__[(r_dim1 << 1) + 1];
	tempr = csq[0] * r__[r_dim1 + 1];
	tempi = -csq[1] * r__[r_dim1 + 1];
	p2[0] = csp[0] * temp[0] - csp[1] * temp[1] + snp * tempr;
	p2[1] = -csp[0] * temp[1] - csp[1] * temp[0] - snp * tempi;
	p3r = csp[0] * tempr + csp[1] * tempi - snp * temp[0];
	p3i = csp[0] * tempi - csp[1] * tempr - snp * temp[1];
    } else {

/*        Case op(M) = M. */

/*        Compute the cos and sine that define  Phat. */

	temp[0] = csq[0] * r__[r_dim1 + 1] + snq * r__[(r_dim1 << 1) + 1];
	temp[1] = csq[1] * r__[r_dim1 + 1];
	d__1 = snq * r__[(r_dim1 << 1) + 2];
	sb03ov_(temp, &d__1, csp, &snp);

/*        Compute p1, p2 and p3 of (6.11). */

	p1 = temp[0];
	temp[0] = csq[0] * r__[(r_dim1 << 1) + 1] - snq * r__[r_dim1 + 1];
	temp[1] = csq[1] * r__[(r_dim1 << 1) + 1];
	tempr = csq[0] * r__[(r_dim1 << 1) + 2];
	tempi = csq[1] * r__[(r_dim1 << 1) + 2];
	p2[0] = csp[0] * temp[0] - csp[1] * temp[1] + snp * tempr;
	p2[1] = csp[0] * temp[1] + csp[1] * temp[0] + snp * tempi;
	p3r = csp[0] * tempr + csp[1] * tempi - snp * temp[0];
	p3i = csp[1] * tempr - csp[0] * tempi + snp * temp[1];
    }

/*     Make  p3  real by multiplying by  conjg ( p3 )/abs( p3 )  to give */

/*     p3 := abs( p3 ). */

    if (p3i == 0.) {
	p3 = abs(p3r);
	dp[0] = d_sign(&c_b4, &p3r);
	dp[1] = 0.;
    } else {
	p3 = dlapy2_(&p3r, &p3i);
	dp[0] = p3r / p3;
	dp[1] = -p3i / p3;
    }

/*     Now compute the quantities v1, v2, v3 and y in (6.13) - (6.15), */
/*     or (10.23) - (10.25). Care is taken to avoid overflows. */

    if (*discr) {
	alpha = sqrt((d__1 = 1. - absb, abs(d__1)) * (absb + 1.));
    } else {
	alpha = sqrt((d__1 = e1 * 2., abs(d__1)));
    }

    scaloc = 1.;
    if (alpha < smin) {
	alpha = smin;
	*info = 1;
    }
    abst = abs(p1);
    if (alpha < 1. && abst > 1.) {
	if (abst > bignum * alpha) {
	    scaloc = 1. / abst;
	}
    }
    if (scaloc != 1.) {
	p1 = scaloc * p1;
	p2[0] = scaloc * p2[0];
	p2[1] = scaloc * p2[1];
	p3 = scaloc * p3;
	*scale = scaloc * *scale;
    }
    v1 = p1 / alpha;

    if (*discr) {
/* Computing 2nd power */
	d__1 = e2;
	g[0] = (1. - e1) * (e1 + 1.) + d__1 * d__1;
	g[1] = e1 * -2. * e2;
	absg = dlapy2_(g, &g[1]);
	scaloc = 1.;
	if (absg < smin) {
	    absg = smin;
	    *info = 1;
	}
	temp[0] = sgn * alpha * p2[0] + v1 * (e1 * t[0] - e2 * t[1]);
	temp[1] = sgn * alpha * p2[1] + v1 * (e1 * t[1] + e2 * t[0]);
/* Computing MAX */
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
	abst = max(d__1,d__2);
	if (absg < 1. && abst > 1.) {
	    if (abst > bignum * absg) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    v1 = scaloc * v1;
	    temp[0] = scaloc * temp[0];
	    temp[1] = scaloc * temp[1];
	    p1 = scaloc * p1;
	    p2[0] = scaloc * p2[0];
	    p2[1] = scaloc * p2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	temp[0] /= absg;
	temp[1] /= absg;

	scaloc = 1.;
	v2[0] = g[0] * temp[0] + g[1] * temp[1];
	v2[1] = g[0] * temp[1] - g[1] * temp[0];
/* Computing MAX */
	d__1 = abs(v2[0]), d__2 = abs(v2[1]);
	abst = max(d__1,d__2);
	if (absg < 1. && abst > 1.) {
	    if (abst > bignum * absg) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    v1 = scaloc * v1;
	    v2[0] = scaloc * v2[0];
	    v2[1] = scaloc * v2[1];
	    p1 = scaloc * p1;
	    p2[0] = scaloc * p2[0];
	    p2[1] = scaloc * p2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	v2[0] /= absg;
	v2[1] /= absg;

	scaloc = 1.;
	temp[0] = p1 * t[0] - e2 * 2. * p2[1];
	temp[1] = p1 * t[1] + e2 * 2. * p2[0];
/* Computing MAX */
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
	abst = max(d__1,d__2);
	if (absg < 1. && abst > 1.) {
	    if (abst > bignum * absg) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    temp[0] = scaloc * temp[0];
	    temp[1] = scaloc * temp[1];
	    v1 = scaloc * v1;
	    v2[0] = scaloc * v2[0];
	    v2[1] = scaloc * v2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	temp[0] /= absg;
	temp[1] /= absg;

	scaloc = 1.;
	y[0] = -(g[0] * temp[0] + g[1] * temp[1]);
	y[1] = -(g[0] * temp[1] - g[1] * temp[0]);
/* Computing MAX */
	d__1 = abs(y[0]), d__2 = abs(y[1]);
	abst = max(d__1,d__2);
	if (absg < 1. && abst > 1.) {
	    if (abst > bignum * absg) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    y[0] = scaloc * y[0];
	    y[1] = scaloc * y[1];
	    v1 = scaloc * v1;
	    v2[0] = scaloc * v2[0];
	    v2[1] = scaloc * v2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	y[0] /= absg;
	y[1] /= absg;
    } else {

	scaloc = 1.;
	if (absb < smin) {
	    absb = smin;
	    *info = 1;
	}
	temp[0] = sgn * alpha * p2[0] + v1 * t[0];
	temp[1] = sgn * alpha * p2[1] + v1 * t[1];
/* Computing MAX */
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
	abst = max(d__1,d__2);
	if (absb < 1. && abst > 1.) {
	    if (abst > bignum * absb) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    v1 = scaloc * v1;
	    temp[0] = scaloc * temp[0];
	    temp[1] = scaloc * temp[1];
	    p2[0] = scaloc * p2[0];
	    p2[1] = scaloc * p2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	temp[0] /= absb * 2.;
	temp[1] /= absb * 2.;
	scaloc = 1.;
	v2[0] = -(e1 * temp[0] + e2 * temp[1]);
	v2[1] = -(e1 * temp[1] - e2 * temp[0]);
/* Computing MAX */
	d__1 = abs(v2[0]), d__2 = abs(v2[1]);
	abst = max(d__1,d__2);
	if (absb < 1. && abst > 1.) {
	    if (abst > bignum * absb) {
		scaloc = 1. / abst;
	    }
	}
	if (scaloc != 1.) {
	    v1 = scaloc * v1;
	    v2[0] = scaloc * v2[0];
	    v2[1] = scaloc * v2[1];
	    p2[0] = scaloc * p2[0];
	    p2[1] = scaloc * p2[1];
	    p3 = scaloc * p3;
	    *scale = scaloc * *scale;
	}
	v2[0] /= absb;
	v2[1] /= absb;
	y[0] = p2[0] - alpha * v2[0];
	y[1] = p2[1] - alpha * v2[1];
    }

    scaloc = 1.;
    v3 = dlapy3_(&p3, y, &y[1]);
    if (alpha < 1. && v3 > 1.) {
	if (v3 > bignum * alpha) {
	    scaloc = 1. / v3;
	}
    }
    if (scaloc != 1.) {
	v1 = scaloc * v1;
	v2[0] = scaloc * v2[0];
	v2[1] = scaloc * v2[1];
	v3 = scaloc * v3;
	p3 = scaloc * p3;
	*scale = scaloc * *scale;
    }
    v3 /= alpha;

    if (*ltrans) {

/*        Case op(M) = M'. */

/*        Form  X = conjg( Qhat' )*v11. */

	x11[0] = csq[0] * v3;
	x11[1] = csq[1] * v3;
	x21[0] = snq * v3;
	x12[0] = csq[0] * v2[0] + csq[1] * v2[1] - snq * v1;
	x12[1] = -csq[0] * v2[1] + csq[1] * v2[0];
	x22[0] = csq[0] * v1 + snq * v2[0];
	x22[1] = -csq[1] * v1 - snq * v2[1];

/*        Obtain u11 from the RQ-factorization of X. The conjugate of */
/*        X22 should be taken. */

	x22[1] = -x22[1];
	sb03ov_(x22, x21, cst, &snt);
	r__[(r_dim1 << 1) + 2] = x22[0];
	r__[(r_dim1 << 1) + 1] = cst[0] * x12[0] - cst[1] * x12[1] + snt * 
		x11[0];
	tempr = cst[0] * x11[0] + cst[1] * x11[1] - snt * x12[0];
	tempi = cst[0] * x11[1] - cst[1] * x11[0] - snt * x12[1];
	if (tempi == 0.) {
	    r__[r_dim1 + 1] = abs(tempr);
	    dt[0] = d_sign(&c_b4, &tempr);
	    dt[1] = 0.;
	} else {
	    r__[r_dim1 + 1] = dlapy2_(&tempr, &tempi);
	    dt[0] = tempr / r__[r_dim1 + 1];
	    dt[1] = -tempi / r__[r_dim1 + 1];
	}
    } else {

/*        Case op(M) = M. */

/*        Now form  X = v11*conjg( Qhat' ). */

	x11[0] = csq[0] * v1 - snq * v2[0];
	x11[1] = -csq[1] * v1 + snq * v2[1];
	x21[0] = -snq * v3;
	x12[0] = csq[0] * v2[0] + csq[1] * v2[1] + snq * v1;
	x12[1] = -csq[0] * v2[1] + csq[1] * v2[0];
	x22[0] = csq[0] * v3;
	x22[1] = csq[1] * v3;

/*        Obtain u11 from the QR-factorization of X. */

	sb03ov_(x11, x21, cst, &snt);
	r__[r_dim1 + 1] = x11[0];
	r__[(r_dim1 << 1) + 1] = cst[0] * x12[0] + cst[1] * x12[1] + snt * 
		x22[0];
	tempr = cst[0] * x22[0] - cst[1] * x22[1] - snt * x12[0];
	tempi = cst[0] * x22[1] + cst[1] * x22[0] - snt * x12[1];
	if (tempi == 0.) {
	    r__[(r_dim1 << 1) + 2] = abs(tempr);
	    dt[0] = d_sign(&c_b4, &tempr);
	    dt[1] = 0.;
	} else {
	    r__[(r_dim1 << 1) + 2] = dlapy2_(&tempr, &tempi);
	    dt[0] = tempr / r__[(r_dim1 << 1) + 2];
	    dt[1] = -tempi / r__[(r_dim1 << 1) + 2];
	}
    }

/*     The computations below are not needed when B and A are not */
/*     useful. Compute delta, eta and gamma as in (6.21) or (10.26). */

    if (y[0] == 0. && y[1] == 0.) {
	delta[0] = 0.;
	delta[1] = 0.;
	gamma[0] = 0.;
	gamma[1] = 0.;
	eta = alpha;
    } else {
	delta[0] = y[0] / v3;
	delta[1] = y[1] / v3;
	gamma[0] = -alpha * delta[0];
	gamma[1] = -alpha * delta[1];
	eta = p3 / v3;
	if (*discr) {
	    tempr = e1 * delta[0] - e2 * delta[1];
	    delta[1] = e1 * delta[1] + e2 * delta[0];
	    delta[0] = tempr;
	}
    }

    if (*ltrans) {

/*        Case op(M) = M'. */

/*        Find  X = conjg( That' )*( inv( v11 )*s11hat*v11 ). */
/*        ( Defer the scaling.) */

	x11[0] = cst[0] * e1 + cst[1] * e2;
	x11[1] = -cst[0] * e2 + cst[1] * e1;
	x21[0] = snt * e1;
	x21[1] = -snt * e2;
	x12[0] = sgn * (cst[0] * gamma[0] + cst[1] * gamma[1]) - snt * e1;
	x12[1] = sgn * (-cst[0] * gamma[1] + cst[1] * gamma[0]) - snt * e2;
	x22[0] = cst[0] * e1 + cst[1] * e2 + sgn * snt * gamma[0];
	x22[1] = cst[0] * e2 - cst[1] * e1 - sgn * snt * gamma[1];

/*        Now find  B = X*That. ( Include the scaling here.) */

	s[s_dim1 + 1] = cst[0] * x11[0] + cst[1] * x11[1] - snt * x12[0];
	tempr = cst[0] * x21[0] + cst[1] * x21[1] - snt * x22[0];
	tempi = cst[0] * x21[1] - cst[1] * x21[0] - snt * x22[1];
	s[s_dim1 + 2] = dt[0] * tempr - dt[1] * tempi;
	tempr = cst[0] * x12[0] - cst[1] * x12[1] + snt * x11[0];
	tempi = cst[0] * x12[1] + cst[1] * x12[0] + snt * x11[1];
	s[(s_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
	s[(s_dim1 << 1) + 2] = cst[0] * x22[0] - cst[1] * x22[1] + snt * x21[
		0];

/*        Form  X = ( inv( v11 )*p11 )*conjg( Phat' ). */

	tempr = dp[0] * eta;
	tempi = -dp[1] * eta;
	x11[0] = csp[0] * tempr - csp[1] * tempi + snp * delta[0];
	x11[1] = csp[0] * tempi + csp[1] * tempr - snp * delta[1];
	x21[0] = snp * alpha;
	x12[0] = -snp * tempr + csp[0] * delta[0] - csp[1] * delta[1];
	x12[1] = -snp * tempi - csp[0] * delta[1] - csp[1] * delta[0];
	x22[0] = csp[0] * alpha;
	x22[1] = -csp[1] * alpha;

/*        Finally form  A = conjg( That' )*X. */

	tempr = cst[0] * x11[0] - cst[1] * x11[1] - snt * x21[0];
	tempi = cst[0] * x22[1] + cst[1] * x22[0];
	a[a_dim1 + 1] = dt[0] * tempr + dt[1] * tempi;
	tempr = cst[0] * x12[0] - cst[1] * x12[1] - snt * x22[0];
	tempi = cst[0] * x12[1] + cst[1] * x12[0] - snt * x22[0];
	a[(a_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
	a[a_dim1 + 2] = 0.;
	a[(a_dim1 << 1) + 2] = cst[0] * x22[0] + cst[1] * x22[1] + snt * x12[
		0];
    } else {

/*        Case op(M) = M. */

/*        Find  X = That*( v11*s11hat*inv( v11 ) ). ( Defer the scaling.) */

	x11[0] = cst[0] * e1 + cst[1] * e2;
	x11[1] = cst[0] * e2 - cst[1] * e1;
	x21[0] = -snt * e1;
	x21[1] = -snt * e2;
	x12[0] = sgn * (cst[0] * gamma[0] - cst[1] * gamma[1]) + snt * e1;
	x12[1] = sgn * (-cst[0] * gamma[1] - cst[1] * gamma[0]) - snt * e2;
	x22[0] = cst[0] * e1 + cst[1] * e2 - sgn * snt * gamma[0];
	x22[1] = -cst[0] * e2 + cst[1] * e1 + sgn * snt * gamma[1];

/*        Now find  B = X*conjg( That' ). ( Include the scaling here.) */

	s[s_dim1 + 1] = cst[0] * x11[0] - cst[1] * x11[1] + snt * x12[0];
	tempr = cst[0] * x21[0] - cst[1] * x21[1] + snt * x22[0];
	tempi = cst[0] * x21[1] + cst[1] * x21[0] + snt * x22[1];
	s[s_dim1 + 2] = dt[0] * tempr - dt[1] * tempi;
	tempr = cst[0] * x12[0] + cst[1] * x12[1] - snt * x11[0];
	tempi = cst[0] * x12[1] - cst[1] * x12[0] - snt * x11[1];
	s[(s_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
	s[(s_dim1 << 1) + 2] = cst[0] * x22[0] + cst[1] * x22[1] - snt * x21[
		0];

/*        Form  X = Phat*( p11*inv( v11 ) ). */

	tempr = dp[0] * eta;
	tempi = -dp[1] * eta;
	x11[0] = csp[0] * alpha;
	x11[1] = csp[1] * alpha;
	x21[0] = snp * alpha;
	x12[0] = csp[0] * delta[0] + csp[1] * delta[1] - snp * tempr;
	x12[1] = -csp[0] * delta[1] + csp[1] * delta[0] - snp * tempi;
	x22[0] = csp[0] * tempr + csp[1] * tempi + snp * delta[0];
	x22[1] = csp[0] * tempi - csp[1] * tempr - snp * delta[1];

/*        Finally form  A = X*conjg( That' ). */

	a[a_dim1 + 1] = cst[0] * x11[0] - cst[1] * x11[1] + snt * x12[0];
	a[a_dim1 + 2] = 0.;
	a[(a_dim1 << 1) + 1] = cst[0] * x12[0] + cst[1] * x12[1] - snt * x11[
		0];
	tempr = cst[0] * x22[0] + cst[1] * x22[1] - snt * x21[0];
	tempi = cst[0] * x22[1] - cst[1] * x22[0];
	a[(a_dim1 << 1) + 2] = dt[0] * tempr + dt[1] * tempi;
    }

    if (*scale != 1.) {
	a[a_dim1 + 1] = *scale * a[a_dim1 + 1];
	a[(a_dim1 << 1) + 1] = *scale * a[(a_dim1 << 1) + 1];
	a[(a_dim1 << 1) + 2] = *scale * a[(a_dim1 << 1) + 2];
    }

    return 0;
/* *** Last line of SB03OY *** */
} /* sb03oy_ */

