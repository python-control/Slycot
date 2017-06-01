/* SB03OR.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;

/* Subroutine */ int sb03or_(logical *discr, logical *ltrans, integer *n, 
	integer *m, doublereal *s, integer *lds, doublereal *a, integer *lda, 
	doublereal *c__, integer *ldc, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, s_dim1, s_offset, i__1, i__2;

    /* Local variables */
    static integer j, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer l1, l2;
    static doublereal g11, g12, g21, g22;
    static integer dl;
    static doublereal at[4]	/* was [2][2] */, vec[4]	/* was [2][2] 
	    */;
    static integer l2p1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer isgn;
    static logical tbyt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer infom;
    extern /* Subroutine */ int sb04px_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlasy2_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To compute the solution of the Sylvester equations */

/*        op(S)'*X + X*op(A) = scale*C, if DISCR = .FALSE.  or */

/*        op(S)'*X*op(A) - X = scale*C, if DISCR = .TRUE. */

/*     where op(K) = K or K' (i.e., the transpose of the matrix K), S is */
/*     an N-by-N block upper triangular matrix with one-by-one and */
/*     two-by-two blocks on the diagonal, A is an M-by-M matrix (M = 1 or */
/*     M = 2), X and C are each N-by-M matrices, and scale is an output */
/*     scale factor, set less than or equal to 1 to avoid overflow in X. */
/*     The solution X is overwritten on C. */

/*     SB03OR  is a service routine for the Lyapunov solver  SB03OT. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the equation to be solved: */
/*             = .FALSE.:  op(S)'*X + X*op(A) = scale*C; */
/*             = .TRUE. :  op(S)'*X*op(A) - X = scale*C. */

/*     LTRANS  LOGICAL */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = .FALSE.:  op(K) = K    (No transpose); */
/*             = .TRUE. :  op(K) = K**T (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix  S  and also the number of rows of */
/*             matrices  X and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix  A  and also the number of columns */
/*             of matrices  X and C.  M = 1 or M = 2. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading  N-by-N  upper Hessenberg part of the array  S */
/*             must contain the block upper triangular matrix. The */
/*             elements below the upper Hessenberg part of the array  S */
/*             are not referenced.  The array  S  must not contain */
/*             diagonal blocks larger than two-by-two and the two-by-two */
/*             blocks must only correspond to complex conjugate pairs of */
/*             eigenvalues, not to real eigenvalues. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDS,M) */
/*             The leading  M-by-M  part of this array must contain a */
/*             given matrix, where M = 1 or M = 2. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= M. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             On entry, C must contain an N-by-M matrix, where M = 1 or */
/*             M = 2. */
/*             On exit, C contains the N-by-M matrix X, the solution of */
/*             the Sylvester equation. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if DISCR = .FALSE., and S and -A have common */
/*                   eigenvalues, or if DISCR = .TRUE., and S and A have */
/*                   eigenvalues whose product is equal to unity; */
/*                   a solution has been computed using slightly */
/*                   perturbed values. */

/*     METHOD */

/*     The LAPACK scheme for solving Sylvester equations is adapted. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(N M) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routines SB03CW and SB03CX by */
/*     Sven Hammarling, NAG Ltd, United Kingdom, Oct. 1986. */
/*     Partly based on routine PLYAP4 by A. Varga, University of Bochum, */
/*     May 1992. */

/*     REVISIONS */

/*     December 1997, April 1998, May 1999, April 2000. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -3;
    } else if (! (*m == 1 || *m == 2)) {
	*info = -4;
    } else if (*lds < max(1,*n)) {
	*info = -6;
    } else if (*lda < *m) {
	*info = -8;
    } else if (*ldc < max(1,*n)) {
	*info = -10;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB03OR", &i__1, (ftnlen)6);
	return 0;
    }

    *scale = 1.;

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    isgn = 1;
    tbyt = *m == 2;
    infom = 0;

/*     Construct A'. */

    at[0] = a[a_dim1 + 1];
    if (tbyt) {
	at[2] = a[a_dim1 + 2];
	at[1] = a[(a_dim1 << 1) + 1];
	at[3] = a[(a_dim1 << 1) + 2];
    }

    if (*ltrans) {

/*        Start row loop (index = L). */
/*        L1 (L2) : row index of the first (last) row of X(L). */

	lnext = *n;

	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L20;
	    }
	    l1 = l;
	    l2 = l;
	    if (l > 1) {
		if (s[l + (l - 1) * s_dim1] != 0.) {
		    --l1;
		}
		lnext = l1 - 1;
	    }
	    dl = l2 - l1 + 1;
/* Computing MIN */
	    i__1 = l2 + 1;
	    l2p1 = min(i__1,*n);

	    if (*discr) {

/*              Solve  S*X*A' - X = scale*C. */

/*              The L-th block of X is determined from */

/*              S(L,L)*X(L)*A' - X(L) = C(L) - R(L), */

/*              where */

/*                      N */
/*              R(L) = SUM [S(L,J)*X(J)] * A' . */
/*                    J=L+1 */

		i__1 = *n - l2;
		g11 = -ddot_(&i__1, &s[l1 + l2p1 * s_dim1], lds, &c__[l2p1 + 
			c_dim1], &c__1);
		if (tbyt) {
		    i__1 = *n - l2;
		    g12 = -ddot_(&i__1, &s[l1 + l2p1 * s_dim1], lds, &c__[
			    l2p1 + (c_dim1 << 1)], &c__1);
		    vec[0] = c__[l1 + c_dim1] + g11 * at[0] + g12 * at[1];
		    vec[2] = c__[l1 + (c_dim1 << 1)] + g11 * at[2] + g12 * at[
			    3];
		} else {
		    vec[0] = c__[l1 + c_dim1] + g11 * at[0];
		}
		if (dl != 1) {
		    i__1 = *n - l2;
		    g21 = -ddot_(&i__1, &s[l2 + l2p1 * s_dim1], lds, &c__[
			    l2p1 + c_dim1], &c__1);
		    if (tbyt) {
			i__1 = *n - l2;
			g22 = -ddot_(&i__1, &s[l2 + l2p1 * s_dim1], lds, &c__[
				l2p1 + (c_dim1 << 1)], &c__1);
			vec[1] = c__[l2 + c_dim1] + g21 * at[0] + g22 * at[1];
			vec[3] = c__[l2 + (c_dim1 << 1)] + g21 * at[2] + g22 *
				 at[3];
		    } else {
			vec[1] = c__[l2 + c_dim1] + g21 * at[0];
		    }
		}
		i__1 = -isgn;
		sb04px_(&c_false, &c_false, &i__1, &dl, m, &s[l1 + l1 * 
			s_dim1], lds, at, &c__2, vec, &c__2, &scaloc, x, &
			c__2, &xnorm, info);
	    } else {

/*              Solve  S*X + X*A' = scale*C. */

/*              The L-th block of X is determined from */

/*              S(L,L)*X(L) + X(L)*A' = C(L) - R(L), */

/*              where */
/*                       N */
/*              R(L) =  SUM S(L,J)*X(J) . */
/*                     J=L+1 */

		i__1 = *n - l2;
		vec[0] = c__[l1 + c_dim1] - ddot_(&i__1, &s[l1 + l2p1 * 
			s_dim1], lds, &c__[l2p1 + c_dim1], &c__1);
		if (tbyt) {
		    i__1 = *n - l2;
		    vec[2] = c__[l1 + (c_dim1 << 1)] - ddot_(&i__1, &s[l1 + 
			    l2p1 * s_dim1], lds, &c__[l2p1 + (c_dim1 << 1)], &
			    c__1);
		}

		if (dl != 1) {
		    i__1 = *n - l2;
		    vec[1] = c__[l2 + c_dim1] - ddot_(&i__1, &s[l2 + l2p1 * 
			    s_dim1], lds, &c__[l2p1 + c_dim1], &c__1);
		    if (tbyt) {
			i__1 = *n - l2;
			vec[3] = c__[l2 + (c_dim1 << 1)] - ddot_(&i__1, &s[l2 
				+ l2p1 * s_dim1], lds, &c__[l2p1 + (c_dim1 << 
				1)], &c__1);
		    }
		}
		dlasy2_(&c_false, &c_false, &isgn, &dl, m, &s[l1 + l1 * 
			s_dim1], lds, at, &c__2, vec, &c__2, &scaloc, x, &
			c__2, &xnorm, info);
	    }
	    infom = max(*info,infom);
	    if (scaloc != 1.) {

		i__1 = *m;
		for (j = 1; j <= i__1; ++j) {
		    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L10: */
		}

		*scale *= scaloc;
	    }
	    c__[l1 + c_dim1] = x[0];
	    if (tbyt) {
		c__[l1 + (c_dim1 << 1)] = x[2];
	    }
	    if (dl != 1) {
		c__[l2 + c_dim1] = x[1];
		if (tbyt) {
		    c__[l2 + (c_dim1 << 1)] = x[3];
		}
	    }
L20:
	    ;
	}

    } else {

/*        Start row loop (index = L). */
/*        L1 (L2) : row index of the first (last) row of X(L). */

	lnext = 1;

	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    if (l < lnext) {
		goto L40;
	    }
	    l1 = l;
	    l2 = l;
	    if (l < *n) {
		if (s[l + 1 + l * s_dim1] != 0.) {
		    ++l2;
		}
		lnext = l2 + 1;
	    }
	    dl = l2 - l1 + 1;

	    if (*discr) {

/*              Solve  A'*X'*S - X' = scale*C'. */

/*              The L-th block of X is determined from */

/*              A'*X(L)'*S(L,L) - X(L)' = C(L)' - R(L), */

/*              where */

/*                          L-1 */
/*              R(L) = A' * SUM [X(J)'*S(J,L)] . */
/*                          J=1 */

		i__2 = l1 - 1;
		g11 = -ddot_(&i__2, &c__[c_offset], &c__1, &s[l1 * s_dim1 + 1]
			, &c__1);
		if (tbyt) {
		    i__2 = l1 - 1;
		    g21 = -ddot_(&i__2, &c__[(c_dim1 << 1) + 1], &c__1, &s[l1 
			    * s_dim1 + 1], &c__1);
		    vec[0] = c__[l1 + c_dim1] + at[0] * g11 + at[2] * g21;
		    vec[1] = c__[l1 + (c_dim1 << 1)] + at[1] * g11 + at[3] * 
			    g21;
		} else {
		    vec[0] = c__[l1 + c_dim1] + at[0] * g11;
		}
		if (dl != 1) {
		    i__2 = l1 - 1;
		    g12 = -ddot_(&i__2, &c__[c_offset], &c__1, &s[l2 * s_dim1 
			    + 1], &c__1);
		    if (tbyt) {
			i__2 = l1 - 1;
			g22 = -ddot_(&i__2, &c__[(c_dim1 << 1) + 1], &c__1, &
				s[l2 * s_dim1 + 1], &c__1);
			vec[2] = c__[l2 + c_dim1] + at[0] * g12 + at[2] * g22;
			vec[3] = c__[l2 + (c_dim1 << 1)] + at[1] * g12 + at[3]
				 * g22;
		    } else {
			vec[2] = c__[l2 + c_dim1] + at[0] * g12;
		    }
		}
		i__2 = -isgn;
		sb04px_(&c_false, &c_false, &i__2, m, &dl, at, &c__2, &s[l1 + 
			l1 * s_dim1], lds, vec, &c__2, &scaloc, x, &c__2, &
			xnorm, info);
	    } else {

/*              Solve  A'*X' + X'*S = scale*C'. */

/*              The L-th block of X is determined from */

/*              A'*X(L)' + X(L)'*S(L,L) = C(L)' - R(L), */

/*              where */
/*                     L-1 */
/*              R(L) = SUM [X(J)'*S(J,L)]. */
/*                     J=1 */

		i__2 = l1 - 1;
		vec[0] = c__[l1 + c_dim1] - ddot_(&i__2, &c__[c_offset], &
			c__1, &s[l1 * s_dim1 + 1], &c__1);
		if (tbyt) {
		    i__2 = l1 - 1;
		    vec[1] = c__[l1 + (c_dim1 << 1)] - ddot_(&i__2, &c__[(
			    c_dim1 << 1) + 1], &c__1, &s[l1 * s_dim1 + 1], &
			    c__1);
		}

		if (dl != 1) {
		    i__2 = l1 - 1;
		    vec[2] = c__[l2 + c_dim1] - ddot_(&i__2, &c__[c_offset], &
			    c__1, &s[l2 * s_dim1 + 1], &c__1);
		    if (tbyt) {
			i__2 = l1 - 1;
			vec[3] = c__[l2 + (c_dim1 << 1)] - ddot_(&i__2, &c__[(
				c_dim1 << 1) + 1], &c__1, &s[l2 * s_dim1 + 1],
				 &c__1);
		    }
		}
		dlasy2_(&c_false, &c_false, &isgn, m, &dl, at, &c__2, &s[l1 + 
			l1 * s_dim1], lds, vec, &c__2, &scaloc, x, &c__2, &
			xnorm, info);
	    }
	    infom = max(*info,infom);
	    if (scaloc != 1.) {

		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
		    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
/* L30: */
		}

		*scale *= scaloc;
	    }
	    c__[l1 + c_dim1] = x[0];
	    if (tbyt) {
		c__[l1 + (c_dim1 << 1)] = x[1];
	    }
	    if (dl != 1) {
		c__[l2 + c_dim1] = x[2];
		if (tbyt) {
		    c__[l2 + (c_dim1 << 1)] = x[3];
		}
	    }
L40:
	    ;
	}
    }

    *info = infom;
    return 0;
/* *** Last line of SB03OR *** */
} /* sb03or_ */

