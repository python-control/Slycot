/* MB04OW.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04ow_(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *t, integer *ldt, doublereal *x, integer *
	incx, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *incd)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, t_dim1, 
	    t_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal ci;
    static integer mn;
    static doublereal si;
    static integer ix;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


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

/*     To perform the QR factorization */

/*        ( U  ) = Q*( R ),  where  U = ( U1  U2 ),  R = ( R1  R2 ), */
/*        ( x' )     ( 0 )              ( 0   T  )       ( 0   R3 ) */

/*     where U and R are (m+n)-by-(m+n) upper triangular matrices, x is */
/*     an m+n element vector, U1 is m-by-m, T is n-by-n, stored */
/*     separately, and Q is an (m+n+1)-by-(m+n+1) orthogonal matrix. */

/*     The matrix ( U1 U2 ) must be supplied in the m-by-(m+n) upper */
/*     trapezoidal part of the array A and this is overwritten by the */
/*     corresponding part ( R1 R2 ) of R. The remaining upper triangular */
/*     part of R, R3, is overwritten on the array T. */

/*     The transformations performed are also applied to the (m+n+1)-by-p */
/*     matrix ( B' C' d )' (' denotes transposition), where B, C, and d' */
/*     are m-by-p, n-by-p, and 1-by-p matrices, respectively. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M      (input) INTEGER */
/*            The number of rows of the matrix ( U1  U2 ).  M >= 0. */

/*     N      (input) INTEGER */
/*            The order of the matrix T.  N >= 0. */

/*     P      (input) INTEGER */
/*            The number of columns of the matrices B and C.  P >= 0. */

/*     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*            On entry, the leading M-by-(M+N) upper trapezoidal part of */
/*            this array must contain the upper trapezoidal matrix */
/*            ( U1 U2 ). */
/*            On exit, the leading M-by-(M+N) upper trapezoidal part of */
/*            this array contains the upper trapezoidal matrix ( R1 R2 ). */
/*            The strict lower triangle of A is not referenced. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,M). */

/*     T      (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*            On entry, the leading N-by-N upper triangular part of this */
/*            array must contain the upper triangular matrix T. */
/*            On exit, the leading N-by-N upper triangular part of this */
/*            array contains the upper triangular matrix R3. */
/*            The strict lower triangle of T is not referenced. */

/*     LDT    INTEGER */
/*            The leading dimension of the array T.  LDT >= max(1,N). */

/*     X      (input/output) DOUBLE PRECISION array, dimension */
/*            (1+(M+N-1)*INCX), if M+N > 0, or dimension (0), if M+N = 0. */
/*            On entry, the incremented array X must contain the */
/*            vector x. On exit, the content of X is changed. */

/*     INCX   (input) INTEGER */
/*            Specifies the increment for the elements of X.  INCX > 0. */

/*     B      (input/output) DOUBLE PRECISION array, dimension (LDB,P) */
/*            On entry, the leading M-by-P part of this array must */
/*            contain the matrix B. */
/*            On exit, the leading M-by-P part of this array contains */
/*            the transformed matrix B. */
/*            If M = 0 or P = 0, this array is not referenced. */

/*     LDB    INTEGER */
/*            The leading dimension of the array B. */
/*            LDB >= max(1,M), if P > 0; */
/*            LDB >= 1,        if P = 0. */

/*     C      (input/output) DOUBLE PRECISION array, dimension (LDC,P) */
/*            On entry, the leading N-by-P part of this array must */
/*            contain the matrix C. */
/*            On exit, the leading N-by-P part of this array contains */
/*            the transformed matrix C. */
/*            If N = 0 or P = 0, this array is not referenced. */

/*     LDC    INTEGER */
/*            The leading dimension of the array C. */
/*            LDC >= max(1,N), if P > 0; */
/*            LDC >= 1,        if P = 0. */

/*     D      (input/output) DOUBLE PRECISION array, dimension */
/*            (1+(P-1)*INCD), if P > 0, or dimension (0), if P = 0. */
/*            On entry, the incremented array D must contain the */
/*            vector d. */
/*            On exit, this incremented array contains the transformed */
/*            vector d. */
/*            If P = 0, this array is not referenced. */

/*     INCD   (input) INTEGER */
/*            Specifies the increment for the elements of D.  INCD > 0. */

/*     METHOD */

/*     Let q = m+n. The matrix Q is formed as a sequence of plane */
/*     rotations in planes (1, q+1), (2, q+1), ..., (q, q+1), the */
/*     rotation in the (j, q+1)th plane, Q(j), being chosen to */
/*     annihilate the jth element of x. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0((M+N)*(M+N+P)) operations and is backward */
/*     stable. */

/*     FURTHER COMMENTS */

/*     For P = 0, this routine produces the same result as SLICOT Library */
/*     routine MB04OX, but matrix T may not be stored in the array A. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --x;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --d__;

    /* Function Body */
    mn = *m + *n;
    if (*incx > 1) {

/*        Code for increment INCX > 1. */

	ix = 1;
	if (*m > 0) {

	    i__1 = *m - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dlartg_(&a[i__ + i__ * a_dim1], &x[ix], &ci, &si, &temp);
		a[i__ + i__ * a_dim1] = temp;
		ix += *incx;
		i__2 = mn - i__;
		drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &x[ix], incx, 
			&ci, &si);
		if (*p > 0) {
		    drot_(p, &b[i__ + b_dim1], ldb, &d__[1], incd, &ci, &si);
		}
/* L10: */
	    }

	    dlartg_(&a[*m + *m * a_dim1], &x[ix], &ci, &si, &temp);
	    a[*m + *m * a_dim1] = temp;
	    ix += *incx;
	    if (*n > 0) {
		drot_(n, &a[*m + (*m + 1) * a_dim1], lda, &x[ix], incx, &ci, &
			si);
	    }
	    if (*p > 0) {
		drot_(p, &b[*m + b_dim1], ldb, &d__[1], incd, &ci, &si);
	    }
	}

	if (*n > 0) {

	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dlartg_(&t[i__ + i__ * t_dim1], &x[ix], &ci, &si, &temp);
		t[i__ + i__ * t_dim1] = temp;
		ix += *incx;
		i__2 = *n - i__;
		drot_(&i__2, &t[i__ + (i__ + 1) * t_dim1], ldt, &x[ix], incx, 
			&ci, &si);
		if (*p > 0) {
		    drot_(p, &c__[i__ + c_dim1], ldc, &d__[1], incd, &ci, &si)
			    ;
		}
/* L20: */
	    }

	    dlartg_(&t[*n + *n * t_dim1], &x[ix], &ci, &si, &temp);
	    t[*n + *n * t_dim1] = temp;
	    if (*p > 0) {
		drot_(p, &c__[*n + c_dim1], ldc, &d__[1], incd, &ci, &si);
	    }
	}

    } else if (*incx == 1) {

/*        Code for increment INCX = 1. */

	if (*m > 0) {

	    i__1 = *m - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dlartg_(&a[i__ + i__ * a_dim1], &x[i__], &ci, &si, &temp);
		a[i__ + i__ * a_dim1] = temp;
		i__2 = mn - i__;
		drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &x[i__ + 1], &
			c__1, &ci, &si);
		if (*p > 0) {
		    drot_(p, &b[i__ + b_dim1], ldb, &d__[1], incd, &ci, &si);
		}
/* L30: */
	    }

	    dlartg_(&a[*m + *m * a_dim1], &x[*m], &ci, &si, &temp);
	    a[*m + *m * a_dim1] = temp;
	    if (*n > 0) {
		drot_(n, &a[*m + (*m + 1) * a_dim1], lda, &x[*m + 1], &c__1, &
			ci, &si);
	    }
	    if (*p > 0) {
		drot_(p, &b[*m + b_dim1], ldb, &d__[1], incd, &ci, &si);
	    }
	}

	if (*n > 0) {
	    ix = *m + 1;

	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dlartg_(&t[i__ + i__ * t_dim1], &x[ix], &ci, &si, &temp);
		t[i__ + i__ * t_dim1] = temp;
		++ix;
		i__2 = *n - i__;
		drot_(&i__2, &t[i__ + (i__ + 1) * t_dim1], ldt, &x[ix], &c__1,
			 &ci, &si);
		if (*p > 0) {
		    drot_(p, &c__[i__ + c_dim1], ldc, &d__[1], incd, &ci, &si)
			    ;
		}
/* L40: */
	    }

	    dlartg_(&t[*n + *n * t_dim1], &x[ix], &ci, &si, &temp);
	    t[*n + *n * t_dim1] = temp;
	    if (*p > 0) {
		drot_(p, &c__[*n + c_dim1], ldc, &d__[1], incd, &ci, &si);
	    }
	}
    }

    return 0;
/* *** Last line of MB04OW *** */
} /* mb04ow_ */

