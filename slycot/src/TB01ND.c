/* TB01ND.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;

/* Subroutine */ int tb01nd_(char *jobu, char *uplo, integer *n, integer *p, 
	doublereal *a, integer *lda, doublereal *c__, integer *ldc, 
	doublereal *u, integer *ldu, doublereal *dwork, integer *info, ftnlen 
	jobu_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, u_dim1, u_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer j, n1, p1, ii, nj;
    static doublereal dz;
    static integer par1, par2, par3, par4, par5, par6;
    static logical ljoba, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical luplo;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlatzm_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);


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

/*     To reduce the pair (A,C) to lower or upper observer Hessenberg */
/*     form using (and optionally accumulating) unitary state-space */
/*     transformations. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the unitary state-space transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is initialized to the unit matrix and the */
/*                     unitary transformation matrix U is returned; */
/*             = 'U':  The given matrix U is updated by the unitary */
/*                     transformations used in the reduction. */

/*     UPLO    CHARACTER*1 */
/*             Indicates whether the user wishes the pair (A,C) to be */
/*             reduced to upper or lower observer Hessenberg form as */
/*             follows: */
/*             = 'U':  Upper observer Hessenberg form; */
/*             = 'L':  Lower observer Hessenberg form. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the */
/*             matrix A.  N >= 0. */

/*     P       (input) INTEGER */
/*             The actual output dimension, i.e. the number of rows of */
/*             the matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A to be transformed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state transition matrix U' * A * U. */
/*             The annihilated elements are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C to be transformed. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output matrix C * U. */
/*             The annihilated elements are set to zero. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, if JOBU = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix U (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix U and the state-space transformation */
/*             matrix which reduces the given pair to observer Hessenberg */
/*             form. */
/*             On exit, if JOBU = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated unitary */
/*             similarity transformations which reduces the given pair */
/*             to observer Hessenberg form. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. If JOBU = 'U' or */
/*             JOBU = 'I', LDU >= MAX(1,N); if JOBU = 'N', LDU >= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(N,P-1)) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a unitary state-space transformation U, which */
/*     reduces the pair (A,C) to one of the following observer Hessenberg */
/*     forms: */

/*                                N */
/*                       |*  . . . . . .  *| */
/*                       |.               .| */
/*                       |.               .| */
/*                       |.               .| N */
/*                       |*               .| */
/*            |U'AU|     |   .            .| */
/*            |----|  =  |     .          .| */
/*            |CU  |     |       * . . .  *| */
/*                       ------------------- */
/*                       |         * . .  *| */
/*                       |           .    .| P */
/*                       |             .  .| */
/*                       |                *| */

/*         if UPLO = 'U', or */

/*                               N */
/*                      |*                | */
/*                      |.  .             | */
/*                      |.    .           | P */
/*                      |*  . . *         | */
/*            |CU  |    ------------------- */
/*            |----|  = |*  . . . *       | */
/*            |U'AU|    |.          .     | */
/*                      |.            .   | */
/*                      |.               *| */
/*                      |.               .| N */
/*                      |.               .| */
/*                      |.               .| */
/*                      |*  . . . . . .  *| */

/*     if UPLO = 'L'. */

/*     If P >= N, then the matrix CU is trapezoidal and U'AU is full. */

/*     REFERENCES */

/*     [1] Van Dooren, P. and Verhaegen, M.H.G. */
/*         On the use of unitary state-space transformations. */
/*         In : Contemporary Mathematics on Linear Algebra and its Role */
/*         in Systems Theory, 47, AMS, Providence, 1985. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O((N + P) x N**2) operations and is */
/*     backward stable (see [1]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01BD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     February 1997. */

/*     KEYWORDS */

/*     Controllability, observer Hessenberg form, orthogonal */
/*     transformation, unitary transformation. */

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
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ljobi = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
    ljoba = ljobi || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljoba && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldc < max(1,*p)) {
	*info = -8;
    } else if (! ljoba && *ldu < 1 || ljoba && *ldu < max(1,*n)) {
	*info = -10;
    }

    if (*info != 0) {

/*        Error return */

	i__1 = -(*info);
	xerbla_("TB01ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *p == 0) {
	return 0;
    }

    p1 = *p + 1;
    n1 = *n - 1;

    if (ljobi) {

/*        Initialize U to the identity matrix. */

	dlaset_("Full", n, n, &c_b9, &c_b10, &u[u_offset], ldu, (ftnlen)4);
    }

/*     Perform transformations involving both C and A. */

    i__1 = min(*p,n1);
    for (j = 1; j <= i__1; ++j) {
	nj = *n - j;
	if (luplo) {
	    par1 = *p - j + 1;
	    par2 = nj + 1;
	    par3 = 1;
	    par4 = *p - j;
	    par5 = nj;
	} else {
	    par1 = j;
	    par2 = j;
	    par3 = j + 1;
	    par4 = *p;
	    par5 = *n;
	}

	i__2 = nj + 1;
	dlarfg_(&i__2, &c__[par1 + par2 * c_dim1], &c__[par1 + par3 * c_dim1],
		 ldc, &dz);

/*        Update A. */

	i__2 = nj + 1;
	dlatzm_("Left", &i__2, n, &c__[par1 + par3 * c_dim1], ldc, &dz, &a[
		par2 + a_dim1], &a[par3 + a_dim1], lda, &dwork[1], (ftnlen)4);
	i__2 = nj + 1;
	dlatzm_("Right", n, &i__2, &c__[par1 + par3 * c_dim1], ldc, &dz, &a[
		par2 * a_dim1 + 1], &a[par3 * a_dim1 + 1], lda, &dwork[1], (
		ftnlen)5);

	if (ljoba) {

/*           Update U. */

	    i__2 = nj + 1;
	    dlatzm_("Right", n, &i__2, &c__[par1 + par3 * c_dim1], ldc, &dz, &
		    u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &dwork[
		    1], (ftnlen)5);
	}

	if (j != *p) {

/*           Update C. */

	    i__2 = par4 - par3 + 1;
	    i__3 = nj + 1;
	    dlatzm_("Right", &i__2, &i__3, &c__[par1 + par3 * c_dim1], ldc, &
		    dz, &c__[par3 + par2 * c_dim1], &c__[par3 + par3 * c_dim1]
		    , ldc, &dwork[1], (ftnlen)5);
	}

	i__2 = par5;
	for (ii = par3; ii <= i__2; ++ii) {
	    c__[par1 + ii * c_dim1] = 0.;
/* L10: */
	}

/* L20: */
    }

    i__1 = n1;
    for (j = p1; j <= i__1; ++j) {

/*        Perform next transformations only involving A. */

	nj = *n - j;
	if (luplo) {
	    par1 = *n + p1 - j;
	    par2 = nj + 1;
	    par3 = 1;
	    par4 = nj;
	    par5 = 1;
	    par6 = *n + *p - j;
	} else {
	    par1 = j - *p;
	    par2 = j;
	    par3 = j + 1;
	    par4 = *n;
	    par5 = j - *p + 1;
	    par6 = *n;
	}

	if (nj > 0) {

	    i__2 = nj + 1;
	    dlarfg_(&i__2, &a[par1 + par2 * a_dim1], &a[par1 + par3 * a_dim1],
		     lda, &dz);

/*           Update A. */

	    i__2 = nj + 1;
	    dlatzm_("Left", &i__2, n, &a[par1 + par3 * a_dim1], lda, &dz, &a[
		    par2 + a_dim1], &a[par3 + a_dim1], lda, &dwork[1], (
		    ftnlen)4);
	    i__2 = par6 - par5 + 1;
	    i__3 = nj + 1;
	    dlatzm_("Right", &i__2, &i__3, &a[par1 + par3 * a_dim1], lda, &dz,
		     &a[par5 + par2 * a_dim1], &a[par5 + par3 * a_dim1], lda, 
		    &dwork[1], (ftnlen)5);

	    if (ljoba) {

/*              Update U. */

		i__2 = nj + 1;
		dlatzm_("Right", n, &i__2, &a[par1 + par3 * a_dim1], lda, &dz,
			 &u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &
			dwork[1], (ftnlen)5);
	    }

	    i__2 = par4;
	    for (ii = par3; ii <= i__2; ++ii) {
		a[par1 + ii * a_dim1] = 0.;
/* L30: */
	    }

	}

/* L40: */
    }

    return 0;
/* *** Last line of TB01ND *** */
} /* tb01nd_ */

