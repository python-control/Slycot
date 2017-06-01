/* MB03YA.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb03ya_(logical *wantt, logical *wantq, logical *wantz, 
	integer *n, integer *ilo, integer *ihi, integer *iloq, integer *ihiq, 
	integer *pos, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;

    /* Local variables */
    static integer j, i1, i2;
    static doublereal cs;
    static integer nq;
    static doublereal sn, temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
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

/*     To annihilate one or two entries on the subdiagonal of the */
/*     Hessenberg matrix A for dealing with zero elements on the diagonal */
/*     of the triangular matrix B. */

/*     MB03YA is an auxiliary routine called by SLICOT Library routines */
/*     MB03XP and MB03YD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WANTT   LOGICAL */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = .TRUE. :  Compute the full Schur form; */
/*             = .FALSE.:  compute the eigenvalues only. */

/*     WANTQ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = .TRUE. :  The matrix Q is updated; */
/*             = .FALSE.:  the matrix Q is not required. */

/*     WANTZ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = .TRUE. :  The matrix Z is updated; */
/*             = .FALSE.:  the matrix Z is not required. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B. N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that the matrices A and B are already */
/*             (quasi) upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N. The routine works primarily with the submatrices */
/*             in rows and columns ILO to IHI, but applies the */
/*             transformations to all the rows and columns of the */
/*             matrices A and B, if WANTT = .TRUE.. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     ILOQ    (input) INTEGER */
/*     IHIQ    (input) INTEGER */
/*             Specify the rows of Q and Z to which transformations */
/*             must be applied if WANTQ = .TRUE. and WANTZ = .TRUE., */
/*             respectively. */
/*             1 <= ILOQ <= ILO; IHI <= IHIQ <= N. */

/*     POS     (input) INTEGER */
/*             The position of the zero element on the diagonal of B. */
/*             ILO <= POS <= IHI. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper Hessenberg matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the updated matrix A where A(POS,POS-1) = 0, if POS > ILO, */
/*             and A(POS+1,POS) = 0, if POS < IHI. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain an upper triangular matrix B with B(POS,POS) = 0. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the updated upper triangular matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Q of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Q updated in the */
/*             submatrix Q(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTQ = .FALSE., Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= 1. */
/*             If WANTQ = .TRUE., LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Z of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Z updated in the */
/*             submatrix Z(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTZ = .FALSE., Z is not referenced. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= 1. */
/*             If WANTZ = .TRUE., LDZ >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The method is illustrated by Wilkinson diagrams for N = 5, */
/*     POS = 3: */

/*           [ x x x x x ]       [ x x x x x ] */
/*           [ x x x x x ]       [ o x x x x ] */
/*       A = [ o x x x x ],  B = [ o o o x x ]. */
/*           [ o o x x x ]       [ o o o x x ] */
/*           [ o o o x x ]       [ o o o o x ] */

/*     First, a QR factorization is applied to A(1:3,1:3) and the */
/*     resulting nonzero in the updated matrix B is immediately */
/*     annihilated by a Givens rotation acting on columns 1 and 2: */

/*           [ x x x x x ]       [ x x x x x ] */
/*           [ x x x x x ]       [ o x x x x ] */
/*       A = [ o o x x x ],  B = [ o o o x x ]. */
/*           [ o o x x x ]       [ o o o x x ] */
/*           [ o o o x x ]       [ o o o o x ] */

/*     Secondly, an RQ factorization is applied to A(4:5,4:5) and the */
/*     resulting nonzero in the updated matrix B is immediately */
/*     annihilated by a Givens rotation acting on rows 4 and 5: */

/*           [ x x x x x ]       [ x x x x x ] */
/*           [ x x x x x ]       [ o x x x x ] */
/*       A = [ o o x x x ],  B = [ o o o x x ]. */
/*           [ o o o x x ]       [ o o o x x ] */
/*           [ o o o x x ]       [ o o o o x ] */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P. */
/*         The periodic Schur decomposition: Algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**2) floating point operations and is */
/*     backward stable. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLADFB). */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    *info = 0;
    nq = *ihiq - *iloq + 1;
    if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -5;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -6;
    } else if (*iloq < 1 || *iloq > *ilo) {
	*info = -7;
    } else if (*ihiq < *ihi || *ihiq > *n) {
	*info = -8;
    } else if (*pos < *ilo || *pos > *ihi) {
	*info = -9;
    } else if (*lda < max(1,*n)) {
	*info = -11;
    } else if (*ldb < max(1,*n)) {
	*info = -13;
    } else if (*ldq < 1 || *wantq && *ldq < *n) {
	*info = -15;
    } else if (*ldz < 1 || *wantz && *ldz < *n) {
	*info = -17;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03YA", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    } else {
	i1 = *ilo;
	i2 = *ihi;
    }

/*     Apply a zero-shifted QR step. */

    i__1 = *pos - 1;
    for (j = *ilo; j <= i__1; ++j) {
	temp = a[j + j * a_dim1];
	dlartg_(&temp, &a[j + 1 + j * a_dim1], &cs, &sn, &a[j + j * a_dim1]);
	a[j + 1 + j * a_dim1] = 0.;
	i__2 = i2 - j;
	drot_(&i__2, &a[j + (j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * 
		a_dim1], lda, &cs, &sn);
/* Computing MIN */
	i__3 = j, i__4 = *pos - 2;
	i__2 = min(i__3,i__4) - i1 + 2;
	drot_(&i__2, &b[i1 + j * b_dim1], &c__1, &b[i1 + (j + 1) * b_dim1], &
		c__1, &cs, &sn);
	if (*wantq) {
	    drot_(&nq, &q[*iloq + j * q_dim1], &c__1, &q[*iloq + (j + 1) * 
		    q_dim1], &c__1, &cs, &sn);
	}
/* L10: */
    }
    i__1 = *pos - 2;
    for (j = *ilo; j <= i__1; ++j) {
	temp = b[j + j * b_dim1];
	dlartg_(&temp, &b[j + 1 + j * b_dim1], &cs, &sn, &b[j + j * b_dim1]);
	b[j + 1 + j * b_dim1] = 0.;
	i__2 = i2 - j;
	drot_(&i__2, &b[j + (j + 1) * b_dim1], ldb, &b[j + 1 + (j + 1) * 
		b_dim1], ldb, &cs, &sn);
	i__2 = j - i1 + 2;
	drot_(&i__2, &a[i1 + j * a_dim1], &c__1, &a[i1 + (j + 1) * a_dim1], &
		c__1, &cs, &sn);
	if (*wantz) {
	    drot_(&nq, &z__[*iloq + j * z_dim1], &c__1, &z__[*iloq + (j + 1) *
		     z_dim1], &c__1, &cs, &sn);
	}
/* L20: */
    }

/*     Apply a zero-shifted RQ step. */

    i__1 = *pos + 1;
    for (j = *ihi; j >= i__1; --j) {
	temp = a[j + j * a_dim1];
	dlartg_(&temp, &a[j + (j - 1) * a_dim1], &cs, &sn, &a[j + j * a_dim1])
		;
	a[j + (j - 1) * a_dim1] = 0.;
	sn = -sn;
	i__2 = j - i1;
	drot_(&i__2, &a[i1 + (j - 1) * a_dim1], &c__1, &a[i1 + j * a_dim1], &
		c__1, &cs, &sn);
/* Computing MAX */
	i__3 = j - 1, i__4 = *pos + 1;
	i__2 = i2 - max(i__3,i__4) + 1;
/* Computing MAX */
	i__5 = j - 1, i__6 = *pos + 1;
/* Computing MAX */
	i__7 = j - 1, i__8 = *pos + 1;
	drot_(&i__2, &b[j - 1 + max(i__5,i__6) * b_dim1], ldb, &b[j + max(
		i__7,i__8) * b_dim1], ldb, &cs, &sn);
	if (*wantz) {
	    drot_(&nq, &z__[*iloq + (j - 1) * z_dim1], &c__1, &z__[*iloq + j *
		     z_dim1], &c__1, &cs, &sn);
	}
/* L30: */
    }
    i__1 = *pos + 2;
    for (j = *ihi; j >= i__1; --j) {
	temp = b[j + j * b_dim1];
	dlartg_(&temp, &b[j + (j - 1) * b_dim1], &cs, &sn, &b[j + j * b_dim1])
		;
	b[j + (j - 1) * b_dim1] = 0.;
	sn = -sn;
	i__2 = j - i1;
	drot_(&i__2, &b[i1 + (j - 1) * b_dim1], &c__1, &b[i1 + j * b_dim1], &
		c__1, &cs, &sn);
	i__2 = i2 - j + 2;
	drot_(&i__2, &a[j - 1 + (j - 1) * a_dim1], lda, &a[j + (j - 1) * 
		a_dim1], lda, &cs, &sn);
	if (*wantq) {
	    drot_(&nq, &q[*iloq + (j - 1) * q_dim1], &c__1, &q[*iloq + j * 
		    q_dim1], &c__1, &cs, &sn);
	}
/* L40: */
    }
    return 0;
/* *** Last line of MB03YA *** */
} /* mb03ya_ */

