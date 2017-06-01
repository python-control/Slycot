/* MB01PD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb01pd_(char *scun, char *type__, integer *m, integer *n,
	 integer *kl, integer *ku, doublereal *anrm, integer *nbl, integer *
	nrows, doublereal *a, integer *lda, integer *info, ftnlen scun_len, 
	ftnlen type_len)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, mn, isum;
    extern /* Subroutine */ int mb01qd_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itype;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical lscale;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


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

/*     To scale a matrix or undo scaling.  Scaling is performed, if */
/*     necessary, so that the matrix norm will be in a safe range of */
/*     representable numbers. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SCUN    CHARACTER*1 */
/*             SCUN indicates the operation to be performed. */
/*             = 'S':  scale the matrix. */
/*             = 'U':  undo scaling of the matrix. */

/*     TYPE    CHARACTER*1 */
/*             TYPE indicates the storage type of the input matrix. */
/*             = 'G':  A is a full matrix. */
/*             = 'L':  A is a (block) lower triangular matrix. */
/*             = 'U':  A is an (block) upper triangular matrix. */
/*             = 'H':  A is an (block) upper Hessenberg matrix. */
/*             = 'B':  A is a symmetric band matrix with lower bandwidth */
/*                     KL and upper bandwidth KU and with the only the */
/*                     lower half stored. */
/*             = 'Q':  A is a symmetric band matrix with lower bandwidth */
/*                     KL and upper bandwidth KU and with the only the */
/*                     upper half stored. */
/*             = 'Z':  A is a band matrix with lower bandwidth KL and */
/*                     upper bandwidth KU. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A. M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A. N >= 0. */

/*     KL      (input) INTEGER */
/*             The lower bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     KU      (input) INTEGER */
/*             The upper bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     ANRM    (input) DOUBLE PRECISION */
/*             The norm of the initial matrix A.  ANRM >= 0. */
/*             When  ANRM = 0  then an immediate return is effected. */
/*             ANRM should be preserved between the call of the routine */
/*             with SCUN = 'S' and the corresponding one with SCUN = 'U'. */

/*     NBL     (input) INTEGER */
/*             The number of diagonal blocks of the matrix A, if it has a */
/*             block structure.  To specify that matrix A has no block */
/*             structure, set NBL = 0.  NBL >= 0. */

/*     NROWS   (input) INTEGER array, dimension max(1,NBL) */
/*             NROWS(i) contains the number of rows and columns of the */
/*             i-th diagonal block of matrix A.  The sum of the values */
/*             NROWS(i),  for  i = 1: NBL,  should be equal to min(M,N). */
/*             The elements of the array  NROWS  are not referenced if */
/*             NBL = 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M by N part of this array must */
/*             contain the matrix to be scaled/unscaled. */
/*             On exit, the leading M by N part of A will contain */
/*             the modified matrix. */
/*             The storage mode of A is specified by TYPE. */

/*     LDA     (input) INTEGER */
/*             The leading dimension of the array A.  LDA  >= max(1,M). */

/*     Error Indicator */

/*     INFO    (output) INTEGER */
/*             = 0:  successful exit */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Denote by ANRM the norm of the matrix, and by SMLNUM and BIGNUM, */
/*     two positive numbers near the smallest and largest safely */
/*     representable numbers, respectively.  The matrix is scaled, if */
/*     needed, such that the norm of the result is in the range */
/*     [SMLNUM, BIGNUM].  The scaling factor is represented as a ratio */
/*     of two numbers, one of them being ANRM, and the other one either */
/*     SMLNUM or BIGNUM, depending on ANRM being less than SMLNUM or */
/*     larger than BIGNUM, respectively.  For undoing the scaling, the */
/*     norm is again compared with SMLNUM or BIGNUM, and the reciprocal */
/*     of the previous scaling factor is used. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*     REVISIONS */

/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save statement .. */
/*     .. Data statements .. */
    /* Parameter adjustments */
    --nrows;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    *info = 0;
    lscale = lsame_(scun, "S", (ftnlen)1, (ftnlen)1);
    if (lsame_(type__, "G", (ftnlen)1, (ftnlen)1)) {
	itype = 0;
    } else if (lsame_(type__, "L", (ftnlen)1, (ftnlen)1)) {
	itype = 1;
    } else if (lsame_(type__, "U", (ftnlen)1, (ftnlen)1)) {
	itype = 2;
    } else if (lsame_(type__, "H", (ftnlen)1, (ftnlen)1)) {
	itype = 3;
    } else if (lsame_(type__, "B", (ftnlen)1, (ftnlen)1)) {
	itype = 4;
    } else if (lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
	itype = 5;
    } else if (lsame_(type__, "Z", (ftnlen)1, (ftnlen)1)) {
	itype = 6;
    } else {
	itype = -1;
    }

    mn = min(*m,*n);

    isum = 0;
    if (*nbl > 0) {
	i__1 = *nbl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    isum += nrows[i__];
/* L10: */
	}
    }

    if (! lscale && ! lsame_(scun, "U", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (itype == -1) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0 || (itype == 4 || itype == 5) && *n != *m) {
	*info = -4;
    } else if (*anrm < 0.) {
	*info = -7;
    } else if (*nbl < 0) {
	*info = -8;
    } else if (*nbl > 0 && isum != mn) {
	*info = -9;
    } else if (itype <= 3 && *lda < max(1,*m)) {
	*info = -11;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > max(i__1,0)) {
	    *info = -5;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
		*info = -6;
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
		*info = -11;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB01PD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (mn == 0 || *anrm == 0.) {
	return 0;
    }

    if (first) {

/*        Get machine parameters. */

	smlnum = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
	bignum = 1. / smlnum;
	dlabad_(&smlnum, &bignum);
	first = FALSE_;
    }

    if (lscale) {

/*        Scale A, if its norm is outside range [SMLNUM,BIGNUM]. */

	if (*anrm < smlnum) {

/*           Scale matrix norm up to SMLNUM. */

	    mb01qd_(type__, m, n, kl, ku, anrm, &smlnum, nbl, &nrows[1], &a[
		    a_offset], lda, info, (ftnlen)1);
	} else if (*anrm > bignum) {

/*           Scale matrix norm down to BIGNUM. */

	    mb01qd_(type__, m, n, kl, ku, anrm, &bignum, nbl, &nrows[1], &a[
		    a_offset], lda, info, (ftnlen)1);
	}

    } else {

/*        Undo scaling. */

	if (*anrm < smlnum) {
	    mb01qd_(type__, m, n, kl, ku, &smlnum, anrm, nbl, &nrows[1], &a[
		    a_offset], lda, info, (ftnlen)1);
	} else if (*anrm > bignum) {
	    mb01qd_(type__, m, n, kl, ku, &bignum, anrm, nbl, &nrows[1], &a[
		    a_offset], lda, info, (ftnlen)1);
	}
    }

    return 0;
/* *** Last line of MB01PD *** */
} /* mb01pd_ */

