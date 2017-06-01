/* MB01QD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb01qd_(char *type__, integer *m, integer *n, integer *
	kl, integer *ku, doublereal *cfrom, doublereal *cto, integer *nbl, 
	integer *nrows, doublereal *a, integer *lda, integer *info, ftnlen 
	type_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, k1, k2, k3, k4;
    static doublereal mul, cto1;
    static logical done;
    static integer ifin, jfin;
    static doublereal ctoc;
    static integer jini;
    static logical noblc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itype;
    static doublereal cfrom1;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cfromc, bignum, smlnum;


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

/*     To multiply the M by N real matrix A by the real scalar CTO/CFROM. */
/*     This is done without over/underflow as long as the final result */
/*     CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that */
/*     A may be full, (block) upper triangular, (block) lower triangular, */
/*     (block) upper Hessenberg, or banded. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             TYPE indices the storage type of the input matrix. */
/*             = 'G':  A is a full matrix. */
/*             = 'L':  A is a (block) lower triangular matrix. */
/*             = 'U':  A is a (block) upper triangular matrix. */
/*             = 'H':  A is a (block) upper Hessenberg matrix. */
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
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     KL      (input) INTEGER */
/*             The lower bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     KU      (input) INTEGER */
/*             The upper bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     CFROM   (input) DOUBLE PRECISION */
/*     CTO     (input) DOUBLE PRECISION */
/*             The matrix A is multiplied by CTO/CFROM. A(I,J) is */
/*             computed without over/underflow if the final result */
/*             CTO*A(I,J)/CFROM can be represented without over/ */
/*             underflow.  CFROM must be nonzero. */

/*     NBL     (input) INTEGER */
/*             The number of diagonal blocks of the matrix A, if it has a */
/*             block structure.  To specify that matrix A has no block */
/*             structure, set NBL = 0.  NBL >= 0. */

/*     NROWS   (input) INTEGER array, dimension max(1,NBL) */
/*             NROWS(i) contains the number of rows and columns of the */
/*             i-th diagonal block of matrix A.  The sum of the values */
/*             NROWS(i),  for  i = 1: NBL,  should be equal to min(M,N). */
/*             The array  NROWS  is not referenced if NBL = 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The matrix to be multiplied by CTO/CFROM.  See TYPE for */
/*             the storage type. */

/*     LDA     (input) INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             Not used in this implementation. */

/*     METHOD */

/*     Matrix A is multiplied by the real scalar CTO/CFROM, taking into */
/*     account the specified storage mode of the matrix. */
/*     MB01QD is a version of the LAPACK routine DLASCL, modified for */
/*     dealing with block triangular, or block Hessenberg matrices. */
/*     For efficiency, no tests of the input scalar parameters are */
/*     performed. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --nrows;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
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
    } else {
	itype = 6;
    }

/*     Quick return if possible. */

    if (min(*m,*n) == 0) {
	return 0;
    }

/*     Get machine parameters. */

    smlnum = dlamch_("S", (ftnlen)1);
    bignum = 1. / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (abs(cto1) > abs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    noblc = *nbl == 0;

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

	if (noblc) {

/*           Lower triangular matrix */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = j; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] *= mul;
/* L40: */
		}
/* L50: */
	    }

	} else {

/*           Block lower triangular matrix */

	    jfin = 0;
	    i__1 = *nbl;
	    for (k = 1; k <= i__1; ++k) {
		jini = jfin + 1;
		jfin += nrows[k];
		i__2 = jfin;
		for (j = jini; j <= i__2; ++j) {
		    i__3 = *m;
		    for (i__ = jini; i__ <= i__3; ++i__) {
			a[i__ + j * a_dim1] *= mul;
/* L60: */
		    }
/* L70: */
		}
/* L80: */
	    }
	}

    } else if (itype == 2) {

	if (noblc) {

/*           Upper triangular matrix */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = min(j,*m);
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] *= mul;
/* L90: */
		}
/* L100: */
	    }

	} else {

/*           Block upper triangular matrix */

	    jfin = 0;
	    i__1 = *nbl;
	    for (k = 1; k <= i__1; ++k) {
		jini = jfin + 1;
		jfin += nrows[k];
		if (k == *nbl) {
		    jfin = *n;
		}
		i__2 = jfin;
		for (j = jini; j <= i__2; ++j) {
		    i__3 = min(jfin,*m);
		    for (i__ = 1; i__ <= i__3; ++i__) {
			a[i__ + j * a_dim1] *= mul;
/* L110: */
		    }
/* L120: */
		}
/* L130: */
	    }
	}

    } else if (itype == 3) {

	if (noblc) {

/*           Upper Hessenberg matrix */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__3 = j + 1;
		i__2 = min(i__3,*m);
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] *= mul;
/* L140: */
		}
/* L150: */
	    }

	} else {

/*           Block upper Hessenberg matrix */

	    jfin = 0;
	    i__1 = *nbl;
	    for (k = 1; k <= i__1; ++k) {
		jini = jfin + 1;
		jfin += nrows[k];

		if (k == *nbl) {
		    jfin = *n;
		    ifin = *n;
		} else {
		    ifin = jfin + nrows[k + 1];
		}

		i__2 = jfin;
		for (j = jini; j <= i__2; ++j) {
		    i__3 = min(ifin,*m);
		    for (i__ = 1; i__ <= i__3; ++i__) {
			a[i__ + j * a_dim1] *= mul;
/* L160: */
		    }
/* L170: */
		}
/* L180: */
	    }
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= mul;
/* L190: */
	    }
/* L200: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
		a[i__ + j * a_dim1] *= mul;
/* L210: */
	    }
/* L220: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = min(i__4,i__5);
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] *= mul;
/* L230: */
	    }
/* L240: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;
/* *** Last line of MB01QD *** */
} /* mb01qd_ */

