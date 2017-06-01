/* MB01UX.f -- translated by f2c (version 20100827).
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

static doublereal c_b11 = 0.;
static integer c__1 = 1;
static doublereal c_b31 = 1.;
static integer c__0 = 0;

/* Subroutine */ int mb01ux_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *alpha, doublereal *t, integer *ldt, 
	doublereal *a, integer *lda, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen side_len, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, pdw;
    static logical lup;
    static integer noff, xdif, ierr;
    static doublereal temp;
    static integer psav;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char atran[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical ltran;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrmv_(char *, char *, char *, integer *, doublereal *
	    , integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dlaset_(char *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer wrkmin, wrkopt;


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

/*     To compute one of the matrix products */

/*       A : = alpha*op( T ) * A, or A : = alpha*A * op( T ), */

/*     where alpha is a scalar, A is an m-by-n matrix, T is a quasi- */
/*     triangular matrix, and op( T ) is one of */

/*        op( T ) = T   or   op( T ) = T',  the transpose of T. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the upper quasi-triangular matrix H */
/*             appears on the left or right in the matrix product as */
/*             follows: */
/*             = 'L':  A := alpha*op( T ) * A; */
/*             = 'R':  A := alpha*A * op( T ). */

/*     UPLO    CHARACTER*1. */
/*             Specifies whether the matrix T is an upper or lower */
/*             quasi-triangular matrix as follows: */
/*             = 'U':  T is an upper quasi-triangular matrix; */
/*             = 'L':  T is a lower quasi-triangular matrix. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( T ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( T ) = T; */
/*             = 'T':  op( T ) = T'; */
/*             = 'C':  op( T ) = T'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then T is not */
/*             referenced and A need not be set before entry. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,k) */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */
/*             On entry with UPLO = 'U', the leading k-by-k upper */
/*             Hessenberg part of this array must contain the upper */
/*             quasi-triangular matrix T. The elements below the */
/*             subdiagonal are not referenced. */
/*             On entry with UPLO = 'L', the leading k-by-k lower */
/*             Hessenberg part of this array must contain the lower */
/*             quasi-triangular matrix T. The elements above the */
/*             supdiagonal are not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= max(1,k), */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the computed product. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 and ALPHA<>0,  DWORK(1)  returns the */
/*             optimal value of LDWORK. */
/*             On exit, if  INFO = -12,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             This array is not referenced when alpha = 0. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= 1,       if alpha =  0 or MIN(M,N) = 0; */
/*             LDWORK >= 2*(M-1), if SIDE  = 'L'; */
/*             LDWORK >= 2*(N-1), if SIDE  = 'R'. */
/*             For maximal efficiency LDWORK should be at least */
/*             NOFF*N + M - 1,    if SIDE  = 'L'; */
/*             NOFF*M + N - 1,    if SIDE  = 'R'; */
/*             where NOFF is the number of nonzero elements on the */
/*             subdiagonal (if UPLO = 'U') or supdiagonal (if UPLO = 'L') */
/*             of T. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The technique used in this routine is similiar to the technique */
/*     used in the SLICOT [1] subroutine MB01UW developed by Vasile Sima. */
/*     The required matrix product is computed in two steps. In the first */
/*     step, the triangle of T specified by UPLO is used; in the second */
/*     step, the contribution of the sub-/supdiagonal is added. If the */
/*     workspace can accommodate parts of A, a fast BLAS 3 DTRMM */
/*     operation is used in the first step. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., Sima, V., Van Huffel, S., and */
/*         Varga, A. */
/*         SLICOT - A subroutine library in systems and control theory. */
/*         In: Applied and computational control, signals, and circuits, */
/*         Vol. 1, pp. 499-539, Birkhauser, Boston, 1999. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DTRQML). */

/*     KEYWORDS */

/*     Elementary matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and test the input scalar arguments. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    lup = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ltran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);
    if (lside) {
	k = *m;
    } else {
	k = *n;
    }
    wrkmin = k - 1 << 1;

    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lup && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! ltran && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*ldt < max(1,k)) {
	*info = -8;
    } else if (*lda < max(1,*m)) {
	*info = -10;
    } else if (*ldwork < 0 || *alpha != 0. && min(*m,*n) > 0 && *ldwork < 
	    wrkmin) {
	dwork[1] = (doublereal) wrkmin;
	*info = -12;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB01UX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return, if possible. */

    if (min(*m,*n) == 0) {
	return 0;
    }

    if (*alpha == 0.) {

/*        Set A to zero and return. */

	dlaset_("Full", m, n, &c_b11, &c_b11, &a[a_offset], lda, (ftnlen)4);
	return 0;
    }

/*     Save and count off-diagonal entries of T. */

    if (lup) {
	i__1 = k - 1;
	i__2 = *ldt + 1;
	dcopy_(&i__1, &t[t_dim1 + 2], &i__2, &dwork[1], &c__1);
    } else {
	i__1 = k - 1;
	i__2 = *ldt + 1;
	dcopy_(&i__1, &t[(t_dim1 << 1) + 1], &i__2, &dwork[1], &c__1);
    }
    noff = 0;
    i__1 = k - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dwork[i__] != 0.) {
	    ++noff;
	}
/* L5: */
    }

/*     Compute optimal workspace. */

    if (lside) {
	wrkopt = noff * *n + *m - 1;
    } else {
	wrkopt = noff * *m + *n - 1;
    }
    psav = k;
    if (! ltran) {
	xdif = 0;
    } else {
	xdif = 1;
    }
    if (! lup) {
	xdif = 1 - xdif;
    }
    if (! lside) {
	xdif = 1 - xdif;
    }

    if (*ldwork >= wrkopt) {

/*        Enough workspace for a fast BLAS 3 calculation. */
/*        Save relevant parts of A in the workspace and compute one of */
/*        the matrix products */
/*          A : = alpha*op( triu( T ) ) * A, or */
/*          A : = alpha*A * op( triu( T ) ), */
/*        involving the upper/lower triangle of T. */

	pdw = psav;
	if (lside) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (dwork[i__] != 0.) {
			dwork[pdw] = a[i__ + xdif + j * a_dim1];
			++pdw;
		    }
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		if (dwork[j] != 0.) {
		    dcopy_(m, &a[(j + xdif) * a_dim1 + 1], &c__1, &dwork[pdw],
			     &c__1);
		    pdw += *m;
		}
/* L30: */
	    }
	}
	dtrmm_(side, uplo, trans, "Non-unit", m, n, alpha, &t[t_offset], ldt, 
		&a[a_offset], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)8)
		;

/*        Add the contribution of the offdiagonal of T. */

	pdw = psav;
	xdif = 1 - xdif;
	if (lside) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = dwork[i__];
		    if (temp != 0.) {
			a[i__ + xdif + j * a_dim1] += *alpha * temp * dwork[
				pdw];
			++pdw;
		    }
/* L40: */
		}
/* L50: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		temp = dwork[j] * *alpha;
		if (temp != 0.) {
		    daxpy_(m, &temp, &dwork[pdw], &c__1, &a[(j + xdif) * 
			    a_dim1 + 1], &c__1);
		    pdw += *m;
		}
/* L60: */
	    }
	}
    } else {

/*        Use a BLAS 2 calculation. */

	if (lside) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

/*              Compute the contribution of the offdiagonal of T to */
/*              the j-th column of the product. */

		i__2 = *m - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[psav + i__ - 1] = dwork[i__] * a[i__ + xdif + j * 
			    a_dim1];
/* L70: */
		}

/*              Multiply the triangle of T by the j-th column of A, */
/*              and add to the above result. */

		dtrmv_(uplo, trans, "Non-unit", m, &t[t_offset], ldt, &a[j * 
			a_dim1 + 1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
		i__2 = *m - 1;
		daxpy_(&i__2, &c_b31, &dwork[psav], &c__1, &a[2 - xdif + j * 
			a_dim1], &c__1);
/* L80: */
	    }
	} else {
	    if (ltran) {
		*(unsigned char *)atran = 'N';
	    } else {
		*(unsigned char *)atran = 'T';
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Compute the contribution of the offdiagonal of T to */
/*              the i-th row of the product. */

		i__2 = *n - 1;
		for (j = 1; j <= i__2; ++j) {
		    dwork[psav + j - 1] = a[i__ + (j + xdif) * a_dim1] * 
			    dwork[j];
/* L90: */
		}

/*              Multiply the i-th row of A by the triangle of T, */
/*              and add to the above result. */

		dtrmv_(uplo, atran, "Non-unit", n, &t[t_offset], ldt, &a[i__ 
			+ a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)8);
		i__2 = *n - 1;
		daxpy_(&i__2, &c_b31, &dwork[psav], &c__1, &a[i__ + (2 - xdif)
			 * a_dim1], lda);
/* L100: */
	    }
	}

/*        Scale the result by alpha. */

	if (*alpha != 1.) {
	    dlascl_("General", &c__0, &c__0, &c_b31, alpha, m, n, &a[a_offset]
		    , lda, &ierr, (ftnlen)7);
	}
    }
    dwork[1] = (doublereal) max(wrkmin,wrkopt);
    return 0;
/* *** Last line of MB01UX *** */
} /* mb01ux_ */

