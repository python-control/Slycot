/* MB04WU.f -- translated by f2c (version 20100827).
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

static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb04wu_(char *tranq1, char *tranq2, integer *m, integer *
	n, integer *k, doublereal *q1, integer *ldq1, doublereal *q2, integer 
	*ldq2, doublereal *cs, doublereal *tau, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen tranq1_len, ftnlen tranq2_len)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal nu;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltrq1, ltrq2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To generate a matrix Q with orthogonal columns (spanning an */
/*     isotropic subspace), which is defined as the first n columns */
/*     of a product of symplectic reflectors and Givens rotators, */

/*         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*             diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                               .... */
/*             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     The matrix Q is returned in terms of its first 2*M rows */

/*                      [  op( Q1 )   op( Q2 ) ] */
/*                  Q = [                      ]. */
/*                      [ -op( Q2 )   op( Q1 ) ] */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANQ1  CHARACTER*1 */
/*             Specifies the form of op( Q1 ) as follows: */
/*             = 'N':  op( Q1 ) = Q1; */
/*             = 'T':  op( Q1 ) = Q1'; */
/*             = 'C':  op( Q1 ) = Q1'. */

/*     TRANQ2  CHARACTER*1 */
/*             Specifies the form of op( Q2 ) as follows: */
/*             = 'N':  op( Q2 ) = Q2; */
/*             = 'T':  op( Q2 ) = Q2'; */
/*             = 'C':  op( Q2 ) = Q2'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices Q1 and Q2. M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices Q1 and Q2. */
/*             M >= N >= 0. */

/*     K       (input) INTEGER */
/*             The number of symplectic Givens rotators whose product */
/*             partly defines the matrix Q. N >= K >= 0. */

/*     Q1      (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDQ1,N) if TRANQ1 = 'N', */
/*                     (LDQ1,M) if TRANQ1 = 'T' or TRANQ1 = 'C' */
/*             On entry with TRANQ1 = 'N', the leading M-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector F(i). */
/*             On entry with TRANQ1 = 'T' or TRANQ1 = 'C', the leading */
/*             K-by-M part of this array must contain in its i-th row */
/*             the vector which defines the elementary reflector F(i). */
/*             On exit with TRANQ1 = 'N', the leading M-by-N part of this */
/*             array contains the matrix Q1. */
/*             On exit with TRANQ1 = 'T' or TRANQ1 = 'C', the leading */
/*             N-by-M part of this array contains the matrix Q1'. */

/*     LDQ1    INTEGER */
/*             The leading dimension of the array Q1. */
/*             LDQ1 >= MAX(1,M),  if TRANQ1 = 'N'; */
/*             LDQ1 >= MAX(1,N),  if TRANQ1 = 'T' or TRANQ1 = 'C'. */

/*     Q2      (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDQ2,N) if TRANQ2 = 'N', */
/*                     (LDQ2,M) if TRANQ2 = 'T' or TRANQ2 = 'C' */
/*             On entry with TRANQ2 = 'N', the leading M-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector H(i) and, on the */
/*             diagonal, the scalar factor of H(i). */
/*             On entry with TRANQ2 = 'T' or TRANQ2 = 'C', the leading */
/*             K-by-M part of this array must contain in its i-th row the */
/*             vector which defines the elementary reflector H(i) and, on */
/*             the diagonal, the scalar factor of H(i). */
/*             On exit with TRANQ2 = 'N', the leading M-by-N part of this */
/*             array contains the matrix Q2. */
/*             On exit with TRANQ2 = 'T' or TRANQ2 = 'C', the leading */
/*             N-by-M part of this array contains the matrix Q2'. */

/*     LDQ2    INTEGER */
/*             The leading dimension of the array Q2. */
/*             LDQ2 >= MAX(1,M),  if TRANQ2 = 'N'; */
/*             LDQ2 >= MAX(1,N),  if TRANQ2 = 'T' or TRANQ2 = 'C'. */

/*     CS      (input) DOUBLE PRECISION array, dimension (2*K) */
/*             On entry, the first 2*K elements of this array must */
/*             contain the cosines and sines of the symplectic Givens */
/*             rotators G(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             On entry, the first K elements of this array must */
/*             contain the scalar factors of the elementary reflectors */
/*             F(i). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -13,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,M+N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Bunse-Gerstner, A. */
/*         Matrix factorizations for symplectic QR-like methods. */
/*         Linear Algebra Appl., 83, pp. 49-77, 1986. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSQ). */

/*     KEYWORDS */

/*     Elementary matrix operations, orthogonal symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    q1_dim1 = *ldq1;
    q1_offset = 1 + q1_dim1;
    q1 -= q1_offset;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --cs;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    ltrq1 = lsame_(tranq1, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq1, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrq2 = lsame_(tranq2, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq2, "C", (
	    ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (ltrq1 || lsame_(tranq1, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (ltrq2 || lsame_(tranq2, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0 || *n > *m) {
	*info = -4;
    } else if (*k < 0 || *k > *n) {
	*info = -5;
    } else if (ltrq1 && *ldq1 < max(1,*n) || ! ltrq1 && *ldq1 < max(1,*m)) {
	*info = -7;
    } else if (ltrq2 && *ldq2 < max(1,*n) || ! ltrq2 && *ldq2 < max(1,*m)) {
	*info = -9;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m + *n;
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
	    i__1 = 1, i__2 = *m + *n;
	    dwork[1] = (doublereal) max(i__1,i__2);
	    *info = -13;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04WU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Initialize columns K+1:N to columns of the unit matrix. */

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q1[i__ + j * q1_dim1] = 0.;
/* L10: */
	}
	q1[j + j * q1_dim1] = 1.;
/* L20: */
    }
    i__1 = *n - *k;
    dlaset_("All", m, &i__1, &c_b12, &c_b12, &q2[(*k + 1) * q2_dim1 + 1], 
	    ldq2, (ftnlen)3);

    if (ltrq1 && ltrq2) {
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I+1:N,I:M) and Q2(I+1:N,I:M) from the */
/*           right. */

	    i__1 = *m - i__ + 1;
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[1], &c__1);
	    if (i__ < *n) {
		q1[i__ + i__ * q1_dim1] = 1.;
		i__1 = *n - i__;
		i__2 = *m - i__ + 1;
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q1[i__ + 1 + i__ * q1_dim1], ldq1, &dwork[
			*m + 1], (ftnlen)5);
		i__1 = *n - i__;
		i__2 = *m - i__ + 1;
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q2[i__ + 1 + i__ * q2_dim1], ldq2, &dwork[
			*m + 1], (ftnlen)5);
	    }
	    if (i__ < *m) {
		i__1 = *m - i__;
		d__1 = -tau[i__];
		dscal_(&i__1, &d__1, &q1[i__ + (i__ + 1) * q1_dim1], ldq1);
	    }
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(I,1:I-1) and Q2(I,1:M) to zero. */

	    i__1 = i__ - 1;
	    for (j = 1; j <= i__1; ++j) {
		q1[i__ + j * q1_dim1] = 0.;
/* L30: */
	    }
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		q2[i__ + j * q2_dim1] = 0.;
/* L40: */
	    }

/*           Apply G(I) to Q1(I:N,I) and Q2(I:N,I) from the right. */

	    i__1 = *n - i__ + 1;
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], &c__1, &q2[i__ + i__ * 
		    q2_dim1], &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:N,I:M) and Q2(I:N,I:M) from the right. */

	    nu = dwork[1];
	    dwork[1] = 1.;
	    i__1 = *n - i__ + 1;
	    i__2 = *m - i__ + 1;
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + 
		    i__ * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)5);
	    i__1 = *n - i__ + 1;
	    i__2 = *m - i__ + 1;
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + 
		    i__ * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)5);
/* L50: */
	}
    } else if (ltrq1) {
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I+1:N,I:M) from the right and to */
/*           Q2(I:M,I+1:N) from the left. */

	    i__1 = *m - i__ + 1;
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], &c__1, &dwork[1], &c__1);
	    if (i__ < *n) {
		q1[i__ + i__ * q1_dim1] = 1.;
		i__1 = *n - i__;
		i__2 = *m - i__ + 1;
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q1[i__ + 1 + i__ * q1_dim1], ldq1, &dwork[
			*m + 1], (ftnlen)5);
		i__1 = *m - i__ + 1;
		i__2 = *n - i__;
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, &
			tau[i__], &q2[i__ + (i__ + 1) * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)4);
	    }
	    if (i__ < *m) {
		i__1 = *m - i__;
		d__1 = -tau[i__];
		dscal_(&i__1, &d__1, &q1[i__ + (i__ + 1) * q1_dim1], ldq1);
	    }
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(I,1:I-1) and Q2(1:M,I) to zero. */

	    i__1 = i__ - 1;
	    for (j = 1; j <= i__1; ++j) {
		q1[i__ + j * q1_dim1] = 0.;
/* L60: */
	    }
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		q2[j + i__ * q2_dim1] = 0.;
/* L70: */
	    }

/*           Apply G(I) to Q1(I:N,I) from the right and to Q2(I,I:N) */
/*           from the left. */

	    i__1 = *n - i__ + 1;
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], &c__1, &q2[i__ + i__ * 
		    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:N,I:M) from the right and to Q2(I:M,I:N) */
/*           from the left. */

	    nu = dwork[1];
	    dwork[1] = 1.;
	    i__1 = *n - i__ + 1;
	    i__2 = *m - i__ + 1;
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + 
		    i__ * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)5);
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__ + 1;
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + i__ 
		    * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)4);
/* L80: */
	}
    } else if (ltrq2) {
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I:M,I+1:N) from the left and to */
/*           Q2(I+1:N,I:M) from the right. */

	    i__1 = *m - i__ + 1;
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[1], &c__1);
	    if (i__ < *n) {
		q1[i__ + i__ * q1_dim1] = 1.;
		i__1 = *m - i__ + 1;
		i__2 = *n - i__;
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q1[i__ + (i__ + 1) * q1_dim1], ldq1, &
			dwork[*m + 1], (ftnlen)4);
		i__1 = *n - i__;
		i__2 = *m - i__ + 1;
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1,
			 &tau[i__], &q2[i__ + 1 + i__ * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)5);
	    }
	    if (i__ < *m) {
		i__1 = *m - i__;
		d__1 = -tau[i__];
		dscal_(&i__1, &d__1, &q1[i__ + 1 + i__ * q1_dim1], &c__1);
	    }
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(1:I-1,I) and Q2(I,1:M) to zero. */

	    i__1 = i__ - 1;
	    for (j = 1; j <= i__1; ++j) {
		q1[j + i__ * q1_dim1] = 0.;
/* L90: */
	    }
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		q2[i__ + j * q2_dim1] = 0.;
/* L100: */
	    }

/*           Apply G(I) to Q1(I,I:N) from the left and to Q2(I:N,I) */
/*           from the right. */

	    i__1 = *n - i__ + 1;
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
		    q2_dim1], &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:M,I:N) from the left and to Q2(I:N,I:M) */
/*           from the left. */

	    nu = dwork[1];
	    dwork[1] = 1.;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__ + 1;
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + i__ 
		    * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)4);
	    i__1 = *n - i__ + 1;
	    i__2 = *m - i__ + 1;
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + 
		    i__ * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)5);
/* L110: */
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I:M,I+1:N) and Q2(I:M,I+1:N) from the left. */

	    i__1 = *m - i__ + 1;
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], &c__1, &dwork[1], &c__1);
	    if (i__ < *n) {
		q1[i__ + i__ * q1_dim1] = 1.;
		i__1 = *m - i__ + 1;
		i__2 = *n - i__;
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q1[i__ + (i__ + 1) * q1_dim1], ldq1, &
			dwork[*m + 1], (ftnlen)4);
		i__1 = *m - i__ + 1;
		i__2 = *n - i__;
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q2[i__ + (i__ + 1) * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)4);
	    }
	    if (i__ < *m) {
		i__1 = *m - i__;
		d__1 = -tau[i__];
		dscal_(&i__1, &d__1, &q1[i__ + 1 + i__ * q1_dim1], &c__1);
	    }
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(1:I-1,I) and Q2(1:M,I) to zero. */

	    i__1 = i__ - 1;
	    for (j = 1; j <= i__1; ++j) {
		q1[j + i__ * q1_dim1] = 0.;
/* L120: */
	    }
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		q2[j + i__ * q2_dim1] = 0.;
/* L130: */
	    }

/*           Apply G(I) to Q1(I,I:N) and Q2(I,I:N) from the left. */

	    i__1 = *n - i__ + 1;
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
		    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:M,I:N) and Q2(I:M,I:N) from the left. */

	    nu = dwork[1];
	    dwork[1] = 1.;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__ + 1;
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + i__ 
		    * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)4);
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__ + 1;
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + i__ 
		    * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)4);
/* L140: */
	}
    }
/* Computing MAX */
    i__1 = 1, i__2 = *m + *n;
    dwork[1] = (doublereal) max(i__1,i__2);
/* *** Last line of MB04WU *** */
    return 0;
} /* mb04wu_ */

