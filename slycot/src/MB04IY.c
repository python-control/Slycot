/* MB04IY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04iy_(char *side, char *trans, integer *n, integer *m, 
	integer *k, integer *p, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal aii;
    static logical left, tran;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dormqr_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal wrkopt;


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

/*     To overwrite the real n-by-m matrix  C  with  Q' * C,  Q * C, */
/*     C * Q',  or  C * Q,  according to the following table */

/*                     SIDE = 'L'     SIDE = 'R' */
/*     TRANS = 'N':      Q * C          C * Q */
/*     TRANS = 'T':      Q'* C          C * Q' */

/*     where  Q  is a real orthogonal matrix defined as the product of */
/*     k elementary reflectors */

/*        Q = H(1) H(2) . . . H(k) */

/*     as returned by SLICOT Library routine MB04ID.  Q  is of order n */
/*     if  SIDE = 'L'  and of order m if  SIDE = 'R'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specify if  Q  or  Q'  is applied from the left or right, */
/*             as follows: */
/*             = 'L':  apply  Q  or  Q'  from the left; */
/*             = 'R':  apply  Q  or  Q'  from the right. */

/*     TRANS   CHARACTER*1 */
/*             Specify if  Q  or  Q'  is to be applied, as follows: */
/*             = 'N':  apply  Q   (No transpose); */
/*             = 'T':  apply  Q'  (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix C.  M >= 0. */

/*     K       (input) INTEGER */
/*             The number of elementary reflectors whose product defines */
/*             the matrix Q. */
/*             N >= K >= 0,  if  SIDE = 'L'; */
/*             M >= K >= 0,  if  SIDE = 'R'. */

/*     P       (input) INTEGER */
/*             The order of the zero triagle (or the number of rows of */
/*             the zero trapezoid) in the matrix triangularized by SLICOT */
/*             Library routine MB04ID.  P >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,K) */
/*             On input, the elements in the rows  i+1:min(n,n-p-1+i)  of */
/*             the  i-th  column, and  TAU(i),  represent the orthogonal */
/*             reflector  H(i),  so that matrix  Q  is the product of */
/*             elementary reflectors:  Q = H(1) H(2) . . . H(k). */
/*             A is modified by the routine but restored on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array  A. */
/*             LDA >= max(1,N),  if  SIDE = 'L'; */
/*             LDA >= max(1,M),  if  SIDE = 'R'. */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             The scalar factors of the elementary reflectors. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix  C. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array  C.  LDC >= max(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,M),  if  SIDE = 'L'; */
/*             LDWORK >= MAX(1,N),  if  SIDE = 'R'. */
/*             For optimum performance LDWORK >= M*NB if SIDE = 'L', */
/*             or LDWORK >= N*NB if SIDE = 'R', where NB is the optimal */
/*             block size. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If  SIDE = 'L',  each elementary reflector  H(i)  modifies */
/*     n-p  elements of each column of  C,  for  i = 1:p+1,  and */
/*     n-i+1  elements, for  i = p+2:k. */
/*     If  SIDE = 'R',  each elementary reflector  H(i)  modifies */
/*     m-p  elements of each row of  C,  for  i = 1:p+1,  and */
/*     m-i+1  elements, for  i = p+2:k. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix operations, QR decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the scalar input arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! tran && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*k < 0 || left && *k > *n || ! left && *k > *m) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (left && *lda < max(1,*n) || ! left && *lda < max(1,*m)) {
	*info = -8;
    } else if (*ldc < max(1,*n)) {
	*info = -11;
    } else if (left && *ldwork < max(1,*m) || ! left && *ldwork < max(1,*n)) {
	*info = -13;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04IY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *k == 0 || left && *n < *p || ! left && *m < *p)
	     {
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    if (left) {
	wrkopt = (doublereal) (*m);
	if (tran) {

	    i__1 = min(*k,*p);
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Apply H(i) to C(i:i+n-p-1,1:m), from the left. */
/*              Workspace: need M. */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__2 = *n - *p;
		dlarf_(side, &i__2, m, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)1);
		a[i__ + i__ * a_dim1] = aii;
/* L10: */
	    }

	    if (*p <= min(*n,*k)) {

/*              Apply H(i) to C, i = p+1:k, from the left. */
/*              Workspace: need M;  prefer M*NB. */

		i__1 = *n - *p;
		i__2 = *k - *p;
		dormqr_(side, trans, &i__1, m, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[*p + 1 + c_dim1], 
			ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)1);
		wrkopt = max(wrkopt,dwork[1]);
	    }

	} else {

	    if (*p <= min(*n,*k)) {

/*              Apply H(i) to C, i = k:p+1:-1, from the left. */
/*              Workspace: need M;  prefer M*NB. */

		i__1 = *n - *p;
		i__2 = *k - *p;
		dormqr_(side, trans, &i__1, m, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[*p + 1 + c_dim1], 
			ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)1);
		wrkopt = max(wrkopt,dwork[1]);
	    }

	    for (i__ = min(*k,*p); i__ >= 1; --i__) {

/*              Apply H(i) to C(i:i+n-p-1,1:m), from the left. */
/*              Workspace: need M. */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__1 = *n - *p;
		dlarf_(side, &i__1, m, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)1);
		a[i__ + i__ * a_dim1] = aii;
/* L20: */
	    }
	}

    } else {

	wrkopt = (doublereal) (*n);
	if (tran) {

	    if (*p <= min(*m,*k)) {

/*              Apply H(i) to C, i = k:p+1:-1, from the right. */
/*              Workspace: need N;  prefer N*NB. */

		i__1 = *m - *p;
		i__2 = *k - *p;
		dormqr_(side, trans, n, &i__1, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[(*p + 1) * c_dim1 + 
			1], ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)
			1);
		wrkopt = max(wrkopt,dwork[1]);
	    }

	    for (i__ = min(*k,*p); i__ >= 1; --i__) {

/*              Apply H(i) to C(1:n,i:i+m-p-1), from the right. */
/*              Workspace: need N. */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__1 = *m - *p;
		dlarf_(side, n, &i__1, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (ftnlen)
			1);
		a[i__ + i__ * a_dim1] = aii;
/* L30: */
	    }

	} else {

	    i__1 = min(*k,*p);
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Apply H(i) to C(1:n,i:i+m-p-1), from the right. */
/*              Workspace: need N. */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__2 = *m - *p;
		dlarf_(side, n, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (ftnlen)
			1);
		a[i__ + i__ * a_dim1] = aii;
/* L40: */
	    }

	    if (*p <= min(*m,*k)) {

/*              Apply H(i) to C, i = p+1:k, from the right. */
/*              Workspace: need N;  prefer N*NB. */

		i__1 = *m - *p;
		i__2 = *k - *p;
		dormqr_(side, trans, n, &i__1, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[(*p + 1) * c_dim1 + 
			1], ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)
			1);
		wrkopt = max(wrkopt,dwork[1]);
	    }

	}
    }

    dwork[1] = wrkopt;
    return 0;

/* *** Last line of MB04IY *** */
} /* mb04iy_ */

