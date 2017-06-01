/* MB04QU.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04qu_(char *tranc, char *trand, char *tranq, char *
	storev, char *storew, integer *m, integer *n, integer *k, doublereal *
	v, integer *ldv, doublereal *w, integer *ldw, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *cs, 
	doublereal *tau, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen tranc_len, ftnlen trand_len, ftnlen tranq_len, ftnlen 
	storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, d_dim1, d_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal nu;
    static logical ltrc, ltrd;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltrq;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lcolv, lcolw;
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

/*     To overwrite general real m-by-n matrices C and D, or their */
/*     transposes, with */

/*               [ op(C) ] */
/*         Q  *  [       ]   if TRANQ = 'N', or */
/*               [ op(D) ] */

/*          T    [ op(C) ] */
/*         Q  *  [       ]   if TRANQ = 'T', */
/*               [ op(D) ] */

/*     where Q is defined as the product of symplectic reflectors and */
/*     Givens rotators, */

/*         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*             diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                               .... */
/*             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     Unblocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANC   CHARACTER*1 */
/*             Specifies the form of op( C ) as follows: */
/*             = 'N':  op( C ) = C; */
/*             = 'T':  op( C ) = C'; */
/*             = 'C':  op( C ) = C'. */

/*     STOREV  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder reflectors contained in V are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     STOREW  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder reflectors contained in W are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     TRAND   CHARACTER*1 */
/*             Specifies the form of op( D ) as follows: */
/*             = 'N':  op( D ) = D; */
/*             = 'T':  op( D ) = D'; */
/*             = 'C':  op( D ) = D'. */

/*     TRANQ   CHARACTER*1 */
/*             = 'N':  apply Q; */
/*             = 'T':  apply Q'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices op(C) and op(D). */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices op(C) and op(D). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The number of elementary reflectors whose product defines */
/*             the matrix Q.  M >= K >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,M) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflectors F(i). */
/*             On entry with STOREV = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflectors F(i). */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M),  if STOREV = 'C'; */
/*             LDV >= MAX(1,K),  if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,M) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflectors H(i). */
/*             On entry with STOREW = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflectors H(i). */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,M),  if STOREW = 'C'; */
/*             LDW >= MAX(1,K),  if STOREW = 'R'. */

/*     C       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDC,N) if TRANC = 'N', */
/*                     (LDC,M) if TRANC = 'T' or TRANC = 'C' */
/*             On entry with TRANC = 'N', the leading M-by-N part of */
/*             this array must contain the matrix C. */
/*             On entry with TRANC = 'C' or TRANC = 'T', the leading */
/*             N-by-M part of this array must contain the transpose of */
/*             the matrix C. */
/*             On exit with TRANC = 'N', the leading M-by-N part of */
/*             this array contains the updated matrix C. */
/*             On exit with TRANC = 'C' or TRANC = 'T', the leading */
/*             N-by-M part of this array contains the transpose of the */
/*             updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= MAX(1,M),  if TRANC = 'N'; */
/*             LDC >= MAX(1,N),  if TRANC = 'T' or TRANC = 'C'. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDD,N) if TRAND = 'N', */
/*                     (LDD,M) if TRAND = 'T' or TRAND = 'C' */
/*             On entry with TRAND = 'N', the leading M-by-N part of */
/*             this array must contain the matrix D. */
/*             On entry with TRAND = 'C' or TRAND = 'T', the leading */
/*             N-by-M part of this array must contain the transpose of */
/*             the matrix D. */
/*             On exit with TRAND = 'N', the leading M-by-N part of */
/*             this array contains the updated matrix D. */
/*             On exit with TRAND = 'C' or TRAND = 'T', the leading */
/*             N-by-M part of this array contains the transpose of the */
/*             updated matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= MAX(1,M),  if TRAND = 'N'; */
/*             LDD >= MAX(1,N),  if TRAND = 'T' or TRAND = 'C'. */

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
/*             On exit, if  INFO = -20,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSMSQ). */

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
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --cs;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);
    ltrc = lsame_(tranc, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranc, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrd = lsame_(trand, "T", (ftnlen)1, (ftnlen)1) || lsame_(trand, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrq = lsame_(tranq, "T", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (ltrc || lsame_(tranc, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (ltrd || lsame_(trand, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (ltrq || lsame_(tranq, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (lcolv || lsame_(storev, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (! (lcolw || lsame_(storew, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0) {
	*info = -7;
    } else if (*k < 0 || *k > *m) {
	*info = -8;
    } else if (lcolv && *ldv < max(1,*m) || ! lcolv && *ldv < max(1,*k)) {
	*info = -10;
    } else if (lcolw && *ldw < max(1,*m) || ! lcolw && *ldw < max(1,*k)) {
	*info = -12;
    } else if (ltrc && *ldc < max(1,*n) || ! ltrc && *ldc < max(1,*m)) {
	*info = -14;
    } else if (ltrd && *ldd < max(1,*n) || ! ltrd && *ldd < max(1,*m)) {
	*info = -16;
    } else if (*ldwork < max(1,*n)) {
	dwork[1] = (doublereal) max(1,*n);
	*info = -20;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04QU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*k,*m);
    if (min(i__1,*n) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    if (ltrq) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Apply H(I) to C(I:M,:) and D(I:M,:) from the left. */

	    nu = w[i__ + i__ * w_dim1];
	    w[i__ + i__ * w_dim1] = 1.;
	    if (lcolw) {
		if (ltrc) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4)
			    ;
		}
		if (ltrd) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4)
			    ;
		}
	    } else {
		if (ltrc) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4);
		}
		if (ltrd) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4);
		}
	    }
	    w[i__ + i__ * w_dim1] = nu;

/*           Apply G(i) to C(I,:) and D(I,:) from the left. */

	    if (ltrc && ltrd) {
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ * d_dim1 + 1]
			, &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
	    } else if (ltrc) {
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ + d_dim1], 
			ldd, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
	    } else if (ltrd) {
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ * d_dim1 + 1], &
			c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
	    } else {
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ + d_dim1], ldd, &
			cs[(i__ << 1) - 1], &cs[i__ * 2]);
	    }

/*           Apply F(I) to C(I:M,:) and D(I:M,:) from the left. */

	    nu = v[i__ + i__ * v_dim1];
	    v[i__ + i__ * v_dim1] = 1.;
	    if (lcolv) {
		if (ltrc) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
		}
		if (ltrd) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
		}
	    } else {
		if (ltrc) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
		}
		if (ltrd) {
		    i__2 = *m - i__ + 1;
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__2 = *m - i__ + 1;
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
		}
	    }
	    v[i__ + i__ * v_dim1] = nu;
/* L10: */
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to C(I:M,:) and D(I:M,:) from the left. */

	    nu = v[i__ + i__ * v_dim1];
	    v[i__ + i__ * v_dim1] = 1.;
	    if (lcolv) {
		if (ltrc) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
		}
		if (ltrd) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
		}
	    } else {
		if (ltrc) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
		}
		if (ltrd) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
		}
	    }
	    v[i__ + i__ * v_dim1] = nu;

/*           Apply G(i) to C(I,:) and D(I,:) from the left. */

	    if (ltrc && ltrd) {
		d__1 = -cs[i__ * 2];
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ * d_dim1 + 1]
			, &c__1, &cs[(i__ << 1) - 1], &d__1);
	    } else if (ltrc) {
		d__1 = -cs[i__ * 2];
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ + d_dim1], 
			ldd, &cs[(i__ << 1) - 1], &d__1);
	    } else if (ltrd) {
		d__1 = -cs[i__ * 2];
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ * d_dim1 + 1], &
			c__1, &cs[(i__ << 1) - 1], &d__1);
	    } else {
		d__1 = -cs[i__ * 2];
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ + d_dim1], ldd, &
			cs[(i__ << 1) - 1], &d__1);
	    }

/*           Apply H(I) to C(I:M,:) and D(I:M,:) from the left. */

	    nu = w[i__ + i__ * w_dim1];
	    w[i__ + i__ * w_dim1] = 1.;
	    if (lcolw) {
		if (ltrc) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4)
			    ;
		}
		if (ltrd) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4)
			    ;
		}
	    } else {
		if (ltrc) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4);
		}
		if (ltrd) {
		    i__1 = *m - i__ + 1;
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__1 = *m - i__ + 1;
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4);
		}
	    }
	    w[i__ + i__ * w_dim1] = nu;
/* L20: */
	}
    }

    dwork[1] = (doublereal) max(1,*n);
/* *** Last line of MB04QU *** */
    return 0;
} /* mb04qu_ */

