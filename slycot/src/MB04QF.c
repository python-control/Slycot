/* MB04QF.f -- translated by f2c (version 20100827).
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
static doublereal c_b11 = 0.;
static doublereal c_b49 = 1.;

/* Subroutine */ int mb04qf_(char *direct, char *storev, char *storew, 
	integer *n, integer *k, doublereal *v, integer *ldv, doublereal *w, 
	integer *ldw, doublereal *cs, doublereal *tau, doublereal *rs, 
	integer *ldrs, doublereal *t, integer *ldt, doublereal *dwork, ftnlen 
	direct_len, ftnlen storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer rs_dim1, rs_offset, t_dim1, t_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k2;
    static doublereal cm1;
    static integer pr1, pr2, pr3, ps1, ps2, ps3, pt11, pt12, pt13, pt21, pt22,
	     pt23, pt31, pt32, pt33;
    static doublereal vii, wii, taui;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical lcolv, lcolw;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);


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

/*     To form the triangular block factors R, S and T of a symplectic */
/*     block reflector SH, which is defined as a product of 2k */
/*     concatenated Householder reflectors and k Givens rotators, */

/*         SH = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*              diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                                .... */
/*              diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     The upper triangular blocks of the matrices */

/*                                 [ S1 ]       [ T11 T12 T13 ] */
/*         R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ], */
/*                                 [ S3 ]       [ T31 T32 T33 ] */

/*     with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular, */
/*     are stored rowwise in the arrays RS and T, respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DIRECT  CHARACTER*1 */
/*             This is a dummy argument, which is reserved for future */
/*             extensions of this subroutine. Not referenced. */

/*     STOREV  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder F(i) reflectors are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     STOREW  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder H(i) reflectors are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the Householder reflectors F(i) and H(i). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The number of Givens rotators.  K >= 1. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,N) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading N-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector F(i). */
/*             On entry with STOREV = 'R', the leading K-by-N part of */
/*             this array must contain in its i-th row the vector */
/*             which defines the elementary reflector F(i). */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,N),  if STOREV = 'C'; */
/*             LDV >= K,         if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,N) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading N-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector H(i). */
/*             On entry with STOREV = 'R', the leading K-by-N part of */
/*             this array must contain in its i-th row the vector */
/*             which defines the elementary reflector H(i). */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,N),  if STOREW = 'C'; */
/*             LDW >= K,         if STOREW = 'R'. */

/*     CS      (input) DOUBLE PRECISION array, dimension (2*K) */
/*             On entry, the first 2*K elements of this array must */
/*             contain the cosines and sines of the symplectic Givens */
/*             rotators G(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             On entry, the first K elements of this array must */
/*             contain the scalar factors of the elementary reflectors */
/*             F(i). */

/*     RS      (output) DOUBLE PRECISION array, dimension (K,6*K) */
/*             On exit, the leading K-by-6*K part of this array contains */
/*             the upper triangular matrices defining the factors R and */
/*             S of the symplectic block reflector SH. The (strictly) */
/*             lower portions of this array are not used. */

/*     LDRS    INTEGER */
/*             The leading dimension of the array RS.  LDRS >= K. */

/*     T       (output) DOUBLE PRECISION array, dimension (K,9*K) */
/*             On exit, the leading K-by-9*K part of this array contains */
/*             the upper triangular matrices defining the factor T of the */
/*             symplectic block reflector SH. The (strictly) lower */
/*             portions of this array are not used. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= K. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*K) */

/*     REFERENCES */

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires ( 4*K - 2 )*K*N + 19/3*K*K*K + 1/2*K*K */
/*     + 43/6*K - 4 floating point operations. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAEST). */

/*     KEYWORDS */

/*     Elementary matrix operations, orthogonal symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --cs;
    --tau;
    rs_dim1 = *ldrs;
    rs_offset = 1 + rs_dim1;
    rs -= rs_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --dwork;

    /* Function Body */
    if (*n == 0) {
	return 0;
    }

    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);

    k2 = *k + *k;
    pr1 = 0;
    pr2 = pr1 + *k;
    pr3 = pr2 + *k;
    ps1 = pr3 + *k;
    ps2 = ps1 + *k;
    ps3 = ps2 + *k;

    pt11 = 0;
    pt12 = pt11 + *k;
    pt13 = pt12 + *k;
    pt21 = pt13 + *k;
    pt22 = pt21 + *k;
    pt23 = pt22 + *k;
    pt31 = pt23 + *k;
    pt32 = pt31 + *k;
    pt33 = pt32 + *k;

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	taui = tau[i__];
	vii = v[i__ + i__ * v_dim1];
	v[i__ + i__ * v_dim1] = 1.;
	wii = w[i__ + i__ * w_dim1];
	w[i__ + i__ * w_dim1] = 1.;
	if (wii == 0.) {
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt11 + i__) * t_dim1] = 0.;
/* L10: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt21 + i__) * t_dim1] = 0.;
/* L20: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt31 + i__) * t_dim1] = 0.;
/* L30: */
	    }
	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		rs[j + (ps1 + i__) * rs_dim1] = 0.;
/* L40: */
	    }
	} else {

/*           Treat first Householder reflection. */

	    if (lcolv && lcolw) {

/*              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i). */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -wii;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &w[i__ + w_dim1], 
			ldw, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[1],
			 &c__1, (ftnlen)9);

/*              Compute t2 = -wii * V(i:n,1:i-1)' * W(i:n,i). */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -wii;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[*k 
			+ 1], &c__1, (ftnlen)9);
	    } else if (lcolv) {

/*              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'. */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -wii;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -wii * V(i:n,1:i-1)' * W(i,i:n)'. */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -wii;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[*k + 
			1], &c__1, (ftnlen)9);
	    } else if (lcolw) {

/*              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i). */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -wii;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &w[i__ + w_dim1], 
			ldw, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[1],
			 &c__1, (ftnlen)9);

/*              Compute t2 = -wii * V(1:i-1,i:n) * W(i:n,i). */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -wii;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &
			dwork[*k + 1], &c__1, (ftnlen)12);
	    } else {

/*              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'. */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -wii;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -wii * V(1:i-1,i:n) * W(i,i:n)'. */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -wii;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[*
			k + 1], &c__1, (ftnlen)12);
	    }

/*           T11(1:i-1,i) := T11(1:i-1,1:i-1)*t1 + T13(1:i-1,1:i-1)*t2 */

	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[1], &c__1, &t[(pt11 + i__) * t_dim1 + 1], &
		    c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt11 + 1) *
		     t_dim1 + 1], ldt, &t[(pt11 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1]
		    , &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) *
		     t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    daxpy_(&i__2, &c_b49, &t[(pt13 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt11 + i__) * t_dim1 + 1], &c__1);
	    t[i__ + (pt11 + i__) * t_dim1] = -wii;

	    if (i__ > 1) {

/*              T21(1:i-1,i) := T21(1:i-1,1:i-1)*t1 + T23(1:i-1,1:i-1)*t2 */

		i__2 = i__ - 2;
		dcopy_(&i__2, &dwork[2], &c__1, &t[(pt21 + i__) * t_dim1 + 1],
			 &c__1);
		i__2 = i__ - 2;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 
			2) * t_dim1 + 1], ldt, &t[(pt21 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		t[i__ - 1 + (pt21 + i__) * t_dim1] = 0.;
		i__2 = i__ - 1;
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * t_dim1 
			+ 1], &c__1);
		i__2 = i__ - 1;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 
			1) * t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		i__2 = i__ - 1;
		daxpy_(&i__2, &c_b49, &t[(pt23 + i__) * t_dim1 + 1], &c__1, &
			t[(pt21 + i__) * t_dim1 + 1], &c__1);

/*              T31(1:i-1,i) := T31(1:i-1,1:i-1)*t1 + T33(1:i-1,1:i-1)*t2 */

		i__2 = i__ - 2;
		dcopy_(&i__2, &dwork[2], &c__1, &t[(pt31 + i__) * t_dim1 + 1],
			 &c__1);
		i__2 = i__ - 2;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 
			2) * t_dim1 + 1], ldt, &t[(pt31 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		t[i__ - 1 + (pt31 + i__) * t_dim1] = 0.;
		i__2 = i__ - 1;
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * t_dim1 
			+ 1], &c__1);
		i__2 = i__ - 1;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 
			1) * t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		i__2 = i__ - 1;
		daxpy_(&i__2, &c_b49, &t[(pt33 + i__) * t_dim1 + 1], &c__1, &
			t[(pt31 + i__) * t_dim1 + 1], &c__1);

/*              S1(1:i-1,i) := S1(1:i-1,1:i-1)*t1 + S3(1:i-1,1:i-1)*t2 */

		i__2 = i__ - 2;
		dcopy_(&i__2, &dwork[2], &c__1, &rs[(ps1 + i__) * rs_dim1 + 1]
			, &c__1);
		i__2 = i__ - 2;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 
			2) * rs_dim1 + 1], ldrs, &rs[(ps1 + i__) * rs_dim1 + 
			1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		rs[i__ - 1 + (ps1 + i__) * rs_dim1] = 0.;
		i__2 = i__ - 1;
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &rs[(ps3 + i__) * 
			rs_dim1 + 1], &c__1);
		i__2 = i__ - 1;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 
			1) * rs_dim1 + 1], ldrs, &rs[(ps3 + i__) * rs_dim1 + 
			1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
		i__2 = i__ - 1;
		daxpy_(&i__2, &c_b49, &rs[(ps3 + i__) * rs_dim1 + 1], &c__1, &
			rs[(ps1 + i__) * rs_dim1 + 1], &c__1);
	    }
	}

/*        Treat Givens rotation. */

	cm1 = cs[(i__ << 1) - 1] - 1.;
	if (lcolw) {
	    dcopy_(&i__, &w[i__ + w_dim1], ldw, &dwork[1], &c__1);
	} else {
	    dcopy_(&i__, &w[i__ * w_dim1 + 1], &c__1, &dwork[1], &c__1);
	}
	if (lcolv) {
	    i__2 = i__ - 1;
	    dcopy_(&i__2, &v[i__ + v_dim1], ldv, &dwork[*k + 1], &c__1);
	} else {
	    i__2 = i__ - 1;
	    dcopy_(&i__2, &v[i__ * v_dim1 + 1], &c__1, &dwork[*k + 1], &c__1);
	}

/*        R1(1:i,i) = T11(1:i,1:i) * dwork(1:i) */
/*                    + [ T13(1:i-1,1:i-1) * dwork(k+1:k+i-1); 0 ] */

	dcopy_(&i__, &dwork[1], &c__1, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1);
	dtrmv_("Upper", "No transpose", "Non-unit", &i__, &t[(pt11 + 1) * 
		t_dim1 + 1], ldt, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1], &
		c__1);
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) * 
		t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	daxpy_(&i__2, &c_b49, &t[(pt13 + i__) * t_dim1 + 1], &c__1, &rs[(pr1 
		+ i__) * rs_dim1 + 1], &c__1);

/*        R2(1:i-1,i) = T21(1:i-1,2:i) * W(i,2:i) */
/*                      + T23(1:i-1,1:i-1) * V(i,1:i-1) */

	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1)
		;
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 2) * 
		t_dim1 + 1], ldt, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * t_dim1 + 1], &
		c__1);
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 1) * 
		t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	daxpy_(&i__2, &c_b49, &t[(pt23 + i__) * t_dim1 + 1], &c__1, &rs[(pr2 
		+ i__) * rs_dim1 + 1], &c__1);

/*        R3(1:i-1,i) = T31(1:i-1,2:i) * dwork(2:i) */
/*                      + T33(1:i-1,1:i-1) * dwork(k+1:k+i-1) */

	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1)
		;
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 2) * 
		t_dim1 + 1], ldt, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * t_dim1 + 1], &
		c__1);
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 1) * 
		t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	daxpy_(&i__2, &c_b49, &t[(pt33 + i__) * t_dim1 + 1], &c__1, &rs[(pr3 
		+ i__) * rs_dim1 + 1], &c__1);

/*        S2(1:i-1,i) = S1(1:i-1,2:i) * dwork(2:i) */
/*                      + S3(1:i-1,1:i-1) * dwork(k+1:k+i-1) */

	i__2 = i__ - 1;
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1)
		;
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 2) * 
		rs_dim1 + 1], ldrs, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 1) * 
		rs_dim1 + 1], ldrs, &dwork[*k + 1], &c__1, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);
	i__2 = i__ - 1;
	daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &rs[(ps2 + i__) * 
		rs_dim1 + 1], &c__1);
	rs[i__ + (ps2 + i__) * rs_dim1] = -cs[i__ * 2];

/*        T12(1:i,i) = [ R1(1:i-1,1:i-1)*S2(1:i-1,i); 0 ] */
/*                     + (c-1) * R1(1:i,i) */

	i__2 = i__ - 1;
	dcopy_(&i__2, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, &t[(pt12 + i__) *
		 t_dim1 + 1], &c__1);
	i__2 = i__ - 1;
	dscal_(&i__2, &cm1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1);
	i__2 = i__ - 1;
	dscal_(&i__2, &cs[i__ * 2], &t[(pt12 + i__) * t_dim1 + 1], &c__1);
	i__2 = i__ - 1;
	dcopy_(&i__2, &t[(pt12 + i__) * t_dim1 + 1], &c__1, &t[(pt22 + i__) * 
		t_dim1 + 1], &c__1);
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(pr1 + 1) * 
		rs_dim1 + 1], ldrs, &t[(pt12 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	t[i__ + (pt12 + i__) * t_dim1] = 0.;
	daxpy_(&i__, &cm1, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1, &t[(pt12 + 
		i__) * t_dim1 + 1], &c__1);

/*        T22(1:i-1,i) = R2(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R2(1:i-1,i) */

	if (i__ > 1) {
	    i__2 = i__ - 2;
	    dcopy_(&i__2, &t[(pt22 + i__) * t_dim1 + 2], &c__1, &t[(pt32 + 
		    i__) * t_dim1 + 1], &c__1);
	}
	i__2 = i__ - 1;
	dtrmv_("Upper", "No transpose", "Unit diagonal", &i__2, &rs[(pr2 + 1) 
		* rs_dim1 + 1], ldrs, &t[(pt22 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)13);
	i__2 = i__ - 1;
	daxpy_(&i__2, &cm1, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1, &t[(pt22 + 
		i__) * t_dim1 + 1], &c__1);
	t[i__ + (pt22 + i__) * t_dim1] = cm1;

/*        T32(1:i-1,i) = R3(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R3(1:i-1,i) */

	if (i__ > 1) {
	    i__2 = i__ - 2;
	    dtrmv_("Upper", "No transpose", "Non-Unit", &i__2, &rs[(pr3 + 2) *
		     rs_dim1 + 1], ldrs, &t[(pt32 + i__) * t_dim1 + 1], &c__1,
		     (ftnlen)5, (ftnlen)12, (ftnlen)8);
	    t[i__ - 1 + (pt32 + i__) * t_dim1] = 0.;
	    i__2 = i__ - 1;
	    daxpy_(&i__2, &cm1, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1, &t[(
		    pt32 + i__) * t_dim1 + 1], &c__1);
	}

	if (taui == 0.) {
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt13 + i__) * t_dim1] = 0.;
/* L50: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt23 + i__) * t_dim1] = 0.;
/* L60: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		t[j + (pt33 + i__) * t_dim1] = 0.;
/* L70: */
	    }
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		rs[j + (ps3 + i__) * rs_dim1] = 0.;
/* L80: */
	    }
	} else {

/*           Treat second Householder reflection. */

	    if (lcolv && lcolw) {

/*              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i:n,i). */

		i__2 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("Transpose", &i__2, &i__, &d__1, &w[i__ + w_dim1], ldw,
			 &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[1], &
			c__1, (ftnlen)9);

/*              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i). */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -taui;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[k2 
			+ 1], &c__1, (ftnlen)9);
	    } else if (lcolv) {

/*              Compute t1 = -tau(i) * W(1:i,i:n) * V(i:n,i). */

		i__2 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("No Transpose", &i__, &i__2, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &
			dwork[1], &c__1, (ftnlen)12);

/*              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i). */

		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		d__1 = -taui;
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[k2 
			+ 1], &c__1, (ftnlen)9);
	    } else if (lcolw) {

/*              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i,i:n)'. */

		i__2 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("Transpose", &i__2, &i__, &d__1, &w[i__ + w_dim1], ldw,
			 &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[1], &
			c__1, (ftnlen)9);

/*              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'. */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			k2 + 1], &c__1, (ftnlen)12);
	    } else {

/*              Compute t1 = -tau(i) * W(1:i,i:n) * V(i,i:n)'. */

		i__2 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("No Transpose", &i__, &i__2, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'. */

		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		d__1 = -taui;
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			k2 + 1], &c__1, (ftnlen)12);
	    }

/*           T13(1:i,i) := T11(1:i,1:i)*t1 - tau(i)*T12(1:i,i) */
/*                                         + [T13(1:i-1,1:i-1)*t2;0] */

	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1]
		    , &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) *
		     t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
	    t[i__ + (pt13 + i__) * t_dim1] = 0.;
	    dcopy_(&i__, &dwork[1], &c__1, &dwork[*k + 1], &c__1);
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__, &t[(pt11 + 1) * 
		    t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
	    daxpy_(&i__, &c_b49, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * 
		    t_dim1 + 1], &c__1);
	    d__1 = -taui;
	    daxpy_(&i__, &d__1, &t[(pt12 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt13 + i__) * t_dim1 + 1], &c__1);

/*           T23(1:i,i) := T21(1:i,1:i)*t1 - tau(i)*T22(1:i,i) */
/*                                         + [T23(1:i-1,1:i-1)*t2;0] */

	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt23 + i__) * t_dim1 + 1]
		    , &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 1) *
		     t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
	    t[i__ + (pt23 + i__) * t_dim1] = 0.;
	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[2], &c__1, &dwork[*k + 1], &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 2) *
		     t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * 
		    t_dim1 + 1], &c__1);
	    d__1 = -taui;
	    daxpy_(&i__, &d__1, &t[(pt22 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt23 + i__) * t_dim1 + 1], &c__1);

/*           T33(1:i,i) := T31(1:i,1:i)*t1 - tau(i)*T32(1:i,i) */
/*                                         + [T33(1:i-1,1:i-1)*t2;0] */

	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt33 + i__) * t_dim1 + 1]
		    , &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 1) *
		     t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[2], &c__1, &dwork[*k + 1], &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 2) *
		     t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * 
		    t_dim1 + 1], &c__1);
	    i__2 = i__ - 1;
	    d__1 = -taui;
	    daxpy_(&i__2, &d__1, &t[(pt32 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt33 + i__) * t_dim1 + 1], &c__1);
	    t[i__ + (pt33 + i__) * t_dim1] = -taui;

/*           S3(1:i,i) := S1(1:i,1:i)*t1 - tau(i)*S2(1:i,i) */
/*                                       + [S3(1:i-1,1:i-1)*t2;0] */

	    i__2 = i__ - 1;
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &rs[(ps3 + i__) * rs_dim1 + 
		    1], &c__1);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 1) *
		     rs_dim1 + 1], ldrs, &rs[(ps3 + i__) * rs_dim1 + 1], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	    i__2 = i__ - 1;
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 2) *
		     rs_dim1 + 1], ldrs, &dwork[2], &c__1, (ftnlen)5, (ftnlen)
		    12, (ftnlen)8);
	    i__2 = i__ - 1;
	    daxpy_(&i__2, &c_b49, &dwork[2], &c__1, &rs[(ps3 + i__) * rs_dim1 
		    + 1], &c__1);
	    rs[i__ + (ps3 + i__) * rs_dim1] = 0.;
	    d__1 = -taui;
	    daxpy_(&i__, &d__1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, &rs[(
		    ps3 + i__) * rs_dim1 + 1], &c__1);
	}
	w[i__ + i__ * w_dim1] = wii;
	v[i__ + i__ * v_dim1] = vii;
/* L90: */
    }

    return 0;
/* *** Last line of MB04QF *** */
} /* mb04qf_ */

