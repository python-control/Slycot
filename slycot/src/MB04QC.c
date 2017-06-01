/* MB04QC.f -- translated by f2c (version 20100827).
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
static doublereal c_b22 = 1.;
static doublereal c_b53 = 0.;
static doublereal c_b367 = -1.;

/* Subroutine */ int mb04qc_(char *struct__, char *trana, char *tranb, char *
	tranq, char *direct, char *storev, char *storew, integer *m, integer *
	n, integer *k, doublereal *v, integer *ldv, doublereal *w, integer *
	ldw, doublereal *rs, integer *ldrs, doublereal *t, integer *ldt, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	dwork, ftnlen struct_len, ftnlen trana_len, ftnlen tranb_len, ftnlen 
	tranq_len, ftnlen direct_len, ftnlen storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, rs_dim1, rs_offset, t_dim1, 
	    t_offset, v_dim1, v_offset, w_dim1, w_offset, i__1;

    /* Local variables */
    static integer i__, pr1, pr2, pr3, ps1, ps2, ps3, pt11, pt12, pt13, pt21, 
	    pt22, pt23, pt31, pt32, pt33, pdw1, pdw2, pdw3, pdw4, pdw5, pdw6, 
	    pdw7, pdw8, pdw9;
    static logical la1b1;
    static doublereal fact;
    static logical ltra, ltrb, ltrq;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itemp;
    static logical lcolv, lcolw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);


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

/*     To apply the orthogonal symplectic block reflector */

/*              [  I+V*T*V'  V*R*S*V'  ] */
/*         Q =  [                      ] */
/*              [ -V*R*S*V'  I+V*T*V'  ] */

/*     or its transpose to a real 2m-by-n matrix [ op(A); op(B) ] from */
/*     the left. */
/*     The k-by-k upper triangular blocks of the matrices */

/*                                 [ S1 ]       [ T11 T12 T13 ] */
/*         R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ], */
/*                                 [ S3 ]       [ T31 T32 T33 ] */

/*     with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular, */
/*     are stored rowwise in the arrays RS and T, respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STRUCT  CHARACTER*1 */
/*             Specifies the structure of the first blocks of A and B: */
/*             = 'Z':  the leading K-by-N submatrices of op(A) and op(B) */
/*                     are (implicitly) assumed to be zero; */
/*             = 'N';  no structure to mention. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op( A ) as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op( B ) as follows: */
/*             = 'N':  op( B ) = B; */
/*             = 'T':  op( B ) = B'; */
/*             = 'C':  op( B ) = B'. */

/*     DIRECT  CHARACTER*1 */
/*             This is a dummy argument, which is reserved for future */
/*             extensions of this subroutine. Not referenced. */

/*     TRANQ   CHARACTER*1 */
/*             = 'N':  apply Q; */
/*             = 'T':  apply Q'. */

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

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices op(A) and op(B). */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices op(A) and op(B). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The order of the triangular matrices defining R, S and T. */
/*             M >= K >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,M) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflector used to form parts of Q. */
/*             On entry with STOREV = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflector used to form parts of Q. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M),  if STOREV = 'C'; */
/*             LDV >= MAX(1,K),  if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,M) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflector used to form parts of Q. */
/*             On entry with STOREW = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflector used to form parts of Q. */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,M),  if STOREW = 'C'; */
/*             LDW >= MAX(1,K),  if STOREW = 'R'. */

/*     RS      (input) DOUBLE PRECISION array, dimension (K,6*K) */
/*             On entry, the leading K-by-6*K part of this array must */
/*             contain the upper triangular matrices defining the factors */
/*             R and S of the symplectic block reflector Q. The */
/*             (strictly) lower portions of this array are not */
/*             referenced. */

/*     LDRS    INTEGER */
/*             The leading dimension of the array RS.  LDRS >= MAX(1,K). */

/*     T       (input) DOUBLE PRECISION array, dimension (K,9*K) */
/*             On entry, the leading K-by-9*K part of this array must */
/*             contain the upper triangular matrices defining the factor */
/*             T of the symplectic block reflector Q. The (strictly) */
/*             lower portions of this array are not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,K). */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDA,N) if TRANA = 'N', */
/*                     (LDA,M) if TRANA = 'C' or TRANA = 'T' */
/*             On entry with TRANA = 'N', the leading M-by-N part of this */
/*             array must contain the matrix A. */
/*             On entry with TRANA = 'T' or TRANA = 'C', the leading */
/*             N-by-M part of this array must contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,M),  if TRANA = 'N'; */
/*             LDA >= MAX(1,N),  if TRANA = 'C' or TRANA = 'T'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDB,N) if TRANB = 'N', */
/*                     (LDB,M) if TRANB = 'C' or TRANB = 'T' */
/*             On entry with TRANB = 'N', the leading M-by-N part of this */
/*             array must contain the matrix B. */
/*             On entry with TRANB = 'T' or TRANB = 'C', the leading */
/*             N-by-M part of this array must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,M),  if TRANB = 'N'; */
/*             LDB >= MAX(1,N),  if TRANB = 'C' or TRANB = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK), where */
/*             LDWORK >= 8*N*K,   if STRUCT = 'Z', */
/*             LDWORK >= 9*N*K,   if STRUCT = 'N'. */

/*     REFERENCES */

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 16*( M - K )*N + ( 26*K - 4 )*K*N floating */
/*     point operations if STRUCT = 'Z' and additional ( 12*K + 2 )*K*N */
/*     floating point operations if STRUCT = 'N'. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAESB). */

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
    rs_dim1 = *ldrs;
    rs_offset = 1 + rs_dim1;
    rs -= rs_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dwork;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }
    la1b1 = lsame_(struct__, "N", (ftnlen)1, (ftnlen)1);
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);
    ltra = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrq = lsame_(tranq, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq, "C", (
	    ftnlen)1, (ftnlen)1);

    pr1 = 1;
    pr2 = pr1 + *k;
    pr3 = pr2 + *k;
    ps1 = pr3 + *k;
    ps2 = ps1 + *k;
    ps3 = ps2 + *k;
    pt11 = 1;
    pt12 = pt11 + *k;
    pt13 = pt12 + *k;
    pt21 = pt13 + *k;
    pt22 = pt21 + *k;
    pt23 = pt22 + *k;
    pt31 = pt23 + *k;
    pt32 = pt31 + *k;
    pt33 = pt32 + *k;
    pdw1 = 1;
    pdw2 = pdw1 + *n * *k;
    pdw3 = pdw2 + *n * *k;
    pdw4 = pdw3 + *n * *k;
    pdw5 = pdw4 + *n * *k;
    pdw6 = pdw5 + *n * *k;
    pdw7 = pdw6 + *n * *k;
    pdw8 = pdw7 + *n * *k;
    pdw9 = pdw8 + *n * *k;

/*     Update the matrix A. */

    if (la1b1) {

/*        NZ1) DW7 := A1' */

	if (ltra) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &a[i__ * a_dim1 + 1], &c__1, &dwork[pdw7 + (i__ - 1)
			 * *n], &c__1);
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(k, &a[i__ * a_dim1 + 1], &c__1, &dwork[pdw7 + i__ - 1],
			 n);
/* L20: */
	    }
	}

/*        NZ2) DW1 := DW7*W1 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw1], &c__1);
	if (lcolw) {
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
	} else {
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
	}

/*        NZ3) DW2 := DW7*V1 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw2], &c__1);
	if (lcolv) {
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw2], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
	} else {
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw2], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
	}
	fact = 1.;
    } else {
	fact = 0.;
    }

/*     1) DW1 := A2'*W2 */

    if (*m > *k) {
	if (ltra && lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[(*
		    k + 1) * a_dim1 + 1], lda, &w[*k + 1 + w_dim1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)12);
	} else if (ltra) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &a[(*k + 
		    1) * a_dim1 + 1], lda, &w[(*k + 1) * w_dim1 + 1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)9);
	} else if (lcolw) {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[*k + 
		    1 + a_dim1], lda, &w[*k + 1 + w_dim1], ldw, &fact, &dwork[
		    pdw1], n, (ftnlen)9, (ftnlen)12);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, &w[(*k + 1) * w_dim1 + 1], ldw, &fact, &
		    dwork[pdw1], n, (ftnlen)9, (ftnlen)9);
	}
    } else if (! la1b1) {
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw1], n, (ftnlen)3);
    }

/*     2) DW2 := A2'*V2 */

    if (*m > *k) {
	if (ltra && lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[(*
		    k + 1) * a_dim1 + 1], lda, &v[*k + 1 + v_dim1], ldv, &
		    fact, &dwork[pdw2], n, (ftnlen)12, (ftnlen)12);
	} else if (ltra) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &a[(*k + 
		    1) * a_dim1 + 1], lda, &v[(*k + 1) * v_dim1 + 1], ldv, &
		    fact, &dwork[pdw2], n, (ftnlen)12, (ftnlen)9);
	} else if (lcolv) {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[*k + 
		    1 + a_dim1], lda, &v[*k + 1 + v_dim1], ldv, &fact, &dwork[
		    pdw2], n, (ftnlen)9, (ftnlen)12);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, &v[(*k + 1) * v_dim1 + 1], ldv, &fact, &
		    dwork[pdw2], n, (ftnlen)9, (ftnlen)9);
	}
    } else if (! la1b1) {
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw2], n, (ftnlen)3);
    }

    if (ltrq) {

/*        3) DW3 := DW1*T11 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        4) DW4 := DW2*T31 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        5) DW3 := DW3 + DW4 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

	if (la1b1) {

/*           NZ4) DW8 := DW7*T21 */

	    i__1 = *n * (*k - 1);
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	    i__1 = *k - 1;
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ5) DW3 := DW3 + DW8 */

	    i__1 = *n * (*k - 1);
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &
		    c__1);
	}

/*        6) DW4 := DW1*T12 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt12 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        7) DW5 := DW2*T32 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        8) DW4 := DW4 + DW5 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4 + *n], &c__1);

	if (la1b1) {

/*           NZ6) DW8 := DW7*T22 */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ7) DW4 := DW4 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
	}

/*        9) DW5 := DW2*T33 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        10) DW6 := DW1*T13 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw6], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        11) DW5 := DW5 + DW6 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);

	if (la1b1) {

/*           NZ8) DW8 := DW7*T23 */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt23 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ9) DW5 := DW5 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);
	}

/*        12) DW1 := DW1*R1 */

	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[pr1 * rs_dim1 + 1], ldrs, &dwork[pdw1], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        13) DW2 := DW2*R3 */

	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        14) DW1 := DW1 + DW2 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw2], &c__1, &dwork[pdw1 + *n], &c__1);

	if (la1b1) {

/*           NZ10) DW7 := DW7*R2 */

	    dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &
		    rs[pr2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);

/*           NZ11) DW1 := DW1 + DW7 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw1], &c__1);
	}

/*        Swap Pointers PDW1 <-> PDW2 */

	itemp = pdw2;
	pdw2 = pdw1;
	pdw1 = itemp;
    } else {

/*        3) DW3 := DW1*T11' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        4) DW4 := DW2*T13' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        5) DW3 := DW3 + DW4 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

	if (la1b1) {

/*           NZ4) DW8 := DW7*T12' */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt12 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ5) DW3 := DW3 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);
	}

/*        6) DW4 := DW2*T23' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt23 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        7) DW5 := DW1*T21' */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw5], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        8) DW4 := DW4 + DW5 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

	if (la1b1) {

/*           NZ6) DW8 := DW7*T22' */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ7) DW4 := DW4 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
	}

/*        9) DW5 := DW2*T33' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        10) DW6 := DW1*T31' */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw6], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        11) DW5 := DW5 + DW6 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);

	if (la1b1) {

/*           NZ8) DW8 := DW7*T32' */

	    i__1 = *n * (*k - 1);
	    dcopy_(&i__1, &dwork[pdw7 + *n], &c__1, &dwork[pdw8], &c__1);
	    i__1 = *k - 1;
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ9) DW5 := DW5 + DW8 */

	    i__1 = *n * (*k - 1);
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);
	}

/*        12) DW1 := DW1*S1' */

	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw1 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        13) DW2 := DW2*S3' */

	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		ps3 * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        14) DW2 := DW1 + DW2 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw1 + *n], &c__1, &dwork[pdw2], &c__1);

	if (la1b1) {

/*           NZ10) DW7 := DW7*S2' */

	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ11) DW2 := DW2 + DW7 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw2], &c__1);
	}
    }

    if (la1b1) {

/*        NZ12) DW9 := B1' */

	if (ltrb) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[pdw9 + (i__ - 1)
			 * *n], &c__1);
/* L30: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(k, &b[i__ * b_dim1 + 1], &c__1, &dwork[pdw9 + i__ - 1],
			 n);
/* L40: */
	    }
	}

/*        NZ13) DW1 := DW9*W1 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw1], &c__1);
	if (lcolw) {
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
	} else {
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
	}

/*        NZ14) DW6 := DW9*V1 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw6], &c__1);
	if (lcolv) {
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
	} else {
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
	}
    }

/*     15) DW1 := B2'*W2 */

    if (*m > *k) {
	if (ltrb && lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[(*
		    k + 1) * b_dim1 + 1], ldb, &w[*k + 1 + w_dim1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)12);
	} else if (ltrb) {

/*           Critical Position */

	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &b[(*k + 
		    1) * b_dim1 + 1], ldb, &w[(*k + 1) * w_dim1 + 1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)9);
	} else if (lcolw) {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[*k + 
		    1 + b_dim1], ldb, &w[*k + 1 + w_dim1], ldw, &fact, &dwork[
		    pdw1], n, (ftnlen)9, (ftnlen)12);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, &w[(*k + 1) * w_dim1 + 1], ldw, &fact, &
		    dwork[pdw1], n, (ftnlen)9, (ftnlen)9);
	}
    } else if (! la1b1) {
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw1], n, (ftnlen)3);
    }

/*     16) DW6 := B2'*V2 */

    if (*m > *k) {
	if (ltrb && lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[(*
		    k + 1) * b_dim1 + 1], ldb, &v[*k + 1 + v_dim1], ldv, &
		    fact, &dwork[pdw6], n, (ftnlen)12, (ftnlen)12);
	} else if (ltrb) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &b[(*k + 
		    1) * b_dim1 + 1], ldb, &v[(*k + 1) * v_dim1 + 1], ldv, &
		    fact, &dwork[pdw6], n, (ftnlen)12, (ftnlen)9);
	} else if (lcolv) {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[*k + 
		    1 + b_dim1], ldb, &v[*k + 1 + v_dim1], ldv, &fact, &dwork[
		    pdw6], n, (ftnlen)9, (ftnlen)12);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, &v[(*k + 1) * v_dim1 + 1], ldv, &fact, &
		    dwork[pdw6], n, (ftnlen)9, (ftnlen)9);
	}
    } else if (! la1b1) {
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw6], n, (ftnlen)3);
    }

    if (ltrq) {

/*        17) DW7 := DW1*R1 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw7], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[pr1 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        18) DW8 := DW6*R3 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw8], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        19) DW7 := DW7 + DW8 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7 + *n], &c__1);

	if (la1b1) {

/*           NZ15) DW8 := DW9*R2 */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &
		    rs[pr2 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);

/*           NZ16) DW7 := DW7 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);
	}

/*        20) DW8 := DW7*S1 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        21) DW3 := DW3 - DW8 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b367, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &c__1);

/*        22) DW8 := DW7*S3 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps3 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        23) DW5 := DW5 - DW8 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b367, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);

/*        24) DW7 := DW7*S2 */

	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b367, &
		rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
    } else {

/*        17) DW7 := DW6*S3' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw7], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		ps3 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        18) DW8 := DW1*S1' */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw8], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        19) DW7 := DW7 + DW8 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);

	if (la1b1) {

/*           NZ15) DW8 := DW9*S2' */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ16) DW7 := DW7 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);
	}

/*        20) DW8 := DW7*R1' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		pr1 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        21) DW3 := DW3 + DW8 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);

/*        22) DW8 := DW7*R3' */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw7 + *n], &c__1, &dwork[pdw8], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        23) DW5 := DW5 + DW8 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);

/*        24) DW7 := DW7*R2' */

	dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &rs[pr2 * 
		rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
    }

/*     25) A2 := A2 + W2*DW3' */

    if (*m > *k) {
	if (ltra && lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw3], n, &w[*k + 1 + w_dim1], ldw, &c_b22, &a[(*k + 1) * 
		    a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
	} else if (ltra) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw3], n, &w[(*k + 1) * w_dim1 + 1], ldw, &c_b22, &
		    a[(*k + 1) * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
	} else if (lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &w[*k + 
		    1 + w_dim1], ldw, &dwork[pdw3], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)12, (ftnlen)9);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &w[(*k + 1) 
		    * w_dim1 + 1], ldw, &dwork[pdw3], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)9, (ftnlen)9);
	}
    }

/*     26) A2 := A2 + V2*DW5' */

    if (*m > *k) {
	if (ltra && lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw5], n, &v[*k + 1 + v_dim1], ldv, &c_b22, &a[(*k + 1) * 
		    a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
	} else if (ltra) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw5], n, &v[(*k + 1) * v_dim1 + 1], ldv, &c_b22, &
		    a[(*k + 1) * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
	} else if (lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &v[*k + 
		    1 + v_dim1], ldv, &dwork[pdw5], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)12, (ftnlen)9);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &v[(*k + 1) 
		    * v_dim1 + 1], ldv, &dwork[pdw5], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)9, (ftnlen)9);
	}
    }

/*     27) DW4 := DW4 + DW7 */

    i__1 = *n * *k;
    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw4], &c__1);

/*     28) DW3 := DW3*W1' */

    if (lcolw) {
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
    } else {
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
    }

/*     29) DW4 := DW4 + DW3 */

    i__1 = *n * *k;
    daxpy_(&i__1, &c_b22, &dwork[pdw3], &c__1, &dwork[pdw4], &c__1);

/*     30) DW5 := DW5*V1' */

    if (lcolv) {
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
    } else {
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
    }

/*     31) DW4 := DW4 + DW5 */

    i__1 = *n * *k;
    daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*     32) A1 := A1 + DW4' */

    if (la1b1) {
	if (ltra) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		daxpy_(n, &c_b22, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &a[
			i__ * a_dim1 + 1], &c__1);
/* L50: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		daxpy_(k, &c_b22, &dwork[pdw4 + i__ - 1], n, &a[i__ * a_dim1 
			+ 1], &c__1);
/* L60: */
	    }
	}
    } else {
	if (ltra) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &a[i__ * 
			a_dim1 + 1], &c__1);
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(k, &dwork[pdw4 + i__ - 1], n, &a[i__ * a_dim1 + 1], &
			c__1);
/* L80: */
	    }
	}
    }

/*     Update the matrix B. */

    if (ltrq) {

/*        33) DW3 := DW1*T11 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        34) DW4 := DW6*T31 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        35) DW3 := DW3 + DW4 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

	if (la1b1) {

/*           NZ17) DW8 := DW9*T21 */

	    i__1 = *n * (*k - 1);
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    i__1 = *k - 1;
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ18) DW3 := DW3 + DW8 */

	    i__1 = *n * (*k - 1);
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &
		    c__1);
	}

/*        36) DW4 := DW2*S1 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw4], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        37) DW3 := DW3 + DW4 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

/*        38) DW4 := DW1*T12 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt12 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        38) DW5 := DW6*T32 */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        40) DW4 := DW4 + DW5 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4 + *n], &c__1);

	if (la1b1) {

/*           NZ19) DW8 := DW9*T22 */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ20) DW4 := DW4 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
	}

/*        41) DW5 := DW2*S2 */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        42) DW4 := DW4 + DW5 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*        43) DW6 := DW6*T33 */

	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        44) DW1 := DW1*T13 */

	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw1], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        45) DW6 := DW6 + DW1 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw1], &c__1, &dwork[pdw6], &c__1);

	if (la1b1) {

/*           NZ19) DW9 := DW9*T23 */

	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt23 * t_dim1 + 1], ldt, &dwork[pdw9], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ20) DW6 := DW6 + DW9 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw9], &c__1, &dwork[pdw6], &c__1);
	}

/*        46) DW2 := DW2*S3 */

	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps3 * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        45) DW6 := DW6 + DW2 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw2], &c__1, &dwork[pdw6], &c__1);
    } else {

/*        33) DW3 := DW1*T11' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        34) DW4 := DW6*T13' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        35) DW3 := DW3 + DW4 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

	if (la1b1) {

/*           NZ17) DW8 := DW9*T12' */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt12 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ18) DW3 := DW3 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);
	}

/*        36) DW4 := DW2*R1' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		pr1 * rs_dim1 + 1], ldrs, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        37) DW3 := DW3 - DW4 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b367, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

/*        38) DW4 := DW6*T23' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt23 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        39) DW5 := DW1*T21' */

	i__1 = *n * (*k - 1);
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw5], &c__1);
	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        40) DW4 := DW4 + DW5 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

	if (la1b1) {

/*           NZ19) DW8 := DW9*T22' */

	    i__1 = *n * *k;
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ20) DW4 := DW4 + DW8 */

	    i__1 = *n * *k;
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
	}

/*        41) DW5 := DW2*R2' */

	i__1 = *n * *k;
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
	dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &rs[pr2 * 
		rs_dim1 + 1], ldrs, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);

/*        42) DW4 := DW4 - DW5 */

	i__1 = *n * *k;
	daxpy_(&i__1, &c_b367, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*        43) DW6 := DW6*T33' */

	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        44) DW1 := DW1*T31' */

	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw1 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        45) DW6 := DW6 + DW1 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b22, &dwork[pdw1 + *n], &c__1, &dwork[pdw6], &c__1);

	if (la1b1) {

/*           NZ19) DW9 := DW9*T32' */

	    i__1 = *k - 1;
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw9 + *n]
		    , n, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ20) DW6 := DW6 + DW9 */

	    i__1 = *n * (*k - 1);
	    daxpy_(&i__1, &c_b22, &dwork[pdw9 + *n], &c__1, &dwork[pdw6], &
		    c__1);
	}

/*        46) DW2 := DW2*R3' */

	i__1 = *k - 1;
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw2 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        45) DW6 := DW6 - DW2 */

	i__1 = *n * (*k - 1);
	daxpy_(&i__1, &c_b367, &dwork[pdw2 + *n], &c__1, &dwork[pdw6], &c__1);
    }

/*     46) B2 := B2 + W2*DW3' */

    if (*m > *k) {
	if (ltrb && lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw3], n, &w[*k + 1 + w_dim1], ldw, &c_b22, &b[(*k + 1) * 
		    b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
	} else if (ltrb) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw3], n, &w[(*k + 1) * w_dim1 + 1], ldw, &c_b22, &
		    b[(*k + 1) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)12);
	} else if (lcolw) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &w[*k + 
		    1 + w_dim1], ldw, &dwork[pdw3], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)9);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &w[(*k + 1) 
		    * w_dim1 + 1], ldw, &dwork[pdw3], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)9, (ftnlen)9);
	}
    }

/*     47) B2 := B2 + V2*DW6' */

    if (*m > *k) {
	if (ltrb && lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw6], n, &v[*k + 1 + v_dim1], ldv, &c_b22, &b[(*k + 1) * 
		    b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
	} else if (ltrb) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw6], n, &v[(*k + 1) * v_dim1 + 1], ldv, &c_b22, &
		    b[(*k + 1) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)12);
	} else if (lcolv) {
	    i__1 = *m - *k;
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &v[*k + 
		    1 + v_dim1], ldv, &dwork[pdw6], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)9);
	} else {
	    i__1 = *m - *k;
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &v[(*k + 1) 
		    * v_dim1 + 1], ldv, &dwork[pdw6], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)9, (ftnlen)9);
	}
    }

/*     48) DW3 := DW3*W1' */

    if (lcolw) {
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
    } else {
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
    }

/*     49) DW4 := DW4 + DW3 */

    i__1 = *n * *k;
    daxpy_(&i__1, &c_b22, &dwork[pdw3], &c__1, &dwork[pdw4], &c__1);

/*     50) DW6 := DW6*V1' */

    if (lcolv) {
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
    } else {
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
    }

/*     51) DW4 := DW4 + DW6 */

    i__1 = *n * *k;
    daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);

/*     52) B1 := B1 + DW4' */

    if (la1b1) {
	if (ltrb) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		daxpy_(n, &c_b22, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &b[
			i__ * b_dim1 + 1], &c__1);
/* L90: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		daxpy_(k, &c_b22, &dwork[pdw4 + i__ - 1], n, &b[i__ * b_dim1 
			+ 1], &c__1);
/* L100: */
	    }
	}
    } else {
	if (ltrb) {
	    i__1 = *k;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &b[i__ * 
			b_dim1 + 1], &c__1);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(k, &dwork[pdw4 + i__ - 1], n, &b[i__ * b_dim1 + 1], &
			c__1);
/* L120: */
	    }
	}
    }

    return 0;
/* *** Last line of MB04QC *** */
} /* mb04qc_ */

