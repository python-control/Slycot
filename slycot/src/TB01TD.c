/* TB01TD.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;

/* Subroutine */ int tb01td_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, integer *low, integer *
	igh, doublereal *scstat, doublereal *scin, doublereal *scout, 
	doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, kold, knew;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), tb01ty_(integer *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *), dgebal_(char *, integer *, doublereal *, integer *,
	     integer *, integer *, doublereal *, integer *, ftnlen);
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal acnorm, arnorm;


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

/*     To reduce a given state-space representation (A,B,C,D) to */
/*     balanced form by means of state permutations and state, input and */
/*     output scalings. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, i.e. the */
/*             order of the original state dynamics matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced state dynamics matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the balanced input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the balanced state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the original direct transmission matrix D. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the scaled direct transmission matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     LOW     (output) INTEGER */
/*             The index of the lower end of the balanced submatrix of A. */

/*     IGH     (output) INTEGER */
/*             The index of the upper end of the balanced submatrix of A. */

/*     SCSTAT  (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the information defining the */
/*             similarity transformations used to permute and balance */
/*             the state dynamics matrix A, as returned from the LAPACK */
/*             library routine DGEBAL. */

/*     SCIN    (output) DOUBLE PRECISION array, dimension (M) */
/*             Contains the scalars used to scale the system inputs so */
/*             that the columns of the final matrix B have norms roughly */
/*             equal to the column sums of the balanced matrix A */
/*             (see FURTHER COMMENTS). */
/*             The j-th input of the balanced state-space representation */
/*             is SCIN(j)*(j-th column of the permuted and balanced */
/*             input/state matrix B). */

/*     SCOUT   (output) DOUBLE PRECISION array, dimension (P) */
/*             Contains the scalars used to scale the system outputs so */
/*             that the rows of the final matrix C have norms roughly */
/*             equal to the row sum of the balanced matrix A. */
/*             The i-th output of the balanced state-space representation */
/*             is SCOUT(i)*(i-th row of the permuted and balanced */
/*             state/ouput matrix C). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Similarity transformations are used to permute the system states */
/*     and balance the corresponding row and column sum norms of a */
/*     submatrix of the state dynamics matrix A. These operations are */
/*     also applied to the input/state matrix B and the system inputs */
/*     are then scaled (see parameter SCIN) so that the columns of the */
/*     final matrix B have norms roughly equal to the column sum norm of */
/*     the balanced matrix A (see FURTHER COMMENTS). */
/*     The above operations are also applied to the matrix C, and the */
/*     system outputs are then scaled (see parameter SCOUT) so that the */
/*     rows of the final matrix C have norms roughly equal to the row sum */
/*     norm of the balanced matrix A (see FURTHER COMMENTS). */
/*     Finally, the (I,J)-th element of the direct transmission matrix D */
/*     is scaled as */
/*          D(I,J) = D(I,J)*(1.0/SCIN(J))*SCOUT(I), where I = 1,2,...,P */
/*     and J = 1,2,...,M. */

/*     Scaling performed to balance the row/column sum norms is by */
/*     integer powers of the machine base so as to avoid introducing */
/*     rounding errors. */

/*     REFERENCES */

/*     [1] Wilkinson, J.H. and Reinsch, C. */
/*         Handbook for Automatic Computation, (Vol II, Linear Algebra). */
/*         Springer-Verlag, 1971, (contribution II/11). */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     The columns (rows) of the final matrix B (matrix C) have norms */
/*     'roughly' equal to the column (row) sum norm of the balanced */
/*     matrix A, i.e. */
/*        size/BASE < abssum <= size */
/*     where */
/*        BASE   = the base of the arithmetic used on the computer, which */
/*                 can be obtained from the LAPACK Library routine */
/*                 DLAMCH; */

/*        size   = column or row sum norm of the balanced matrix A; */
/*        abssum = column sum norm of the balanced matrix B or row sum */
/*                 norm of the balanced matrix C. */

/*     The routine is BASE dependent. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01HD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, October 1982. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balanced form, orthogonal transformation, similarity */
/*     transformation, state-space model, state-space representation. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --scstat;
    --scin;
    --scout;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldc < max(1,*p)) {
	*info = -9;
    } else if (*ldd < max(1,*p)) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB01TD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*n,*m);
    if (max(i__1,*p) == 0) {
	*low = 1;
	*igh = *n;
	return 0;
    }

/*     Permute states, and balance a submatrix of A. */

    dgebal_("Both", n, &a[a_offset], lda, low, igh, &scstat[1], info, (ftnlen)
	    4);

/*     Use the information in SCSTAT on state scalings and reorderings */
/*     to transform B and C. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kold = k;
	if (*low > kold || kold > *igh) {
	    if (kold < *low) {
		kold = *low - kold;
	    }
	    knew = (integer) scstat[kold];
	    if (knew != kold) {

/*              Exchange rows KOLD and KNEW of B. */

		dswap_(m, &b[kold + b_dim1], ldb, &b[knew + b_dim1], ldb);

/*              Exchange columns KOLD and KNEW of C. */

		dswap_(p, &c__[kold * c_dim1 + 1], &c__1, &c__[knew * c_dim1 
			+ 1], &c__1);
	    }
	}
/* L10: */
    }

    if (*igh != *low) {

	i__1 = *igh;
	for (k = *low; k <= i__1; ++k) {
	    scale = scstat[k];

/*           Scale the K-th row of permuted B. */

	    d__1 = 1. / scale;
	    dscal_(m, &d__1, &b[k + b_dim1], ldb);

/*           Scale the K-th column of permuted C. */

	    dscal_(p, &scale, &c__[k * c_dim1 + 1], &c__1);
/* L20: */
	}

    }

/*     Calculate the column and row sum norms of A. */

    acnorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)6);
    arnorm = dlange_("I-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)6);

/*     Scale the columns of B (i.e. inputs) to have norms roughly ACNORM. */

    tb01ty_(&c__1, &c__0, &c__0, n, m, &acnorm, &b[b_offset], ldb, &scin[1]);

/*     Scale the rows of C (i.e. outputs) to have norms roughly ARNORM. */

    tb01ty_(&c__0, &c__0, &c__0, p, n, &arnorm, &c__[c_offset], ldc, &scout[1]
	    );

/*     Finally, apply these input and output scalings to D and set SCIN. */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	scale = scin[j];

	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    d__[i__ + j * d_dim1] *= scale * scout[i__];
/* L30: */
	}

	scin[j] = 1. / scale;
/* L40: */
    }

    return 0;
/* *** Last line of TB01TD *** */
} /* tb01td_ */

