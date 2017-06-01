/* MB03PY.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;

/* Subroutine */ int mb03py_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *rcond, doublereal *svlmax, integer *rank, doublereal 
	*sval, integer *jpvt, doublereal *tau, doublereal *dwork, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1, c2, s1, s2, aii;
    static integer mki, nki, pvt;
    static doublereal smin, temp, smax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static integer itemp, ismin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismax;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal sminpr, smaxpr;


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

/*     To compute a rank-revealing RQ factorization of a real general */
/*     M-by-N matrix  A,  which may be rank-deficient, and estimate its */
/*     effective rank using incremental condition estimation. */

/*     The routine uses a truncated RQ factorization with row pivoting: */
/*                                   [ R11 R12 ] */
/*        P * A = R * Q,  where  R = [         ], */
/*                                   [  0  R22 ] */
/*     with R22 defined as the largest trailing upper triangular */
/*     submatrix whose estimated condition number is less than 1/RCOND. */
/*     The order of R22, RANK, is the effective rank of A.  Condition */
/*     estimation is performed during the RQ factorization process. */
/*     Matrix R11 is full (but of small norm), or empty. */

/*     MB03PY  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the upper triangle of the subarray */
/*             A(M-RANK+1:M,N-RANK+1:N) contains the RANK-by-RANK upper */
/*             triangular matrix R22;  the remaining elements in the last */
/*             RANK  rows, with the array TAU, represent the orthogonal */
/*             matrix Q as a product of  RANK  elementary reflectors */
/*             (see METHOD).  The first  M-RANK  rows contain the result */
/*             of the RQ factorization process used. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest trailing triangular */
/*             submatrix R22 in the RQ factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             0 <= RCOND <= 1. */
/*             NOTE that when SVLMAX > 0, the estimated rank could be */
/*             less than that defined above (see SVLMAX). */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix B, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of B */
/*             (for instance, the Frobenius norm of B).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     RANK    (output) INTEGER */
/*             The effective (estimated) rank of A, i.e., the order of */
/*             the submatrix R22. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R: */
/*             SVAL(1): largest singular value of */
/*                      R(M-RANK+1:M,N-RANK+1:N); */
/*             SVAL(2): smallest singular value of */
/*                      R(M-RANK+1:M,N-RANK+1:N); */
/*             SVAL(3): smallest singular value of R(M-RANK:M,N-RANK:N), */
/*                      if RANK < MIN( M, N ), or of */
/*                      R(M-RANK+1:M,N-RANK+1:N), otherwise. */
/*             If the triangular factorization is a rank-revealing one */
/*             (which will be the case if the trailing rows were well- */
/*             conditioned), then SVAL(1) will also be an estimate for */
/*             the largest singular value of A, and SVAL(2) and SVAL(3) */
/*             will be estimates for the RANK-th and (RANK+1)-st singular */
/*             values of A, respectively. */
/*             By examining these values, one can confirm that the rank */
/*             is well defined with respect to the chosen value of RCOND. */
/*             The ratio SVAL(1)/SVAL(2) is an estimate of the condition */
/*             number of R(M-RANK+1:M,N-RANK+1:N). */

/*     JPVT    (output) INTEGER array, dimension ( M ) */
/*             If JPVT(i) = k, then the i-th row of P*A was the k-th row */
/*             of A. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             The trailing  RANK  elements of TAU contain the scalar */
/*             factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( 3*M-1 ) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated RQ factorization with row */
/*     pivoting of A,  P * A = R * Q,  with  R  defined above, and, */
/*     during this process, finds the largest trailing submatrix whose */
/*     estimated condition number is less than 1/RCOND, taking the */
/*     possible positive value of SVLMAX into account.  This is performed */
/*     using an adaptation of the LAPACK incremental condition estimation */
/*     scheme and a slightly modified rank decision test.  The */
/*     factorization process stops when  RANK  has been determined. */

/*     The matrix Q is represented as a product of elementary reflectors */

/*        Q = H(k-rank+1) H(k-rank+2) . . . H(k), where k = min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit */
/*     in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */

/*     The matrix P is represented in jpvt as follows: If */
/*        jpvt(j) = i */
/*     then the jth row of P is the ith canonical unit vector. */

/*     REFERENCES */

/*     [1] Bischof, C.H. and P. Tang. */
/*         Generalizing Incremental Condition Estimation. */
/*         LAPACK Working Notes 32, Mathematics and Computer Science */
/*         Division, Argonne National Laboratory, UT, CS-91-132, */
/*         May 1991. */

/*     [2] Bischof, C.H. and P. Tang. */
/*         Robust Incremental Condition Estimation. */
/*         LAPACK Working Notes 33, Mathematics and Computer Science */
/*         Division, Argonne National Laboratory, UT, CS-91-133, */
/*         May 1991. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001, */
/*     Jan. 2009. */

/*     KEYWORDS */

/*     Eigenvalue problem, matrix operations, orthogonal transformation, */
/*     singular values. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --sval;
    --jpvt;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*rcond < 0. || *rcond > 1.) {
	*info = -5;
    } else if (*svlmax < 0.) {
	*info = -6;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03PY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    k = min(*m,*n);
    if (k == 0) {
	*rank = 0;
	sval[1] = 0.;
	sval[2] = 0.;
	sval[3] = 0.;
	return 0;
    }

    ismin = *m;
    ismax = ismin + *m;
    jwork = ismax + 1;

/*     Initialize partial row norms and pivoting vector. The first m */
/*     elements of DWORK store the exact row norms. The already used */
/*     trailing part is then overwritten by the condition estimator. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[i__] = dnrm2_(n, &a[i__ + a_dim1], lda);
	dwork[*m + i__] = dwork[i__];
	jpvt[i__] = i__;
/* L10: */
    }

/*     Compute factorization and determine RANK using incremental */
/*     condition estimation. */

    *rank = 0;

L20:
    if (*rank < k) {
	i__ = k - *rank;

/*        Determine ith pivot row and swap if necessary. */

	mki = *m - *rank;
	nki = *n - *rank;
	pvt = idamax_(&mki, &dwork[1], &c__1);

	if (pvt != mki) {
	    dswap_(n, &a[pvt + a_dim1], lda, &a[mki + a_dim1], lda);
	    itemp = jpvt[pvt];
	    jpvt[pvt] = jpvt[mki];
	    jpvt[mki] = itemp;
	    dwork[pvt] = dwork[mki];
	    dwork[*m + pvt] = dwork[*m + mki];
	}

	if (nki > 1) {

/*           Save A(m-k+i,n-k+i) and generate elementary reflector H(i) */
/*           to annihilate A(m-k+i,1:n-k+i-1), k = min(m,n). */

	    aii = a[mki + nki * a_dim1];
	    dlarfg_(&nki, &a[mki + nki * a_dim1], &a[mki + a_dim1], lda, &tau[
		    i__]);
	}

	if (*rank == 0) {

/*           Initialize; exit if matrix is zero (RANK = 0). */

	    smax = (d__1 = a[*m + *n * a_dim1], abs(d__1));
	    if (smax == 0.) {
		sval[1] = 0.;
		sval[2] = 0.;
		sval[3] = 0.;
		return 0;
	    }
	    smin = smax;
	    smaxpr = smax;
	    sminpr = smin;
	    c1 = 1.;
	    c2 = 1.;
	} else {

/*           One step of incremental condition estimation. */

	    dcopy_(rank, &a[mki + (nki + 1) * a_dim1], lda, &dwork[jwork], &
		    c__1);
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &dwork[jwork], &a[mki 
		    + nki * a_dim1], &sminpr, &s1, &c1);
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &dwork[jwork], &a[mki 
		    + nki * a_dim1], &smaxpr, &s2, &c2);
	}

	if (*svlmax * *rcond <= smaxpr) {
	    if (*svlmax * *rcond <= sminpr) {
		if (smaxpr * *rcond <= sminpr) {

		    if (mki > 1) {

/*                    Continue factorization, as rank is at least RANK. */
/*                    Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right. */

			aii = a[mki + nki * a_dim1];
			a[mki + nki * a_dim1] = 1.;
			i__1 = mki - 1;
			dlarf_("Right", &i__1, &nki, &a[mki + a_dim1], lda, &
				tau[i__], &a[a_offset], lda, &dwork[jwork], (
				ftnlen)5);
			a[mki + nki * a_dim1] = aii;

/*                    Update partial row norms. */

			i__1 = mki - 1;
			for (j = 1; j <= i__1; ++j) {
			    if (dwork[j] != 0.) {
/* Computing 2nd power */
				d__2 = (d__1 = a[j + nki * a_dim1], abs(d__1))
					 / dwork[j];
				temp = 1. - d__2 * d__2;
				temp = max(temp,0.);
/* Computing 2nd power */
				d__1 = dwork[j] / dwork[*m + j];
				temp2 = temp * .05 * (d__1 * d__1) + 1.;
				if (temp2 == 1.) {
				    i__2 = nki - 1;
				    dwork[j] = dnrm2_(&i__2, &a[j + a_dim1], 
					    lda);
				    dwork[*m + j] = dwork[j];
				} else {
				    dwork[j] *= sqrt(temp);
				}
			    }
/* L30: */
			}

		    }

		    i__1 = *rank;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 1];
			dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 1];
/* L40: */
		    }

		    if (*rank > 0) {
			--ismin;
			--ismax;
		    }
		    dwork[ismin] = c1;
		    dwork[ismax] = c2;
		    smin = sminpr;
		    smax = smaxpr;
		    ++(*rank);
		    goto L20;
		}
	    }
	}
    }

/*     Restore the changed part of the (M-RANK)-th row and set SVAL. */

    if (*rank < k && nki > 1) {
	i__1 = nki - 1;
	d__1 = -a[mki + nki * a_dim1] * tau[i__];
	dscal_(&i__1, &d__1, &a[mki + a_dim1], lda);
	a[mki + nki * a_dim1] = aii;
    }
    sval[1] = smax;
    sval[2] = smin;
    sval[3] = sminpr;

    return 0;
/* *** Last line of MB03PY *** */
} /* mb03py_ */

