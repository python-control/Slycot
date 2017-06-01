/* MB3OYZ.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb3oyz_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *rcond, doublereal *svlmax, integer *rank, 
	doublereal *sval, integer *jpvt, doublecomplex *tau, doublereal *
	dwork, doublecomplex *zwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublecomplex c1, c2, s1, s2;
    static integer mn;
    static doublecomplex aii;
    static integer pvt;
    static doublereal smin, temp, smax, temp2;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer itemp, ismin;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static integer ismax;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlaic1_(integer *, integer *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *);
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

/*     To compute a rank-revealing QR factorization of a complex general */
/*     M-by-N matrix  A,  which may be rank-deficient, and estimate its */
/*     effective rank using incremental condition estimation. */

/*     The routine uses a truncated QR factorization with column pivoting */
/*                                   [ R11 R12 ] */
/*        A * P = Q * R,  where  R = [         ], */
/*                                   [  0  R22 ] */
/*     with R11 defined as the largest leading upper triangular submatrix */
/*     whose estimated condition number is less than 1/RCOND.  The order */
/*     of R11, RANK, is the effective rank of A.  Condition estimation is */
/*     performed during the QR factorization process.  Matrix R22 is full */
/*     (but of small norm), or empty. */

/*     MB3OYZ  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension ( LDA, N ) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the leading RANK-by-RANK upper triangular part */
/*             of A contains the triangular factor R11, and the elements */
/*             below the diagonal in the first  RANK  columns, with the */
/*             array TAU, represent the unitary matrix Q as a product */
/*             of  RANK  elementary reflectors. */
/*             The remaining  N-RANK  columns contain the result of the */
/*             QR factorization process used. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest leading triangular */
/*             submatrix R11 in the QR factorization with pivoting of A, */
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
/*             the submatrix R11. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R: */
/*             SVAL(1): largest singular value of R(1:RANK,1:RANK); */
/*             SVAL(2): smallest singular value of R(1:RANK,1:RANK); */
/*             SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1), */
/*                      if RANK < MIN( M, N ), or of R(1:RANK,1:RANK), */
/*                      otherwise. */
/*             If the triangular factorization is a rank-revealing one */
/*             (which will be the case if the leading columns were well- */
/*             conditioned), then SVAL(1) will also be an estimate for */
/*             the largest singular value of A, and SVAL(2) and SVAL(3) */
/*             will be estimates for the RANK-th and (RANK+1)-st singular */
/*             values of A, respectively. */
/*             By examining these values, one can confirm that the rank */
/*             is well defined with respect to the chosen value of RCOND. */
/*             The ratio SVAL(1)/SVAL(2) is an estimate of the condition */
/*             number of R(1:RANK,1:RANK). */

/*     JPVT    (output) INTEGER array, dimension ( N ) */
/*             If JPVT(i) = k, then the i-th column of A*P was the k-th */
/*             column of A. */

/*     TAU     (output) COMPLEX*16 array, dimension ( MIN( M, N ) ) */
/*             The leading  RANK  elements of TAU contain the scalar */
/*             factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( 2*N ) */

/*     ZWORK   COMPLEX*16 array, dimension ( 3*N-1 ) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated QR factorization with column */
/*     pivoting of A,  A * P = Q * R,  with  R  defined above, and, */
/*     during this process, finds the largest leading submatrix whose */
/*     estimated condition number is less than 1/RCOND, taking the */
/*     possible positive value of SVLMAX into account.  This is performed */
/*     using the LAPACK incremental condition estimation scheme and a */
/*     slightly modified rank decision test.  The factorization process */
/*     stops when  RANK  has been determined. */

/*     The matrix Q is represented as a product of elementary reflectors */

/*        Q = H(1) H(2) . . . H(k), where k = rank <= min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a complex scalar, and v is a complex vector with */
/*     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in */
/*     A(i+1:m,i), and tau in TAU(i). */

/*     The matrix P is represented in jpvt as follows: If */
/*        jpvt(j) = i */
/*     then the jth column of P is the ith canonical unit vector. */

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
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue problem, matrix operations, unitary transformation, */
/*     singular values. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
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
    --zwork;

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
	xerbla_("MB3OYZ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    mn = min(*m,*n);
    if (mn == 0) {
	*rank = 0;
	sval[1] = 0.;
	sval[2] = 0.;
	sval[3] = 0.;
	return 0;
    }

    ismin = 1;
    ismax = ismin + *n;

/*     Initialize partial column norms and pivoting vector. The first n */
/*     elements of DWORK store the exact column norms. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[i__] = dznrm2_(m, &a[i__ * a_dim1 + 1], &c__1);
	dwork[*n + i__] = dwork[i__];
	jpvt[i__] = i__;
/* L10: */
    }

/*     Compute factorization and determine RANK using incremental */
/*     condition estimation. */

    *rank = 0;

L20:
    if (*rank < mn) {
	i__ = *rank + 1;

/*        Determine ith pivot column and swap if necessary. */

	i__1 = *n - i__ + 1;
	pvt = i__ - 1 + idamax_(&i__1, &dwork[i__], &c__1);

	if (pvt != i__) {
	    zswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
	    itemp = jpvt[pvt];
	    jpvt[pvt] = jpvt[i__];
	    jpvt[i__] = itemp;
	    dwork[pvt] = dwork[i__];
	    dwork[*n + pvt] = dwork[*n + i__];
	}

/*        Save A(I,I) and generate elementary reflector H(i) */
/*        such that H(i)'*[A(i,i);*] = [*;0]. */

	if (i__ < *m) {
	    i__1 = i__ + i__ * a_dim1;
	    aii.r = a[i__1].r, aii.i = a[i__1].i;
	    i__1 = *m - i__ + 1;
	    zlarfg_(&i__1, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1],
		     &c__1, &tau[i__]);
	} else {
	    i__1 = *m;
	    tau[i__1].r = 0., tau[i__1].i = 0.;
	}

	if (*rank == 0) {

/*           Initialize; exit if matrix is zero (RANK = 0). */

	    smax = z_abs(&a[a_dim1 + 1]);
	    if (smax == 0.) {
		sval[1] = 0.;
		sval[2] = 0.;
		sval[3] = 0.;
		return 0;
	    }
	    smin = smax;
	    smaxpr = smax;
	    sminpr = smin;
	    c1.r = 1., c1.i = 0.;
	    c2.r = 1., c2.i = 0.;
	} else {

/*           One step of incremental condition estimation. */

	    zlaic1_(&c__2, rank, &zwork[ismin], &smin, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &sminpr, &s1, &c1);
	    zlaic1_(&c__1, rank, &zwork[ismax], &smax, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &smaxpr, &s2, &c2);
	}

	if (*svlmax * *rcond <= smaxpr) {
	    if (*svlmax * *rcond <= sminpr) {
		if (smaxpr * *rcond <= sminpr) {

/*                 Continue factorization, as rank is at least RANK. */

		    if (i__ < *n) {

/*                    Apply H(i)' to A(i:m,i+1:n) from the left. */

			i__1 = i__ + i__ * a_dim1;
			aii.r = a[i__1].r, aii.i = a[i__1].i;
			i__1 = i__ + i__ * a_dim1;
			a[i__1].r = 1., a[i__1].i = 0.;
			i__1 = *m - i__ + 1;
			i__2 = *n - i__;
			d_cnjg(&z__1, &tau[i__]);
			zlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &
				c__1, &z__1, &a[i__ + (i__ + 1) * a_dim1], 
				lda, &zwork[(*n << 1) + 1], (ftnlen)4);
			i__1 = i__ + i__ * a_dim1;
			a[i__1].r = aii.r, a[i__1].i = aii.i;
		    }

/*                 Update partial column norms. */

		    i__1 = *n;
		    for (j = i__ + 1; j <= i__1; ++j) {
			if (dwork[j] != 0.) {
/* Computing 2nd power */
			    d__1 = z_abs(&a[i__ + j * a_dim1]) / dwork[j];
			    temp = 1. - d__1 * d__1;
			    temp = max(temp,0.);
/* Computing 2nd power */
			    d__1 = dwork[j] / dwork[*n + j];
			    temp2 = temp * .05 * (d__1 * d__1) + 1.;
			    if (temp2 == 1.) {
				if (*m - i__ > 0) {
				    i__2 = *m - i__;
				    dwork[j] = dznrm2_(&i__2, &a[i__ + 1 + j *
					     a_dim1], &c__1);
				    dwork[*n + j] = dwork[j];
				} else {
				    dwork[j] = 0.;
				    dwork[*n + j] = 0.;
				}
			    } else {
				dwork[j] *= sqrt(temp);
			    }
			}
/* L30: */
		    }

		    i__1 = *rank;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			i__2 = ismin + i__ - 1;
			i__3 = ismin + i__ - 1;
			z__1.r = s1.r * zwork[i__3].r - s1.i * zwork[i__3].i, 
				z__1.i = s1.r * zwork[i__3].i + s1.i * zwork[
				i__3].r;
			zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
			i__2 = ismax + i__ - 1;
			i__3 = ismax + i__ - 1;
			z__1.r = s2.r * zwork[i__3].r - s2.i * zwork[i__3].i, 
				z__1.i = s2.r * zwork[i__3].i + s2.i * zwork[
				i__3].r;
			zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L40: */
		    }

		    i__1 = ismin + *rank;
		    zwork[i__1].r = c1.r, zwork[i__1].i = c1.i;
		    i__1 = ismax + *rank;
		    zwork[i__1].r = c2.r, zwork[i__1].i = c2.i;
		    smin = sminpr;
		    smax = smaxpr;
		    ++(*rank);
		    goto L20;
		}
	    }
	}
    }

/*     Restore the changed part of the (RANK+1)-th column and set SVAL. */

    if (*rank < *n) {
	if (i__ < *m) {
	    i__1 = *m - i__;
	    i__2 = i__ + i__ * a_dim1;
	    z__2.r = -a[i__2].r, z__2.i = -a[i__2].i;
	    i__3 = i__;
	    z__1.r = z__2.r * tau[i__3].r - z__2.i * tau[i__3].i, z__1.i = 
		    z__2.r * tau[i__3].i + z__2.i * tau[i__3].r;
	    zscal_(&i__1, &z__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
	    i__1 = i__ + i__ * a_dim1;
	    a[i__1].r = aii.r, a[i__1].i = aii.i;
	}
    }
    if (*rank == 0) {
	smin = 0.;
	sminpr = 0.;
    }
    sval[1] = smax;
    sval[2] = smin;
    sval[3] = sminpr;

    return 0;
/* *** Last line of MB3OYZ *** */
} /* mb3oyz_ */

