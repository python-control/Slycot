/* MB03OD.f -- translated by f2c (version 20100827).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03od_(char *jobqr, integer *m, integer *n, doublereal *
	a, integer *lda, integer *jpvt, doublereal *rcond, doublereal *svlmax,
	 doublereal *tau, integer *rank, doublereal *sval, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobqr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal c1, c2, s1, s2;
    static integer mn;
    static doublereal smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ismin, ismax;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dgeqp3_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical ljobqr;
    static integer minwrk;
    static doublereal sminpr;
    static integer maxwrk;
    static doublereal smaxpr;


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

/*     To compute (optionally) a rank-revealing QR factorization of a */
/*     real general M-by-N matrix  A,  which may be rank-deficient, */
/*     and estimate its effective rank using incremental condition */
/*     estimation. */

/*     The routine uses a QR factorization with column pivoting: */
/*        A * P = Q * R,  where  R = [ R11 R12 ], */
/*                                   [  0  R22 ] */
/*     with R11 defined as the largest leading submatrix whose estimated */
/*     condition number is less than 1/RCOND.  The order of R11, RANK, */
/*     is the effective rank of A. */

/*     MB03OD  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBQR   CHARACTER*1 */
/*             = 'Q':  Perform a QR factorization with column pivoting; */
/*             = 'N':  Do not perform the QR factorization (but assume */
/*                     that it has been done outside). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry with JOBQR = 'Q', the leading M by N part of this */
/*             array must contain the given matrix A. */
/*             On exit with JOBQR = 'Q', the leading min(M,N) by N upper */
/*             triangular part of A contains the triangular factor R, */
/*             and the elements below the diagonal, with the array TAU, */
/*             represent the orthogonal matrix Q as a product of */
/*             min(M,N) elementary reflectors. */
/*             On entry and on exit with JOBQR = 'N', the leading */
/*             min(M,N) by N upper triangular part of A contains the */
/*             triangular factor R, as determined by the QR factorization */
/*             with pivoting.  The elements below the diagonal of A are */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     JPVT    (input/output) INTEGER array, dimension ( N ) */
/*             On entry with JOBQR = 'Q', if JPVT(i) <> 0, the i-th */
/*             column of A is an initial column, otherwise it is a free */
/*             column. Before the QR factorization of A, all initial */
/*             columns are permuted to the leading positions; only the */
/*             remaining free columns are moved as a result of column */
/*             pivoting during the factorization.  For rank determination */
/*             it is preferable that all columns be free. */
/*             On exit with JOBQR = 'Q', if JPVT(i) = k, then the i-th */
/*             column of A*P was the k-th column of A. */
/*             Array JPVT is not referenced when JOBQR = 'N'. */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest leading triangular */
/*             submatrix R11 in the QR factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             RCOND >= 0. */
/*             NOTE that when SVLMAX > 0, the estimated rank could be */
/*             less than that defined above (see SVLMAX). */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix B, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of B */
/*             (for instance, the Frobenius norm of B).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             On exit with JOBQR = 'Q', the leading min(M,N) elements of */
/*             TAU contain the scalar factors of the elementary */
/*             reflectors. */
/*             Array TAU is not referenced when JOBQR = 'N'. */

/*     RANK    (output) INTEGER */
/*             The effective (estimated) rank of A, i.e. the order of */
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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( LDWORK ) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 3*N + 1,                 if JOBQR = 'Q'; */
/*             LDWORK >= max( 1, 2*min( M, N ) ), if JOBQR = 'N'. */
/*             For good performance when JOBQR = 'Q', LDWORK should be */
/*             larger. Specifically, LDWORK >= 2*N + ( N + 1 )*NB, where */
/*             NB is the optimal block size for the LAPACK Library */
/*             routine DGEQP3. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes or uses a QR factorization with column */
/*     pivoting of A,  A * P = Q * R,  with  R  defined above, and then */
/*     finds the largest leading submatrix whose estimated condition */
/*     number is less than 1/RCOND, taking the possible positive value of */
/*     SVLMAX into account.  This is performed using the LAPACK */
/*     incremental condition estimation scheme and a slightly modified */
/*     rank decision test. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --sval;
    --dwork;

    /* Function Body */
    ljobqr = lsame_(jobqr, "Q", (ftnlen)1, (ftnlen)1);
    mn = min(*m,*n);
    ismin = 1;
    ismax = mn + 1;
    if (ljobqr) {
	minwrk = *n * 3 + 1;
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = mn << 1;
	minwrk = max(i__1,i__2);
    }
    maxwrk = minwrk;

/*     Test the input scalar arguments. */

    *info = 0;
    if (! ljobqr && ! lsame_(jobqr, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*rcond < 0.) {
	*info = -7;
    } else if (*svlmax < 0.) {
	*info = -8;
    } else if (*ldwork < minwrk) {
	*info = -13;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (mn == 0) {
	*rank = 0;
	sval[1] = 0.;
	sval[2] = 0.;
	sval[3] = 0.;
	dwork[1] = 1.;
	return 0;
    }

    if (ljobqr) {

/*        Compute QR factorization with column pivoting of A: */
/*           A * P = Q * R */
/*        Workspace need   3*N + 1; */
/*                  prefer 2*N + (N+1)*NB. */
/*        Details of Householder rotations stored in TAU. */

	dgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &tau[1], &dwork[1], ldwork,
		 info);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[1];
	maxwrk = max(i__1,i__2);
    }

/*     Determine RANK using incremental condition estimation */

    dwork[ismin] = 1.;
    dwork[ismax] = 1.;
    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
    smin = smax;
    if (smax == 0. || *svlmax * *rcond > smax) {
	*rank = 0;
	sval[1] = smax;
	sval[2] = 0.;
	sval[3] = 0.;
    } else {
	*rank = 1;
	sminpr = smin;

L10:
	if (*rank < mn) {
	    i__ = *rank + 1;
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &sminpr, &s1, &c1);
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

	    if (*svlmax * *rcond <= smaxpr) {
		if (*svlmax * *rcond <= sminpr) {
		    if (smaxpr * *rcond <= sminpr) {
			i__1 = *rank;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 
				    1];
			    dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 
				    1];
/* L20: */
			}
			dwork[ismin + *rank] = c1;
			dwork[ismax + *rank] = c2;
			smin = sminpr;
			smax = smaxpr;
			++(*rank);
			goto L10;
		    }
		}
	    }
	}
	sval[1] = smax;
	sval[2] = smin;
	sval[3] = sminpr;
    }

    dwork[1] = (doublereal) maxwrk;
    return 0;
/* *** Last line of MB03OD *** */
} /* mb03od_ */

