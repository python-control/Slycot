/* FB01QD.f -- translated by f2c (version 20100827).
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

static doublereal c_b13 = 1.;

/* Subroutine */ int fb01qd_(char *jobk, char *multbq, integer *n, integer *m,
	 integer *p, doublereal *s, integer *lds, doublereal *a, integer *lda,
	 doublereal *b, integer *ldb, doublereal *q, integer *ldq, doublereal 
	*c__, integer *ldc, doublereal *r__, integer *ldr, doublereal *k, 
	integer *ldk, doublereal *tol, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobk_len, ftnlen multbq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, k_dim1, 
	    k_offset, q_dim1, q_offset, r_dim1, r_offset, s_dim1, s_offset, 
	    i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer n1, i12, pn, itau;
    extern /* Subroutine */ int mb04ld_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), mb02od_(char *, char *, char *, char *, char *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical ljobk;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical lmultb;
    static integer wrkopt;


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

/*     To calculate a combined measurement and time update of one */
/*     iteration of the time-varying Kalman filter. This update is given */
/*     for the square root covariance filter, using dense matrices. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBK    CHARACTER*1 */
/*             Indicates whether the user wishes to compute the Kalman */
/*             filter gain matrix K  as follows: */
/*                                 i */
/*             = 'K':  K  is computed and stored in array K; */
/*                      i */
/*             = 'N':  K  is not required. */
/*                      i */

/*     MULTBQ  CHARACTER*1                    1/2 */
/*             Indicates how matrices B  and Q    are to be passed to */
/*                                     i      i */
/*             the routine as follows: */
/*             = 'P':  Array Q is not used and the array B must contain */
/*                                    1/2 */
/*                     the product B Q   ; */
/*                                  i i */
/*             = 'N':  Arrays B and Q must contain the matrices as */
/*                     described below. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices S    and A .  N >= 0. */
/*                       i-1      i */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e., the order of the matrix */
/*              1/2 */
/*             Q   .  M >= 0. */
/*              i */

/*     P       (input) INTEGER */
/*             The actual output dimension, i.e., the order of the matrix */
/*              1/2 */
/*             R   .  P >= 0. */
/*              i */

/*     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             On entry, the leading N-by-N lower triangular part of this */
/*             array must contain S   , the square root (left Cholesky */
/*                                 i-1 */
/*             factor) of the state covariance matrix at instant (i-1). */
/*             On exit, the leading N-by-N lower triangular part of this */
/*             array contains S , the square root (left Cholesky factor) */
/*                             i */
/*             of the state covariance matrix at instant i. */
/*             The strict upper triangular part of this array is not */
/*             referenced. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain A , */
/*                                                                 i */
/*             the state transition matrix of the discrete system at */
/*             instant i. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain B , */
/*                                                        1/2      i */
/*             the input weight matrix (or the product B Q    if */
/*                                                      i i */
/*             MULTBQ = 'P') of the discrete system at instant i. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             If MULTBQ = 'N', then the leading M-by-M lower triangular */
/*                                              1/2 */
/*             part of this array must contain Q   , the square root */
/*                                              i */
/*             (left Cholesky factor) of the input (process) noise */
/*             covariance matrix at instant i. */
/*             The strict upper triangular part of this array is not */
/*             referenced. */
/*             If MULTBQ = 'P', Q is not referenced and can be supplied */
/*             as a dummy array (i.e., set parameter LDQ = 1 and declare */
/*             this array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,M) if MULTBQ = 'N'; */
/*             LDQ >= 1        if MULTBQ = 'P'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain C , the */
/*                                                                 i */
/*             output weight matrix of the discrete system at instant i. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,P) */
/*             On entry, the leading P-by-P lower triangular part of this */
/*                                 1/2 */
/*             array must contain R   , the square root (left Cholesky */
/*                                 i */
/*             factor) of the output (measurement) noise covariance */
/*             matrix at instant i. */
/*             On exit, the leading P-by-P lower triangular part of this */
/*                                    1/2 */
/*             array contains (RINOV )   , the square root (left Cholesky */
/*                                  i */
/*             factor) of the covariance matrix of the innovations at */
/*             instant i. */
/*             The strict upper triangular part of this array is not */
/*             referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,P). */

/*     K       (output) DOUBLE PRECISION array, dimension (LDK,P) */
/*             If JOBK = 'K', and INFO = 0, then the leading N-by-P part */
/*             of this array contains K , the Kalman filter gain matrix */
/*                                     i */
/*             at instant i. */
/*             If JOBK = 'N', or JOBK = 'K' and INFO = 1, then the */
/*             leading N-by-P part of this array contains AK , a matrix */
/*                                                          i */
/*             related to the Kalman filter gain matrix at instant i (see */
/*                                                            -1/2 */
/*             METHOD). Specifically, AK  = A P     C'(RINOV')    . */
/*                                      i    i i|i-1 i      i */

/*     LDK     INTEGER */
/*             The leading dimension of array K.   LDK >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If JOBK = 'K', then TOL is used to test for near */
/*                                               1/2 */
/*             singularity of the matrix (RINOV )   . If the user sets */
/*                                             i */
/*             TOL > 0, then the given value of TOL is used as a */
/*             lower bound for the reciprocal condition number of that */
/*             matrix; a matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. If the user */
/*             sets TOL <= 0, then an implicitly computed, default */
/*             tolerance, defined by TOLDEF = P*P*EPS, is used instead, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             routine DLAMCH). */
/*             Otherwise, TOL is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), */
/*             where LIWORK = P if JOBK = 'K', */
/*             and   LIWORK = 1 otherwise. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK.  If INFO = 0 and JOBK = 'K', DWORK(2) returns */
/*             an estimate of the reciprocal of the condition number */
/*                                        1/2 */
/*             (in the 1-norm) of (RINOV )   . */
/*                                      i */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(P+N)+2*P,N*(N+M+2)),     if JOBK = 'N'; */
/*             LDWORK >= MAX(2,N*(P+N)+2*P,N*(N+M+2),3*P), if JOBK = 'K'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*                                                        1/2 */
/*             = 1:  if JOBK = 'K' and the matrix (RINOV )   is singular, */
/*                                                      i           1/2 */
/*                   i.e., the condition number estimate of (RINOV ) */
/*                                                                i */
/*                   (in the 1-norm) exceeds 1/TOL.  The matrices S, AK , */
/*                               1/2                                   i */
/*                   and (RINOV )    have been computed. */
/*                             i */

/*     METHOD */

/*     The routine performs one recursion of the square root covariance */
/*     filter algorithm, summarized as follows: */

/*      |  1/2                      |     |         1/2          | */
/*      | R      C x S      0       |     | (RINOV )     0     0 | */
/*      |  i      i   i-1           |     |       i              | */
/*      |                      1/2  | T = |                      | */
/*      | 0      A x S    B x Q     |     |     AK       S     0 | */
/*      |         i   i-1  i   i    |     |       i       i      | */

/*          (Pre-array)                      (Post-array) */

/*     where T is an orthogonal transformation triangularizing the */
/*     pre-array. */

/*     The state covariance matrix P    is factorized as */
/*                                  i|i-1 */
/*        P     = S  S' */
/*         i|i-1   i  i */

/*     and one combined time and measurement update for the state X */
/*                                                                 i|i-1 */
/*     is given by */

/*        X     = A X      + K (Y - C X     ), */
/*         i+1|i   i i|i-1    i  i   i i|i-1 */

/*                          -1/2 */
/*     where K = AK (RINOV )     is the Kalman filter gain matrix and Y */
/*            i    i      i                                            i */
/*     is the observed output of the system. */

/*     The triangularization is done entirely via Householder */
/*     transformations exploiting the zero pattern of the pre-array. */

/*     REFERENCES */

/*     [1] Anderson, B.D.O. and Moore, J.B. */
/*         Optimal Filtering. */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1979. */

/*     [2] Verhaegen, M.H.G. and Van Dooren, P. */
/*         Numerical Aspects of Different Kalman Filter Implementations. */
/*         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986. */

/*     [3] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G. */
/*         Algorithm 675: FORTRAN Subroutines for Computing the Square */
/*         Root Covariance Filter and Square Root Information Filter in */
/*         Dense or Hessenberg Forms. */
/*         ACM Trans. Math. Software, 15, pp. 243-256, 1989. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires */

/*           3    2                               2   2 */
/*     (7/6)N  + N  x (5/2 x P + M) + N x (1/2 x M + P ) */

/*     operations and is backward stable (see [2]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine FB01ED by M. Vanbegin, */
/*     P. Van Dooren, and M.H.G. Verhaegen. */

/*     REVISIONS */

/*     February 20, 1998, November 20, 2003. */

/*     KEYWORDS */

/*     Kalman filtering, optimal filtering, orthogonal transformation, */
/*     recursive estimation, square-root covariance filtering, */
/*     square-root filtering. */

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
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    k_dim1 = *ldk;
    k_offset = 1 + k_dim1;
    k -= k_offset;
    --iwork;
    --dwork;

    /* Function Body */
    pn = *p + *n;
    n1 = max(1,*n);
    *info = 0;
    ljobk = lsame_(jobk, "K", (ftnlen)1, (ftnlen)1);
    lmultb = lsame_(multbq, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljobk && ! lsame_(jobk, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lmultb && ! lsame_(multbq, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*lds < n1) {
	*info = -7;
    } else if (*lda < n1) {
	*info = -9;
    } else if (*ldb < n1) {
	*info = -11;
    } else if (*ldq < 1 || ! lmultb && *ldq < *m) {
	*info = -13;
    } else if (*ldc < max(1,*p)) {
	*info = -15;
    } else if (*ldr < max(1,*p)) {
	*info = -17;
    } else if (*ldk < n1) {
	*info = -19;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 2, i__2 = pn * *n + (*p << 1), i__1 = max(i__1,i__2), i__2 = *
		n * (*n + *m + 2), i__1 = max(i__1,i__2), i__2 = *p * 3;
/* Computing MAX */
	i__3 = 1, i__4 = pn * *n + (*p << 1), i__3 = max(i__3,i__4), i__4 = *
		n * (*n + *m + 2);
	if (ljobk && *ldwork < max(i__1,i__2) || ! ljobk && *ldwork < max(
		i__3,i__4)) {
	    *info = -23;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("FB01QD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (ljobk) {
	    dwork[1] = 2.;
	    dwork[2] = 1.;
	} else {
	    dwork[1] = 1.;
	}
	return 0;
    }

/*     Construction of the needed part of the pre-array in DWORK. */
/*     To save workspace, only the blocks (1,2), (2,2), and (2,3) will be */
/*     constructed as shown below. */

/*     Storing A x S and C x S in the (1,1) and (2,1) blocks of DWORK, */
/*     respectively. */
/*     Workspace: need (N+P)*N. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], &pn, (ftnlen)4);
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[*n + 1], &pn, (ftnlen)4)
	    ;
    dtrmm_("Right", "Lower", "No transpose", "Non-unit", &pn, n, &c_b13, &s[
	    s_offset], lds, &dwork[1], &pn, (ftnlen)5, (ftnlen)5, (ftnlen)12, 
	    (ftnlen)8);

/*     Triangularization (2 steps). */

/*     Step 1: annihilate the matrix C x S. */
/*     Workspace: need (N+P)*N + 2*P. */

    itau = pn * *n + 1;
    jwork = itau + *p;

    mb04ld_("Full", p, n, n, &r__[r_offset], ldr, &dwork[*n + 1], &pn, &dwork[
	    1], &pn, &k[k_offset], ldk, &dwork[itau], &dwork[jwork], (ftnlen)
	    4);
    wrkopt = pn * *n + (*p << 1);

/*     Now, the workspace for C x S is no longer needed. */
/*     Adjust the leading dimension of DWORK, to save space for the */
/*     following computations. */

    dlacpy_("Full", n, n, &dwork[1], &pn, &dwork[1], n, (ftnlen)4);
    i12 = *n * *n + 1;

/*     Storing B x Q in the (1,2) block of DWORK. */
/*     Workspace: need N*(N+M). */

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[i12], n, (ftnlen)4);
    if (! lmultb) {
	dtrmm_("Right", "Lower", "No transpose", "Non-unit", n, m, &c_b13, &q[
		q_offset], ldq, &dwork[i12], n, (ftnlen)5, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = *n * (*n + *m);
    wrkopt = max(i__1,i__2);

/*     Step 2: LQ triangularization of the matrix [ A x S  B x Q ], where */
/*     A x S was modified at Step 1. */
/*     Workspace: need N*(N+M+2);  prefer N*(N+M+1)+N*NB. */

    itau = *n * (*n + *m) + 1;
    jwork = itau + *n;

    i__1 = *n + *m;
    i__2 = *ldwork - jwork + 1;
    dgelqf_(n, &i__1, &dwork[1], n, &dwork[itau], &dwork[jwork], &i__2, info);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Output S and K (if needed) and set the optimal workspace */
/*     dimension (and the reciprocal of the condition number estimate). */

    dlacpy_("Lower", n, n, &dwork[1], n, &s[s_offset], lds, (ftnlen)5);

    if (ljobk) {

/*        Compute K. */
/*        Workspace: need 3*P. */

	mb02od_("Right", "Lower", "No transpose", "Non-unit", "1-norm", n, p, 
		&c_b13, &r__[r_offset], ldr, &k[k_offset], ldk, &rcond, tol, &
		iwork[1], &dwork[1], info, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
		ftnlen)8, (ftnlen)6);
	if (*info == 0) {
/* Computing MAX */
	    i__1 = wrkopt, i__2 = *p * 3;
	    wrkopt = max(i__1,i__2);
	    dwork[2] = rcond;
	}
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of FB01QD *** */
} /* fb01qd_ */

