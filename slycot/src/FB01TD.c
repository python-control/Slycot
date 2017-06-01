/* FB01TD.f -- translated by f2c (version 20100827).
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
static doublereal c_b23 = 1.;

/* Subroutine */ int fb01td_(char *jobx, char *multrc, integer *n, integer *m,
	 integer *p, doublereal *sinv, integer *ldsinv, doublereal *ainv, 
	integer *ldainv, doublereal *ainvb, integer *ldainb, doublereal *rinv,
	 integer *ldrinv, doublereal *c__, integer *ldc, doublereal *qinv, 
	integer *ldqinv, doublereal *x, doublereal *rinvy, doublereal *z__, 
	doublereal *e, doublereal *tol, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobx_len, ftnlen multrc_len)
{
    /* System generated locals */
    integer ainv_dim1, ainv_offset, ainvb_dim1, ainvb_offset, c_dim1, 
	    c_offset, qinv_dim1, qinv_offset, rinv_dim1, rinv_offset, 
	    sinv_dim1, sinv_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static integer i__, m1, n1, i12, i13, i23, i32, i33, ii, ij, nm, np, mp1, 
	    ldw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau;
    extern /* Subroutine */ int mb04id_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), mb04kd_(char *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, ftnlen), mb02od_(char *, 
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    static logical ljobx;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static logical lmultr;
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
/*     iteration of the time-invariant Kalman filter. This update is */
/*     given for the square root information filter, using the condensed */
/*     controller Hessenberg form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBX    CHARACTER*1 */
/*             Indicates whether X    is to be computed as follows: */
/*                                i+1 */
/*             = 'X':  X    is computed and stored in array X; */
/*                      i+1 */
/*             = 'N':  X    is not required. */
/*                      i+1 */

/*     MULTRC  CHARACTER*1             -1/2 */
/*             Indicates how matrices R     and C    are to be passed to */
/*                                     i+1       i+1 */
/*             the routine as follows: */
/*             = 'P':  Array RINV is not used and the array C must */
/*                                          -1/2 */
/*                     contain the product R    C   ; */
/*                                          i+1  i+1 */
/*             = 'N':  Arrays RINV and C must contain the matrices */
/*                     as described below. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*                       -1      -1 */
/*             matrices S   and A  .  N >= 0. */
/*                       i */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e., the order of the matrix */
/*              -1/2 */
/*             Q    .  M >= 0. */
/*              i */

/*     P       (input) INTEGER */
/*             The actual output dimension, i.e., the order of the matrix */
/*              -1/2 */
/*             R    .  P >= 0. */
/*              i+1 */

/*     SINV    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDSINV,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*                                 -1 */
/*             array must contain S  , the inverse of the square root */
/*                                 i */
/*             (right Cholesky factor) of the state covariance matrix */
/*             P    (hence the information square root) at instant i. */
/*              i|i */
/*             On exit, the leading N-by-N upper triangular part of this */
/*                             -1 */
/*             array contains S   , the inverse of the square root (right */
/*                             i+1 */
/*             Cholesky factor) of the state covariance matrix P */
/*                                                              i+1|i+1 */
/*             (hence the information square root) at instant i+1. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDSINV  INTEGER */
/*             The leading dimension of array SINV.  LDSINV >= MAX(1,N). */

/*     AINV    (input) DOUBLE PRECISION array, dimension (LDAINV,N) */
/*                                                                 -1 */
/*             The leading N-by-N part of this array must contain A  , */
/*             the inverse of the state transition matrix of the discrete */
/*             system in controller Hessenberg form (e.g., as produced by */
/*             SLICOT Library Routine TB01MD). */

/*     LDAINV  INTEGER */
/*             The leading dimension of array AINV.  LDAINV >= MAX(1,N). */

/*     AINVB   (input) DOUBLE PRECISION array, dimension (LDAINB,M) */
/*                                                                  -1 */
/*             The leading N-by-M part of this array must contain  A  B, */
/*                             -1 */
/*             the product of A   and the input weight matrix B of the */
/*             discrete system, in upper controller Hessenberg form */
/*             (e.g., as produced by SLICOT Library Routine TB01MD). */

/*     LDAINB  INTEGER */
/*             The leading dimension of array AINVB.  LDAINB >= MAX(1,N). */

/*     RINV    (input) DOUBLE PRECISION array, dimension (LDRINV,*) */
/*             If MULTRC = 'N', then the leading P-by-P upper triangular */
/*                                              -1/2 */
/*             part of this array must contain R    , the inverse of the */
/*                                              i+1 */
/*             covariance square root (right Cholesky factor) of the */
/*             output (measurement) noise (hence the information square */
/*             root) at instant i+1. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */
/*             Otherwise, RINV is not referenced and can be supplied as a */
/*             dummy array (i.e., set parameter LDRINV = 1 and declare */
/*             this array to be RINV(1,1) in the calling program). */

/*     LDRINV  INTEGER */
/*             The leading dimension of array RINV. */
/*             LDRINV >= MAX(1,P) if MULTRC = 'N'; */
/*             LDRINV >= 1        if MULTRC = 'P'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain C   , */
/*                                                       -1/2      i+1 */
/*             the output weight matrix (or the product R    C    if */
/*                                                       i+1  i+1 */
/*             MULTRC = 'P') of the discrete system at instant i+1. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     QINV    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQINV,M) */
/*             On entry, the leading M-by-M upper triangular part of this */
/*                                 -1/2 */
/*             array must contain Q    , the inverse of the covariance */
/*                                 i */
/*             square root (right Cholesky factor) of the input (process) */
/*             noise (hence the information square root) at instant i. */
/*             On exit, the leading M-by-M upper triangular part of this */
/*                                    -1/2 */
/*             array contains (QINOV )    , the inverse of the covariance */
/*                                  i */
/*             square root (right Cholesky factor) of the process noise */
/*             innovation (hence the information square root) at */
/*             instant i. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDQINV  INTEGER */
/*             The leading dimension of array QINV.  LDQINV >= MAX(1,M). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain X , the estimated */
/*                                                i */
/*             filtered state at instant i. */
/*             On exit, if JOBX = 'X', and INFO = 0, then this array */
/*             contains X   , the estimated filtered state at */
/*                       i+1 */
/*             instant i+1. */
/*             On exit, if JOBX = 'N', or JOBX = 'X' and INFO = 1, then */
/*                                  -1 */
/*             this array contains S   X   . */
/*                                  i+1 i+1 */

/*     RINVY   (input) DOUBLE PRECISION array, dimension (P) */
/*                                      -1/2 */
/*             This array must contain R    Y   , the product of the */
/*                                      i+1  i+1 */
/*                                      -1/2 */
/*             upper triangular matrix R     and the measured output */
/*                                      i+1 */
/*             vector Y    at instant i+1. */
/*                     i+1 */

/*     Z       (input) DOUBLE PRECISION array, dimension (M) */
/*             This array must contain Z , the mean value of the state */
/*                                      i */
/*             process noise at instant i. */

/*     E       (output) DOUBLE PRECISION array, dimension (P) */
/*             This array contains E   , the estimated error at instant */
/*                                  i+1 */
/*             i+1. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If JOBX = 'X', then TOL is used to test for near */
/*                                        -1 */
/*             singularity of the matrix S   . If the user sets */
/*                                        i+1 */
/*             TOL > 0, then the given value of TOL is used as a */
/*             lower bound for the reciprocal condition number of that */
/*             matrix; a matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. If the user */
/*             sets TOL <= 0, then an implicitly computed, default */
/*             tolerance, defined by TOLDEF = N*N*EPS, is used instead, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             routine DLAMCH). */
/*             Otherwise, TOL is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             where LIWORK = N if JOBX = 'X', */
/*             and   LIWORK = 1 otherwise. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK.  If INFO = 0 and JOBX = 'X', DWORK(2) returns */
/*             an estimate of the reciprocal of the condition number */
/*                                 -1 */
/*             (in the 1-norm) of S   . */
/*                                 i+1 */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1)), */
/*                                 if JOBX = 'N'; */
/*             LDWORK >= MAX(2,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1), */
/*                           3*N), if JOBX = 'X'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value;                        -1 */
/*             = 1:  if JOBX = 'X' and the matrix S    is singular, */
/*                                                 i+1       -1 */
/*                   i.e., the condition number estimate of S    (in the */
/*                                                           i+1 */
/*                                                         -1    -1/2 */
/*                   1-norm) exceeds 1/TOL.  The matrices S   , Q */
/*                                                         i+1   i */
/*                   and E have been computed. */

/*     METHOD */

/*     The routine performs one recursion of the square root information */
/*     filter algorithm, summarized as follows: */

/*       |    -1/2             -1/2    |     |         -1/2             | */
/*       |   Q         0      Q    Z   |     | (QINOV )     *     *     | */
/*       |    i                i    i  |     |       i                  | */
/*       |                             |     |                          | */
/*       |           -1/2      -1/2    |     |             -1    -1     | */
/*     T |    0     R    C    R    Y   |  =  |    0       S     S   X   | */
/*       |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1| */
/*       |                             |     |                          | */
/*       |  -1 -1     -1 -1    -1      |     |                          | */
/*       | S  A  B   S  A     S  X     |     |    0         0     E     | */
/*       |  i         i        i  i    |     |                     i+1  | */

/*                   (Pre-array)                      (Post-array) */

/*     where T is an orthogonal transformation triangularizing the */
/*                        -1/2 */
/*     pre-array, (QINOV )     is the inverse of the covariance square */
/*                      i */
/*     root (right Cholesky factor) of the process noise innovation */
/*                                                            -1  -1 */
/*     (hence the information square root) at instant i and (A  ,A  B) is */
/*     in upper controller Hessenberg form. */

/*     An example of the pre-array is given below (where N = 6, M = 2, */
/*     and P = 3): */

/*         |x x |             | x| */
/*         |  x |             | x| */
/*         _______________________ */
/*         |    | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         _______________________ */
/*         |x x | x x x x x x | x| */
/*         |  x | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         |    |   x x x x x | x| */
/*         |    |     x x x x | x| */
/*         |    |       x x x | x| */

/*     The inverse of the corresponding state covariance matrix P */
/*                                                               i+1|i+1 */
/*     (hence the information matrix I) is then factorized as */

/*                    -1         -1     -1 */
/*         I       = P       = (S   )' S */
/*          i+1|i+1   i+1|i+1    i+1    i+1 */

/*     and one combined time and measurement update for the state is */
/*     given by X   . */
/*               i+1 */

/*     The triangularization is done entirely via Householder */
/*     transformations exploiting the zero pattern of the pre-array. */

/*     REFERENCES */

/*     [1] Anderson, B.D.O. and Moore, J.B. */
/*         Optimal Filtering. */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1979. */

/*     [2] Van Dooren, P. and Verhaegen, M.H.G. */
/*         Condensed Forms for Efficient Time-Invariant Kalman Filtering. */
/*         SIAM J. Sci. Stat. Comp., 9. pp. 516-530, 1988. */

/*     [3] Verhaegen, M.H.G. and Van Dooren, P. */
/*         Numerical Aspects of Different Kalman Filter Implementations. */
/*         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986. */

/*     [4] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G. */
/*         Algorithm 675: FORTRAN Subroutines for Computing the Square */
/*         Root Covariance Filter and Square Root Information Filter in */
/*         Dense or Hessenberg Forms. */
/*         ACM Trans. Math. Software, 15, pp. 243-256, 1989. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*           3    2                           2          3 */
/*     (1/6)N  + N x (3/2 x M + P) + 2 x N x M  + 2/3 x M */

/*     operations and is backward stable (see [3]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine FB01HD by M. Vanbegin, */
/*     P. Van Dooren, and M.H.G. Verhaegen. */

/*     REVISIONS */

/*     February 20, 1998, November 20, 2003, February 14, 2004. */

/*     KEYWORDS */

/*     Controller Hessenberg form, Kalman filtering, optimal filtering, */
/*     orthogonal transformation, recursive estimation, square-root */
/*     filtering, square-root information filtering. */

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
    sinv_dim1 = *ldsinv;
    sinv_offset = 1 + sinv_dim1;
    sinv -= sinv_offset;
    ainv_dim1 = *ldainv;
    ainv_offset = 1 + ainv_dim1;
    ainv -= ainv_offset;
    ainvb_dim1 = *ldainb;
    ainvb_offset = 1 + ainvb_dim1;
    ainvb -= ainvb_offset;
    rinv_dim1 = *ldrinv;
    rinv_offset = 1 + rinv_dim1;
    rinv -= rinv_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    qinv_dim1 = *ldqinv;
    qinv_offset = 1 + qinv_dim1;
    qinv -= qinv_offset;
    --x;
    --rinvy;
    --z__;
    --e;
    --iwork;
    --dwork;

    /* Function Body */
    np = *n + *p;
    nm = *n + *m;
    n1 = max(1,*n);
    m1 = max(1,*m);
    mp1 = *m + 1;
    *info = 0;
    ljobx = lsame_(jobx, "X", (ftnlen)1, (ftnlen)1);
    lmultr = lsame_(multrc, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljobx && ! lsame_(jobx, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lmultr && ! lsame_(multrc, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*ldsinv < n1) {
	*info = -7;
    } else if (*ldainv < n1) {
	*info = -9;
    } else if (*ldainb < n1) {
	*info = -11;
    } else if (*ldrinv < 1 || ! lmultr && *ldrinv < *p) {
	*info = -13;
    } else if (*ldc < max(1,*p)) {
	*info = -15;
    } else if (*ldqinv < m1) {
	*info = -17;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n - 1;
	i__1 = 2, i__2 = *n * (nm + *m) + *m * 3, i__1 = max(i__1,i__2), i__2 
		= np * (*n + 1) + *n + max(i__3,mp1), i__1 = max(i__1,i__2), 
		i__2 = *n * 3;
/* Computing MAX */
/* Computing MAX */
	i__6 = *n - 1;
	i__4 = 1, i__5 = *n * (nm + *m) + *m * 3, i__4 = max(i__4,i__5), i__5 
		= np * (*n + 1) + *n + max(i__6,mp1);
	if (ljobx && *ldwork < max(i__1,i__2) || ! ljobx && *ldwork < max(
		i__4,i__5)) {
	    *info = -25;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("FB01TD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (max(*n,*p) == 0) {
	if (ljobx) {
	    dwork[1] = 2.;
	    dwork[2] = 1.;
	} else {
	    dwork[1] = 1.;
	}
	return 0;
    }

/*     Construction of the needed part of the pre-array in DWORK. */
/*     To save workspace, only the blocks (1,3), (3,1)-(3,3), (2,2), and */
/*     (2,3) will be constructed when needed as shown below. */

/*     Storing SINV x AINVB and SINV x AINV in the (1,1) and (1,2) */
/*     blocks of DWORK, respectively. The upper trapezoidal structure of */
/*     [ AINVB AINV ] is fully exploited. Specifically, if M <= N, the */
/*     following partition is used: */

/*       [ S1  S2 ] [ B1  A1 A3 ] */
/*       [ 0   S3 ] [ 0   A2 A4 ], */

/*     where B1, A3, and S1 are M-by-M matrices, A1 and S2 are */
/*     M-by-(N-M), A2 and S3 are (N-M)-by-(N-M), A4 is (N-M)-by-M, and */
/*     B1, S1, A2, and S3 are upper triangular. The right hand side */
/*     matrix above is stored in the workspace. If M > N, the partition */
/*     is [ SINV ] [ B1 B2  A ], where B1 is N-by-N, B2 is N-by-(M-N), */
/*     and B1 and SINV are upper triangular. */
/*     The variables called Ixy define the starting positions where the */
/*     (x,y) blocks of the pre-array are initially stored in DWORK. */
/*     Workspace: need N*(M+N). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    ldw = n1;
    i32 = *n * *m + 1;

    dlacpy_("Upper", n, m, &ainvb[ainvb_offset], ldainb, &dwork[1], &ldw, (
	    ftnlen)5);
    i__1 = min(*m,*n);
    dlacpy_("Full", &i__1, n, &ainv[ainv_offset], ldainv, &dwork[i32], &ldw, (
	    ftnlen)4);
    if (*n > *m) {
	i__1 = *n - *m;
	dlacpy_("Upper", &i__1, n, &ainv[mp1 + ainv_dim1], ldainv, &dwork[i32 
		+ *m], &ldw, (ftnlen)5);
    }

/*                    [ B1  A1 ] */
/*     Compute SINV x [ 0   A2 ] or SINV x B1 as a product of upper */
/*     triangular matrices. */
/*     Workspace: need N*(M+N+1). */

    ii = 1;
    i13 = *n * nm + 1;
/* Computing MAX */
    i__1 = 1, i__2 = *n * nm + *n;
    wrkopt = max(i__1,i__2);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(&i__, &dwork[ii], &c__1, &dwork[i13], &c__1);
	dtrmv_("Upper", "No transpose", "Non-unit", &i__, &sinv[sinv_offset], 
		ldsinv, &dwork[i13], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	dcopy_(&i__, &dwork[i13], &c__1, &dwork[ii], &c__1);
	ii += *n;
/* L10: */
    }

/*                    [ A3 ] */
/*     Compute SINV x [ A4 ] or SINV x [ B2 A ]. */

    dtrmm_("Left", "Upper", "No transpose", "Non-unit", n, m, &c_b23, &sinv[
	    sinv_offset], ldsinv, &dwork[ii], &ldw, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

/*     Storing the process noise mean value in (1,3) block of DWORK. */
/*     Workspace: need N*(M+N) + M. */

    dcopy_(m, &z__[1], &c__1, &dwork[i13], &c__1);
    dtrmv_("Upper", "No transpose", "Non-unit", m, &qinv[qinv_offset], ldqinv,
	     &dwork[i13], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*     Computing SINV x X in X. */

    dtrmv_("Upper", "No transpose", "Non-unit", n, &sinv[sinv_offset], ldsinv,
	     &x[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*     Triangularization (2 steps). */

/*     Step 1: annihilate the matrix SINV x AINVB. */
/*     Workspace: need N*(N+2*M) + 3*M. */

    i12 = i13 + *m;
    itau = i12 + *m * *n;
    jwork = itau + *m;

    mb04kd_("Upper", m, n, n, &qinv[qinv_offset], ldqinv, &dwork[1], &ldw, &
	    dwork[i32], &ldw, &dwork[i12], &m1, &dwork[itau], &dwork[jwork], (
	    ftnlen)5);
/* Computing MAX */
    i__1 = wrkopt, i__2 = *n * (nm + *m) + *m * 3;
    wrkopt = max(i__1,i__2);

    if (*n == 0) {
	dcopy_(p, &rinvy[1], &c__1, &e[1], &c__1);
	if (ljobx) {
	    dwork[2] = 1.;
	}
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

/*     Apply the transformations to the last column of the pre-array. */
/*     (Only the updated (3,3) block is now needed.) */

    ij = 1;

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = min(i__,*n);
	i__3 = min(i__,*n);
	d__1 = -dwork[itau + i__ - 1] * (dwork[i13 + i__ - 1] + ddot_(&i__3, &
		dwork[ij], &c__1, &x[1], &c__1));
	daxpy_(&i__2, &d__1, &dwork[ij], &c__1, &x[1], &c__1);
	ij += *n;
/* L20: */
    }

/*     Now, the workspace for SINV x AINVB, as well as for the updated */
/*     (1,2) block of the pre-array, are no longer needed. */
/*     Move the computed (3,2) and (3,3) blocks of the pre-array in the */
/*     (1,1) and (1,2) block positions of DWORK, to save space for the */
/*     following computations. */
/*     Then, adjust the implicitly defined leading dimension of DWORK, */
/*     to make space for storing the (2,2) and (2,3) blocks of the */
/*     pre-array. */
/*     Workspace: need (P+N)*(N+1). */

    i__1 = min(*m,*n);
    dlacpy_("Full", &i__1, n, &dwork[i32], &ldw, &dwork[1], &ldw, (ftnlen)4);
    if (*n > *m) {
	i__1 = *n - *m;
	dlacpy_("Upper", &i__1, n, &dwork[i32 + *m], &ldw, &dwork[mp1], &ldw, 
		(ftnlen)5);
    }
    ldw = max(1,np);

    for (i__ = *n; i__ >= 1; --i__) {
/* Computing MIN */
	i__1 = *n, i__2 = i__ + *m;
	for (ij = min(i__1,i__2); ij >= 1; --ij) {
	    dwork[np * (i__ - 1) + *p + ij] = dwork[*n * (i__ - 1) + ij];
/* L30: */
	}
/* L40: */
    }

/*     Copy of RINV x C in the (1,1) block of DWORK. */

    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], &ldw, (ftnlen)4);
    if (! lmultr) {
	dtrmm_("Left", "Upper", "No transpose", "Non-unit", p, n, &c_b23, &
		rinv[rinv_offset], ldrinv, &dwork[1], &ldw, (ftnlen)4, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
    }

/*     Copy the inclusion measurement in the (1,2) block and the updated */
/*     X in the (2,2) block of DWORK. */

    i23 = np * *n + 1;
    i33 = i23 + *p;
    dcopy_(p, &rinvy[1], &c__1, &dwork[i23], &c__1);
    dcopy_(n, &x[1], &c__1, &dwork[i33], &c__1);
/* Computing MAX */
    i__1 = wrkopt, i__2 = np * (*n + 1);
    wrkopt = max(i__1,i__2);

/*     Step 2: QR factorization of the first block column of the matrix */

/*        [ RINV x C     RINV x Y ], */
/*        [ SINV x AINV  SINV x X ] */

/*     where the second block row was modified at Step 1. */
/*     Workspace: need   (P+N)*(N+1) + N + MAX(N-1,M+1); */
/*                prefer (P+N)*(N+1) + N + (M+1)*NB, where NB is the */
/*                       optimal block size for DGEQRF called in MB04ID. */

    itau = i23 + np;
    jwork = itau + *n;

/* Computing MAX */
    i__2 = *n - mp1;
    i__1 = max(i__2,0);
    i__3 = *ldwork - jwork + 1;
    mb04id_(&np, n, &i__1, &c__1, &dwork[1], &ldw, &dwork[i23], &ldw, &dwork[
	    itau], &dwork[jwork], &i__3, info);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Output SINV, X, and E and set the optimal workspace dimension */
/*     (and the reciprocal of the condition number estimate). */

    dlacpy_("Upper", n, n, &dwork[1], &ldw, &sinv[sinv_offset], ldsinv, (
	    ftnlen)5);
    dcopy_(n, &dwork[i23], &c__1, &x[1], &c__1);
    if (*p > 0) {
	dcopy_(p, &dwork[i23 + *n], &c__1, &e[1], &c__1);
    }

    if (ljobx) {

/*        Compute X. */
/*        Workspace: need 3*N. */

	mb02od_("Left", "Upper", "No transpose", "Non-unit", "1-norm", n, &
		c__1, &c_b23, &sinv[sinv_offset], ldsinv, &x[1], n, &rcond, 
		tol, &iwork[1], &dwork[1], info, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8, (ftnlen)6);
	if (*info == 0) {
/* Computing MAX */
	    i__1 = wrkopt, i__2 = *n * 3;
	    wrkopt = max(i__1,i__2);
	    dwork[2] = rcond;
	}
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of FB01TD*** */
} /* fb01td_ */

