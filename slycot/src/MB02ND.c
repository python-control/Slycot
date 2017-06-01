/* MB02ND.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 0.;
static doublereal c_b5 = 1.;
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b85 = -1.;

/* Subroutine */ int mb02nd_(integer *m, integer *n, integer *l, integer *
	rank, doublereal *theta, doublereal *c__, integer *ldc, doublereal *x,
	 integer *ldx, doublereal *q, logical *inul, doublereal *tol, 
	doublereal *reltol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *iwarn, integer *info)
{
    /* System generated locals */
    integer c_dim1, c_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, p, i1, j1, n1, jf, kf, mc, ij;
    static doublereal hh;
    static integer kj;
    static doublereal cs;
    static integer mj, nj, nl, jv;
    static doublereal sn;
    static integer lw, ldf, mnl;
    static doublereal eps;
    static integer ioff;
    static doublereal temp;
    static integer ifail;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), mb04yd_(char *, char *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iwarm;
    static doublereal fnorm;
    static integer itaup, itauq;
    static doublereal first;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    static doublereal dummy[1];
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dgerqf_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), dlaset_(char *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dtrcon_(char *, char *, char *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal inprod;
    static integer ihoush;
    static logical lfirst;
    extern /* Subroutine */ int dormrq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical sufwrk;
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

/*     To solve the Total Least Squares (TLS) problem using a Partial */
/*     Singular Value Decomposition (PSVD) approach. */
/*     The TLS problem assumes an overdetermined set of linear equations */
/*     AX = B, where both the data matrix A as well as the observation */
/*     matrix B are inaccurate. The routine also solves determined and */
/*     underdetermined sets of equations by computing the minimum norm */
/*     solution. */
/*     It is assumed that all preprocessing measures (scaling, coordinate */
/*     transformations, whitening, ... ) of the data have been performed */
/*     in advance. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in the data matrix A and the */
/*             observation matrix B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in the data matrix A.  N >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the observation matrix B. */
/*             L >= 0. */

/*     RANK    (input/output) INTEGER */
/*             On entry, if RANK < 0, then the rank of the TLS */
/*             approximation [A+DA|B+DB] (r say) is computed by the */
/*             routine. */
/*             Otherwise, RANK must specify the value of r. */
/*             RANK <= min(M,N). */
/*             On exit, if RANK < 0 on entry and INFO = 0, then RANK */
/*             contains the computed rank of the TLS approximation */
/*             [A+DA|B+DB]. */
/*             Otherwise, the user-supplied value of RANK may be */
/*             changed by the routine on exit if the RANK-th and the */
/*             (RANK+1)-th singular values of C = [A|B] are considered */
/*             to be equal, or if the upper triangular matrix F (as */
/*             defined in METHOD) is (numerically) singular. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, if RANK < 0, then the rank of the TLS */
/*             approximation [A+DA|B+DB] is computed using THETA as */
/*             (min(M,N+L) - d), where d is the number of singular */
/*             values of [A|B] <= THETA. THETA >= 0.0. */
/*             Otherwise, THETA is an initial estimate (t say) for */
/*             computing a lower bound on the RANK largest singular */
/*             values of [A|B]. If THETA < 0.0 on entry however, then */
/*             t is computed by the routine. */
/*             On exit, if RANK >= 0 on entry, then THETA contains the */
/*             computed bound such that precisely RANK singular values */
/*             of C = [A|B] are greater than THETA + TOL. */
/*             Otherwise, THETA is unchanged. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N+L) */
/*             On entry, the leading M-by-(N+L) part of this array must */
/*             contain the matrices A and B. Specifically, the first N */
/*             columns must contain the data matrix A and the last L */
/*             columns the observation matrix B (right-hand sides). */
/*             On exit, if INFO = 0, the first N+L components of the */
/*             columns of this array whose index i corresponds with */
/*             INUL(i) = .TRUE., are the possibly transformed (N+L-RANK) */
/*             base vectors of the right singular subspace corresponding */
/*             to the singular values of C = [A|B] which are less than or */
/*             equal to THETA. Specifically, if L = 0, or if RANK = 0 and */
/*             IWARN <> 2, these vectors are indeed the base vectors */
/*             above. Otherwise, these vectors form the matrix V2, */
/*             transformed as described in Step 4 of the PTLS algorithm */
/*             (see METHOD). The TLS solution is computed from these */
/*             vectors. The other columns of array C contain no useful */
/*             information. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= max(1,M,N+L). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,L) */
/*             If INFO = 0, the leading N-by-L part of this array */
/*             contains the solution X to the TLS problem specified by */
/*             A and B. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= max(1,N). */

/*     Q       (output) DOUBLE PRECISION array, dimension */
/*             (max(1,2*min(M,N+L)-1)) */
/*             This array contains the partially diagonalized bidiagonal */
/*             matrix J computed from C, at the moment that the desired */
/*             singular subspace has been found. Specifically, the */
/*             leading p = min(M,N+L) entries of Q contain the diagonal */
/*             elements q(1),q(2),...,q(p) and the entries Q(p+1),Q(p+2), */
/*             ...,Q(2*p-1) contain the superdiagonal elements e(1),e(2), */
/*             ...,e(p-1) of J. */

/*     INUL    (output) LOGICAL array, dimension (N+L) */
/*             The indices of the elements of this array with value */
/*             .TRUE. indicate the columns in C containing the base */
/*             vectors of the right singular subspace of C from which */
/*             the TLS solution has been computed. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             This parameter defines the multiplicity of singular values */
/*             by considering all singular values within an interval of */
/*             length TOL as coinciding. TOL is used in checking how many */
/*             singular values are less than or equal to THETA. Also in */
/*             computing an appropriate upper bound THETA by a bisection */
/*             method, TOL is used as a stopping criterion defining the */
/*             minimum (absolute) subinterval width. TOL is also taken */
/*             as an absolute tolerance for negligible elements in the */
/*             QR/QL iterations. If the user sets TOL to be less than or */
/*             equal to 0, then the tolerance is taken as specified in */
/*             SLICOT Library routine MB04YD document. */

/*     RELTOL  DOUBLE PRECISION */
/*             This parameter specifies the minimum relative width of an */
/*             interval. When an interval is narrower than TOL, or than */
/*             RELTOL times the larger (in magnitude) endpoint, then it */
/*             is considered to be sufficiently small and bisection has */
/*             converged. If the user sets RELTOL to be less than */
/*             BASE * EPS, where BASE is machine radix and EPS is machine */
/*             precision (see LAPACK Library routine DLAMCH), then the */
/*             tolerance is taken as BASE * EPS. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+2*L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and DWORK(2) returns the reciprocal of the */
/*             condition number of the matrix F. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max(2, max(M,N+L) + 2*min(M,N+L), */
/*                          min(M,N+L) + LW + max(6*(N+L)-5, */
/*                                                L*L+max(N+L,3*L)), */
/*             where */
/*             LW = (N+L)*(N+L-1)/2,  if M >= N+L, */
/*             LW = M*(N+L-(M-1)/2),  if M <  N+L. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (N+L) */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warnings; */
/*             = 1:  if the rank of matrix C has been lowered because a */
/*                   singular value of multiplicity greater than 1 was */
/*                   found; */
/*             = 2:  if the rank of matrix C has been lowered because the */
/*                   upper triangular matrix F is (numerically) singular. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the maximum number of QR/QL iteration steps */
/*                   (30*MIN(M,N)) has been exceeded; */
/*             = 2:  if the computed rank of the TLS approximation */
/*                   [A+DA|B+DB] exceeds MIN(M,N). Try increasing the */
/*                   value of THETA or set the value of RANK to min(M,N). */

/*     METHOD */

/*     The method used is the Partial Total Least Squares (PTLS) approach */
/*     proposed by Van Huffel and Vandewalle [5]. */

/*     Let C = [A|B] denote the matrix formed by adjoining the columns of */
/*     B to the columns of A on the right. */

/*     Total Least Squares (TLS) definition: */
/*     ------------------------------------- */

/*       Given matrices A and B, find a matrix X satisfying */

/*            (A + DA) X = B + DB, */

/*       where A and DA are M-by-N matrices, B and DB are M-by-L matrices */
/*       and X is an N-by-L matrix. */
/*       The solution X must be such that the Frobenius norm of [DA|DB] */
/*       is a minimum and each column of B + DB is in the range of */
/*       A + DA. Whenever the solution is not unique, the routine singles */
/*       out the minimum norm solution X. */

/*     Let V denote the right singular subspace of C. Since the TLS */
/*     solution can be computed from any orthogonal basis of the subspace */
/*     of V corresponding to the smallest singular values of C, the */
/*     Partial Singular Value Decomposition (PSVD) can be used instead of */
/*     the classical SVD. The dimension of this subspace of V may be */
/*     determined by the rank of C or by an upper bound for those */
/*     smallest singular values. */

/*     The PTLS algorithm proceeds as follows (see [2 - 5]): */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*      (a) If M is large enough than N + L, transform C into upper */
/*          triangular form R by Householder transformations. */
/*      (b) Transform C (or R) into upper bidiagonal form */
/*          (p = min(M,N+L)): */

/*                     |q(1) e(1)  0   ...  0   | */
/*                (0)  | 0   q(2) e(2)      .   | */
/*               J   = | .                  .   | */
/*                     | .                e(p-1)| */
/*                     | 0             ... q(p) | */

/*          if M >= N + L, or lower bidiagonal form: */

/*                     |q(1)  0    0   ...  0     0   | */
/*                (0)  |e(1) q(2)  0        .     .   | */
/*               J   = | .                  .     .   | */
/*                     | .                 q(p)   .   | */
/*                     | 0             ... e(p-1) q(p)| */

/*          if M < N + L, using Householder transformations. */
/*          In the second case, transform the matrix to the upper */
/*          bidiagonal form by applying Givens rotations. */
/*      (c) Initialize the right singular base matrix with the identity */
/*          matrix. */

/*     Step 2: Partial diagonalization phase */
/*             ----------------------------- */
/*     If the upper bound THETA is not given, then compute THETA such */
/*     that precisely p - RANK singular values (p=min(M,N+L)) of the */
/*     bidiagonal matrix are less than or equal to THETA, using a */
/*     bisection method [5]. Diagonalize the given bidiagonal matrix J */
/*     partially, using either QL iterations (if the upper left diagonal */
/*     element of the considered bidiagonal submatrix is smaller than the */
/*     lower right diagonal element) or QR iterations, such that J is */
/*     split into unreduced bidiagonal submatrices whose singular values */
/*     are either all larger than THETA or are all less than or equal */
/*     to THETA. Accumulate the Givens rotations in V. */

/*     Step 3: Back transformation phase */
/*             ------------------------- */
/*     Apply the Householder transformations of Step 1(b) onto the base */
/*     vectors of V associated with the bidiagonal submatrices with all */
/*     singular values less than or equal to THETA. */

/*     Step 4: Computation of F and Y */
/*             ---------------------- */
/*     Let V2 be the matrix of the columns of V corresponding to the */
/*     (N + L - RANK) smallest singular values of C. */
/*     Compute with Householder transformations the matrices F and Y */
/*     such that: */

/*                       |VH   Y| */
/*              V2 x Q = |      | */
/*                       |0    F| */

/*     where Q is an orthogonal matrix, VH is an N-by-(N-RANK) matrix, */
/*     Y is an N-by-L matrix and F is an L-by-L upper triangular matrix. */
/*     If F is singular, then reduce the value of RANK by one and repeat */
/*     Steps 2, 3 and 4. */

/*     Step 5: Computation of the TLS solution */
/*             ------------------------------- */
/*     If F is non-singular then the solution X is obtained by solving */
/*     the following equations by forward elimination: */

/*              X F = -Y. */

/*     Notes: */
/*     If RANK is lowered in Step 4, some additional base vectors must */
/*     be computed in Step 2. The additional computations are kept to */
/*     a minimum. */
/*     If RANK is lowered in Step 4 but the multiplicity of the RANK-th */
/*     singular value is larger than 1, then the value of RANK is further */
/*     lowered with its multiplicity defined by the parameter TOL. This */
/*     is done at the beginning of Step 2 by calling SLICOT Library */
/*     routine MB03MD (from MB04YD), which estimates THETA using a */
/*     bisection method. If F in Step 4 is singular, then the computed */
/*     solution is infinite and hence does not satisfy the second TLS */
/*     criterion (see TLS definition). For these cases, Golub and */
/*     Van Loan [1] claim that the TLS problem has no solution. The */
/*     properties of these so-called nongeneric problems are described */
/*     in [6] and the TLS computations are generalized in order to solve */
/*     them. As proven in [6], the proposed generalization satisfies the */
/*     TLS criteria for any number L of observation vectors in B provided */
/*     that, in addition, the solution | X| is constrained to be */
/*                                     |-I| */
/*     orthogonal to all vectors of the form |w| which belong to the */
/*                                           |0| */
/*     space generated by the columns of the submatrix |Y|. */
/*                                                     |F| */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         An Analysis of the Total Least-Squares Problem. */
/*         SIAM J. Numer. Anal., 17, pp. 883-893, 1980. */

/*     [2] Van Huffel, S., Vandewalle, J. and Haegemans, A. */
/*         An Efficient and Reliable Algorithm for Computing the */
/*         Singular Subspace of a Matrix Associated with its Smallest */
/*         Singular Values. */
/*         J. Comput. and Appl. Math., 19, pp. 313-330, 1987. */

/*     [3] Van Huffel, S. */
/*         Analysis of the Total Least Squares Problem and its Use in */
/*         Parameter Estimation. */
/*         Doctoral dissertation, Dept. of Electr. Eng., Katholieke */
/*         Universiteit Leuven, Belgium, June 1987. */

/*     [4] Chan, T.F. */
/*         An Improved Algorithm for Computing the Singular Value */
/*         Decomposition. */
/*         ACM TOMS, 8, pp. 72-83, 1982. */

/*     [5] Van Huffel, S. and Vandewalle, J. */
/*         The Partial Total Least Squares Algorithm. */
/*         J. Comput. Appl. Math., 21, pp. 333-341, 1988. */

/*     [6] Van Huffel, S. and Vandewalle, J. */
/*         Analysis and Solution of the Nongeneric Total Least Squares */
/*         Problem. */
/*         SIAM J. Matr. Anal. and Appl., 9, pp. 360-372, 1988. */

/*     NUMERICAL ASPECTS */

/*     The computational efficiency of the PTLS algorithm compared with */
/*     the classical TLS algorithm (see [2 - 5]) is obtained by making */
/*     use of PSVD (see [1]) instead of performing the entire SVD. */
/*     Depending on the gap between the RANK-th and the (RANK+1)-th */
/*     singular values of C, the number (N + L - RANK) of base vectors to */
/*     be computed with respect to the column dimension (N + L) of C and */
/*     the desired accuracy RELTOL, the algorithm used by this routine is */
/*     approximately twice as fast as the classical TLS algorithm at the */
/*     expense of extra storage requirements, namely: */
/*       (N + L) x (N + L - 1)/2  if M >= N + L or */
/*       M x (N + L - (M - 1)/2)  if M <  N + L. */
/*     This is because the Householder transformations performed on the */
/*     rows of C in the bidiagonalization phase (see Step 1) must be kept */
/*     until the end (Step 5). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB02BD by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     June 30, 1997, Oct. 19, 2003, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Least-squares approximation, singular subspace, singular value */
/*     decomposition, singular values, total least-squares. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --q;
    --inul;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    *iwarn = 0;
    *info = 0;
    nl = *n + *l;
    k = max(*m,nl);
    p = min(*m,nl);
    if (*m >= nl) {
	lw = nl * (nl - 1) / 2;
    } else {
	lw = *m * nl - *m * (*m - 1) / 2;
    }
/* Computing MAX */
/* Computing MAX */
    i__3 = nl, i__4 = *l * 3;
    i__1 = nl * 6 - 5, i__2 = *l * *l + max(i__3,i__4);
    jv = p + lw + max(i__1,i__2);

/*     Test the input scalar arguments. */

    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*rank > min(*m,*n)) {
	*info = -4;
    } else if (*rank < 0 && *theta < 0.) {
	*info = -5;
    } else if (*ldc < max(1,k)) {
	*info = -7;
    } else if (*ldx < max(1,*n)) {
	*info = -9;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 2, i__2 = k + (p << 1), i__1 = max(i__1,i__2);
	if (*ldwork < max(i__1,jv)) {
	    *info = -16;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB02ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*m,nl) == 0) {
	if (*m == 0) {
	    dlaset_("Full", &nl, &nl, &c_b4, &c_b5, &c__[c_offset], ldc, (
		    ftnlen)4);
	    dlaset_("Full", n, l, &c_b4, &c_b4, &x[x_offset], ldx, (ftnlen)4);

	    i__1 = nl;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		inul[i__] = TRUE_;
/* L10: */
	    }

	}
	if (*rank >= 0) {
	    *theta = 0.;
	}
	*rank = 0;
	dwork[1] = 2.;
	dwork[2] = 1.;
	return 0;
    }

    wrkopt = 2;
    n1 = *n + 1;

    eps = dlamch_("Precision", (ftnlen)9);
    lfirst = TRUE_;

/*     Initializations. */

    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	inul[i__] = FALSE_;
	bwork[i__] = FALSE_;
/* L20: */
    }

    i__1 = nl;
    for (i__ = p + 1; i__ <= i__1; ++i__) {
	inul[i__] = TRUE_;
	bwork[i__] = FALSE_;
/* L40: */
    }

/*     Subroutine MB02ND solves a set of linear equations by a Total */
/*     Least Squares Approximation, based on the Partial SVD. */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*     1.a): If M is large enough than N+L, transform C into upper */
/*           triangular form R by Householder transformations. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/* Computing MAX */
    i__1 = nl, i__2 = ilaenv_(&c__6, "DGESVD", "NN", m, &nl, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)2);
    if (*m >= max(i__1,i__2)) {

/*        Workspace: need   2*(N+L), */
/*                   prefer N+L + (N+L)*NB. */

	itauq = 1;
	jwork = itauq + nl;
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(m, &nl, &c__[c_offset], ldc, &dwork[itauq], &dwork[jwork], &
		i__1, &ifail);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	if (nl > 1) {
	    i__1 = nl - 1;
	    i__2 = nl - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b4, &c_b4, &c__[c_dim1 + 2], 
		    ldc, (ftnlen)5);
	}
	mnl = nl;
    } else {
	mnl = *m;
    }

/*     1.b): Transform C (or R) into bidiagonal form Q using Householder */
/*           transformations. */
/*     Workspace: need   2*min(M,N+L) + max(M,N+L), */
/*                prefer 2*min(M,N+L) + (M+N+L)*NB. */

    itaup = 1;
    itauq = itaup + p;
    jwork = itauq + p;
    i__1 = *ldwork - jwork + 1;
    dgebrd_(&mnl, &nl, &c__[c_offset], ldc, &q[1], &q[p + 1], &dwork[itauq], &
	    dwork[itaup], &dwork[jwork], &i__1, &ifail);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     If the matrix is lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left. */

    if (*m < nl) {
	ioff = 0;

	i__1 = p - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dlartg_(&q[i__], &q[p + i__], &cs, &sn, &temp);
	    q[i__] = temp;
	    q[p + i__] = sn * q[i__ + 1];
	    q[i__ + 1] = cs * q[i__ + 1];
/* L60: */
	}

    } else {
	ioff = 1;
    }

/*     Store the Householder transformations performed onto the rows of C */
/*     in the extra storage locations DWORK(IHOUSH). */
/*     Workspace: need   LDW = min(M,N+L) + (N+L)*(N+L-1)/2, if M >= N+L, */
/*                       LDW = min(M,N+L) + M*(N+L-(M-1)/2), if M <  N+L; */
/*                prefer LDW = min(M,N+L) + (N+L)**2,        if M >= N+L, */
/*                       LDW = min(M,N+L) + M*(N+L),         if M <  N+L. */

    ihoush = itauq;
    mc = nl - ioff;
    kf = ihoush + p * nl;
/* Computing MAX */
/* Computing 2nd power */
    i__3 = nl;
/* Computing MAX */
    i__4 = nl, i__5 = *l * 3;
    i__1 = (*n + *l) * 6 - 5, i__2 = i__3 * i__3 + max(i__4,i__5) - 1;
    sufwrk = *ldwork >= kf + max(i__1,i__2);
    if (sufwrk) {

/*        Enough workspace for a fast algorithm. */

	dlacpy_("Upper", &p, &nl, &c__[c_offset], ldc, &dwork[ihoush], &p, (
		ftnlen)5);
	kj = kf;
/* Computing MAX */
	i__1 = wrkopt, i__2 = kf - 1;
	wrkopt = max(i__1,i__2);
    } else {

/*        Not enough workspace for a fast algorithm. */

	kj = ihoush;

	i__1 = min(p,mc);
	for (nj = 1; nj <= i__1; ++nj) {
	    j = mc - nj + 1;
	    dcopy_(&j, &c__[nj + (nj + ioff) * c_dim1], ldc, &dwork[kj], &
		    c__1);
	    kj += j;
/* L80: */
	}

    }

/*     1.c): Initialize the right singular base matrix V with the */
/*           identity matrix (V overwrites C). */

    dlaset_("Full", &nl, &nl, &c_b4, &c_b5, &c__[c_offset], ldc, (ftnlen)4);
    jv = kj;
    iwarm = 0;

/*     REPEAT */

/*     Compute the Householder matrix Q and matrices F and Y such that */
/*     F is nonsingular. */

/*     Step 2: Partial diagonalization phase. */
/*             ----------------------------- */
/*     Diagonalize the bidiagonal Q partially until convergence to */
/*     the desired right singular subspace. */
/*     Workspace: LDW + 6*(N+L)-5. */

L100:
    jwork = jv;
    i__1 = *ldwork - jwork + 1;
    mb04yd_("No U", "Update V", &p, &nl, rank, theta, &q[1], &q[p + 1], dummy,
	     &c__1, &c__[c_offset], ldc, &inul[1], tol, reltol, &dwork[jwork],
	     &i__1, iwarn, info, (ftnlen)4, (ftnlen)8);
/* Computing MAX */
    i__1 = wrkopt, i__2 = jwork + nl * 6 - 6;
    wrkopt = max(i__1,i__2);

    *iwarn = max(*iwarn,iwarm);
    if (*info > 0) {
	return 0;
    }

/*     Set pointers to the selected base vectors in the right singular */
/*     matrix of C. */

    k = 0;

    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (inul[i__]) {
	    ++k;
	    iwork[k] = i__;
	}
/* L120: */
    }

    if (k < *l) {

/*        Rank of the TLS approximation is larger than min(M,N). */

	*info = 2;
	return 0;
    }

/*     Step 3: Back transformation phase. */
/*             ------------------------- */
/*     Apply in backward order the Householder transformations (stored */
/*     in DWORK(IHOUSH)) performed onto the rows of C during the */
/*     bidiagonalization phase, to the selected base vectors (specified */
/*     by INUL(I) = .TRUE.). Already transformed vectors are those for */
/*     which BWORK(I) = .TRUE.. */

    kf = k;
    if (sufwrk && lfirst) {

/*        Enough workspace for a fast algorithm and first pass. */

	ij = jv;

	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&nl, &c__[iwork[j] * c_dim1 + 1], &c__1, &dwork[ij], &c__1)
		    ;
	    ij += nl;
/* L140: */
	}

/*        Workspace: need   LDW + (N+L)*K + K, */
/*                   prefer LDW + (N+L)*K + K*NB. */

	ij = jv;
	jwork = ij + nl * k;
	i__1 = *ldwork - jwork + 1;
	dormbr_("P vectors", "Left", "No transpose", &nl, &k, &mnl, &dwork[
		ihoush], &p, &dwork[itaup], &dwork[ij], &nl, &dwork[jwork], &
		i__1, &ifail, (ftnlen)9, (ftnlen)4, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

	i__1 = nl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inul[i__] && ! bwork[i__]) {
		bwork[i__] = TRUE_;
	    }
/* L160: */
	}

    } else {

/*        Not enough workspace for a fast algorithm or subsequent passes. */

	i__1 = nl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inul[i__] && ! bwork[i__]) {
		kj = jv;

		for (nj = min(p,mc); nj >= 1; --nj) {
		    j = mc - nj + 1;
		    kj -= j;
		    first = dwork[kj];
		    dwork[kj] = 1.;
		    dlarf_("Left", &j, &c__1, &dwork[kj], &c__1, &dwork[itaup 
			    + nj - 1], &c__[nj + ioff + i__ * c_dim1], ldc, &
			    dwork[jwork], (ftnlen)4);
		    dwork[kj] = first;
/* L170: */
		}

		bwork[i__] = TRUE_;
	    }
/* L180: */
	}
    }

    if (*rank <= 0) {
	*rank = 0;
    }
    if (min(*rank,*l) == 0) {
	if (sufwrk && lfirst) {
	    dlacpy_("Full", &nl, &k, &dwork[jv], &nl, &c__[c_offset], ldc, (
		    ftnlen)4);
	}
	dwork[1] = (doublereal) wrkopt;
	dwork[2] = 1.;
	return 0;
    }

/*     Step 4: Compute matrices F and Y */
/*             ------------------------ */
/*             using Householder transformation Q. */

/*     Compute the orthogonal matrix Q (in factorized form) and the */
/*     matrices F and Y using RQ factorization. It is assumed that, */
/*     generically, the last L rows of V2 matrix have full rank. */
/*     The code could not be the most efficient when RANK has been */
/*     lowered, because the already created zero pattern of the last */
/*     L rows of V2 matrix is not exploited. */

    if (sufwrk && lfirst) {

/*        Enough workspace for a fast algorithm and first pass. */
/*        Workspace: need   LDW1 + 2*L, */
/*                   prefer LDW1 + L + L*NB, where */
/*                          LDW1 = LDW + (N+L)*K; */

	itauq = jwork;
	jwork = itauq + *l;
	i__1 = *ldwork - jwork + 1;
	dgerqf_(l, &k, &dwork[jv + *n], &nl, &dwork[itauq], &dwork[jwork], &
		i__1, info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Workspace: need   LDW1 + N+L, */
/*                   prefer LDW1 + L + N*NB. */

	i__1 = *ldwork - jwork + 1;
	dormrq_("Right", "Transpose", n, &k, l, &dwork[jv + *n], &nl, &dwork[
		itauq], &dwork[jv], &nl, &dwork[jwork], &i__1, info, (ftnlen)
		5, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

	jf = jv + nl * (k - *l) + *n;
	ldf = nl;
	jwork = jf + ldf * *l - *n;
	i__1 = k - *l;
	dlaset_("Full", l, &i__1, &c_b4, &c_b4, &dwork[jv + *n], &ldf, (
		ftnlen)4);
	if (*l > 1) {
	    i__1 = *l - 1;
	    i__2 = *l - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b4, &c_b4, &dwork[jf + 1], &ldf,
		     (ftnlen)5);
	}
	ij = jv;

	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&nl, &dwork[ij], &c__1, &c__[iwork[j] * c_dim1 + 1], &c__1)
		    ;
	    ij += nl;
/* L200: */
	}

    } else {

/*        Not enough workspace for a fast algorithm or subsequent passes. */
/*        Workspace: LDW2 + N+L, where LDW2 = LDW + L*L. */

	i__ = nl;
	jf = jv;
	ldf = *l;
	jwork = jf + ldf * *l;
/* Computing MAX */
	i__1 = wrkopt, i__2 = jwork + nl - 1;
	wrkopt = max(i__1,i__2);

/*        WHILE ( ( K >= 1 ) .AND. ( I > N ) ) DO */
L220:
	if (k >= 1 && i__ > *n) {

	    i__1 = k;
	    for (j = 1; j <= i__1; ++j) {
		dwork[jwork + j - 1] = c__[i__ + iwork[j] * c_dim1];
/* L240: */
	    }

/*           Compute Householder transformation. */

	    dlarfg_(&k, &dwork[jwork + k - 1], &dwork[jwork], &c__1, &temp);
	    c__[i__ + iwork[k] * c_dim1] = dwork[jwork + k - 1];
	    if (temp != 0.) {

/*              Apply Householder transformation onto the selected base */
/*              vectors. */

		i__1 = i__ - 1;
		for (i1 = 1; i1 <= i__1; ++i1) {
		    inprod = c__[i1 + iwork[k] * c_dim1];

		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			inprod += dwork[jwork + j - 1] * c__[i1 + iwork[j] * 
				c_dim1];
/* L260: */
		    }

		    hh = inprod * temp;
		    c__[i1 + iwork[k] * c_dim1] -= hh;

		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			j1 = iwork[j];
			c__[i1 + j1 * c_dim1] -= dwork[jwork + j - 1] * hh;
			c__[i__ + j1 * c_dim1] = 0.;
/* L280: */
		    }

/* L300: */
		}

	    }
	    i__1 = i__ - *n;
	    dcopy_(&i__1, &c__[n1 + iwork[k] * c_dim1], &c__1, &dwork[jf + (
		    i__ - *n - 1) * *l], &c__1);
	    --k;
	    --i__;
	    goto L220;
	}
/*        END WHILE 220 */
    }

/*     Estimate the reciprocal condition number of the matrix F. */
/*     If F singular, lower the rank of the TLS approximation. */
/*     Workspace: LDW1 + 3*L or */
/*                LDW2 + 3*L. */

    dtrcon_("1-norm", "Upper", "Non-unit", l, &dwork[jf], &ldf, &rcond, &
	    dwork[jwork], &iwork[kf + 1], info, (ftnlen)6, (ftnlen)5, (ftnlen)
	    8);
/* Computing MAX */
    i__1 = wrkopt, i__2 = jwork + *l * 3 - 1;
    wrkopt = max(i__1,i__2);

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	dcopy_(n, &c__[iwork[kf - *l + j] * c_dim1 + 1], &c__1, &x[j * x_dim1 
		+ 1], &c__1);
/* L320: */
    }

    fnorm = dlantr_("1-norm", "Upper", "Non-unit", l, l, &dwork[jf], &ldf, &
	    dwork[jwork], (ftnlen)6, (ftnlen)5, (ftnlen)8);
    if (rcond <= eps * fnorm) {
	--(*rank);
	goto L340;
    }
    if (fnorm <= eps * dlange_("1-norm", n, l, &x[x_offset], ldx, &dwork[
	    jwork], (ftnlen)6)) {
	*rank -= *l;
	goto L340;
    } else {
	goto L400;
    }

L340:
    iwarm = 2;
    *theta = -1.;
    if (sufwrk && lfirst) {

/*           Rearrange the stored Householder transformations for */
/*           subsequent passes, taking care to avoid overwriting. */

	if (p < nl) {
	    kj = ihoush + nl * (nl - 1);
	    mj = ihoush + p * (nl - 1);

	    i__1 = nl;
	    for (nj = 1; nj <= i__1; ++nj) {
		for (j = p - 1; j >= 0; --j) {
		    dwork[kj + j] = dwork[mj + j];
/* L350: */
		}
		kj -= nl;
		mj -= p;
/* L360: */
	    }

	}
	kj = ihoush;
	mj = ihoush + nl * ioff;

	i__1 = min(p,mc);
	for (nj = 1; nj <= i__1; ++nj) {
	    i__2 = mc - nj;
	    for (j = 0; j <= i__2; ++j) {
		dwork[kj] = dwork[mj + j * p];
		++kj;
/* L370: */
	    }
	    mj = mj + nl + 1;
/* L380: */
	}

	jv = kj;
	lfirst = FALSE_;
    }
    goto L100;
/*     UNTIL ( F nonsingular, i.e., RCOND.GT.EPS*FNORM or */
/*                                  FNORM.GT.EPS*norm(Y) ) */
L400:

/*     Step 5: Compute TLS solution. */
/*             -------------------- */
/*     Solve X F = -Y  by forward elimination  (F is upper triangular). */

    dtrsm_("Right", "Upper", "No transpose", "Non-unit", n, l, &c_b85, &dwork[
	    jf], &ldf, &x[x_offset], ldx, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
	    ftnlen)8);

/*     Set the optimal workspace and reciprocal condition number of F. */

    dwork[1] = (doublereal) wrkopt;
    dwork[2] = rcond;

    return 0;
/* *** Last line of MB02ND *** */
} /* mb02nd_ */

