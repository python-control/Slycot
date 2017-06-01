/* AB13BD.f -- translated by f2c (version 20100827).
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

static logical c_false = FALSE_;
static doublereal c_b16 = 1.;

doublereal ab13bd_(char *dico, char *jobn, integer *n, integer *m, integer *p,
	 doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal 
	*c__, integer *ldc, doublereal *d__, integer *ldd, integer *nq, 
	doublereal *tol, doublereal *dwork, integer *ldwork, integer *iwarn, 
	integer *info, ftnlen dico_len, ftnlen jobn_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer ku, nr, kcr, kdr, krw, ktau, mxnp;
    extern /* Subroutine */ int sb08dd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen);
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal s2norm;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal wrkopt;


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

/*     To compute the H2 or L2 norm of the transfer-function matrix G */
/*     of the system (A,B,C,D). G must not have poles on the imaginary */
/*     axis, for a continuous-time system, or on the unit circle, for */
/*     a discrete-time system. If the H2-norm is computed, the system */
/*     must be stable. */

/*     FUNCTION VALUE */

/*     AB13BD   DOUBLE PRECISION */
/*              The H2-norm of G, if JOBN = 'H', or the L2-norm of G, */
/*              if JOBN = 'L' (if INFO = 0). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBN    CHARACTER*1 */
/*             Specifies the norm to be computed as follows: */
/*             = 'H':  the H2-norm; */
/*             = 'L':  the L2-norm. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of the */
/*             matrix B, and the number of columns of the matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices B and D. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrices C and D. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix of the system. */
/*             On exit, the leading NQ-by-NQ part of this array contains */
/*             the state dynamics matrix (in a real Schur form) of the */
/*             numerator factor Q of the right coprime factorization with */
/*             inner denominator of G (see METHOD). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix of the system. */
/*             On exit, the leading NQ-by-M part of this array contains */
/*             the input/state matrix of the numerator factor Q of the */
/*             right coprime factorization with inner denominator of G */
/*             (see METHOD). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix of the system. */
/*             On exit, the leading P-by-NQ part of this array contains */
/*             the state/output matrix of the numerator factor Q of the */
/*             right coprime factorization with inner denominator of G */
/*             (see METHOD). */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix of the system. */
/*             If DICO = 'C', D must be a null matrix. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the input/output matrix of the numerator factor Q of */
/*             the right coprime factorization with inner denominator */
/*             of G (see METHOD). */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NQ      (output) INTEGER */
/*             The order of the resulting numerator Q of the right */
/*             coprime factorization with inner denominator of G (see */
/*             METHOD). */
/*             Generally, NQ = N - NS, where NS is the number of */
/*             uncontrollable unstable eigenvalues. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of */
/*             B are considered zero (used for controllability tests). */
/*             If the user sets TOL <= 0, then an implicitly computed, */
/*             default tolerance, defined by  TOLDEF = N*EPS*NORM(B), */
/*             is used instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH) and NORM(B) denotes */
/*             the 1-norm of B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX( 1, M*(N+M) + MAX( N*(N+5), M*(M+2), 4*P ), */
/*                               N*( MAX( N, P ) + 4 ) + MIN( N, P ) ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = K:  K violations of the numerical stability condition */
/*                   occured during the assignment of eigenvalues in */
/*                   computing the right coprime factorization with inner */
/*                   denominator of G (see the SLICOT subroutine SB08DD). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to a real Schur form failed; */
/*             = 2:  a failure was detected during the reordering of the */
/*                   real Schur form of A, or in the iterative process */
/*                   for reordering the eigenvalues of Z'*(A + B*F)*Z */
/*                   along the diagonal (see SLICOT routine SB08DD); */
/*             = 3:  if DICO = 'C' and the matrix A has a controllable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' */
/*                   and A has a controllable eigenvalue on the unit */
/*                   circle; */
/*             = 4:  the solution of Lyapunov equation failed because */
/*                   the equation is singular; */
/*             = 5:  if DICO = 'C' and D is a nonzero matrix; */
/*             = 6:  if JOBN = 'H' and the system is unstable. */

/*     METHOD */

/*     The subroutine is based on the algorithms proposed in [1] and [2]. */

/*     If the given transfer-function matrix G is unstable, then a right */
/*     coprime factorization with inner denominator of G is first */
/*     computed */
/*               -1 */
/*        G = Q*R  , */

/*     where Q and R are stable transfer-function matrices and R is */
/*     inner. If G is stable, then Q = G and R = I. */
/*     Let (AQ,BQ,CQ,DQ) be the state-space representation of Q. */

/*     If DICO = 'C', then the L2-norm of G is computed as */

/*        NORM2(G) = NORM2(Q) = SQRT(TRACE(BQ'*X*BQ)), */

/*     where X satisfies the continuous-time Lyapunov equation */

/*        AQ'*X + X*AQ + CQ'*CQ = 0. */

/*     If DICO = 'D', then the l2-norm of G is computed as */

/*        NORM2(G) = NORM2(Q) = SQRT(TRACE(BQ'*X*BQ+DQ'*DQ)), */

/*     where X satisfies the discrete-time Lyapunov equation */

/*        AQ'*X*AQ - X + CQ'*CQ = 0. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         On computing 2-norms of transfer-function matrices. */
/*         Proc. 1992 ACC, Chicago, June 1992. */

/*     [2] Varga A. */
/*         A Schur method for computing coprime factorizations with */
/*         inner denominators and applications in model reduction. */
/*         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993. */

/*     NUMERICAL ASPECTS */
/*                                            3 */
/*     The algorithm requires no more than 14N  floating point */
/*     operations. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine SL2NRM. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Jan. 2003, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Coprime factorization, Lyapunov equation, multivariable system, */
/*     state-space model, system norms. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External functions .. */
/*     .. External subroutines .. */
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
    --dwork;

    /* Function Body */
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    *info = 0;
    *iwarn = 0;

/*     Check the scalar input parameters. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (lsame_(jobn, "H", (ftnlen)1, (ftnlen)1) || lsame_(jobn, 
	    "L", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*p)) {
	*info = -11;
    } else if (*ldd < max(1,*p)) {
	*info = -13;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n * (*n + 5), i__4 = *m * (*m + 2), i__3 = max(i__3,i__4), 
		i__4 = *p << 2;
	i__1 = 1, i__2 = *m * (*n + *m) + max(i__3,i__4), i__1 = max(i__1,
		i__2), i__2 = *n * (max(*n,*p) + 4) + min(*n,*p);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -17;
	}
    }
    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB13BD", &i__1, (ftnlen)6);
	return ret_val;
    }

/*     Compute the Frobenius norm of D. */

    s2norm = dlange_("Frobenius", p, m, &d__[d_offset], ldd, &dwork[1], (
	    ftnlen)9);
    if (! discr && s2norm != 0.) {
	*info = 5;
	return ret_val;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nq = 0;
	ret_val = 0.;
	dwork[1] = 1.;
	return ret_val;
    }

    kcr = 1;
    kdr = kcr + *m * *n;
    krw = kdr + *m * *m;

/*     Compute the right coprime factorization with inner denominator */
/*     of G. */

/*     Workspace needed:      M*(N+M); */
/*     Additional workspace:  need MAX( N*(N+5), M*(M+2), 4*M, 4*P ); */
/*                            prefer larger. */

    i__1 = *ldwork - krw + 1;
    sb08dd_(dico, n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, nq, &nr, &dwork[kcr], m, &
	    dwork[kdr], m, tol, &dwork[krw], &i__1, iwarn, info, (ftnlen)1);
    if (*info != 0) {
	return ret_val;
    }

    wrkopt = dwork[krw] + (doublereal) (krw - 1);

/*     Check stability. */

    if (lsame_(jobn, "H", (ftnlen)1, (ftnlen)1) && nr > 0) {
	*info = 6;
	return ret_val;
    }

    if (*nq > 0) {
	ku = 1;
	mxnp = max(*nq,*p);
	ktau = *nq * mxnp + 1;
	krw = ktau + min(*nq,*p);

/*        Find X, the solution of Lyapunov equation. */

/*        Workspace needed:      N*MAX(N,P) + MIN(N,P); */
/*        Additional workspace:  4*N; */
/*                               prefer larger. */

	dlacpy_("Full", p, nq, &c__[c_offset], ldc, &dwork[ku], &mxnp, (
		ftnlen)4);
	i__1 = *ldwork - krw + 1;
	sb03ou_(&discr, &c_false, nq, p, &a[a_offset], lda, &dwork[ku], &mxnp,
		 &dwork[ktau], &dwork[ku], nq, &scale, &dwork[krw], &i__1, 
		info);
	if (*info != 0) {
	    if (*info == 1) {
		*info = 4;
	    } else if (*info == 2) {
		*info = 3;
	    }
	    return ret_val;
	}

/* Computing MAX */
	d__1 = wrkopt, d__2 = dwork[krw] + (doublereal) (krw - 1);
	wrkopt = max(d__1,d__2);

/*        Add the contribution of BQ'*X*BQ. */

/*        Workspace needed:      N*(N+M). */

	ktau = *nq * *nq + 1;
	dlacpy_("Full", nq, m, &b[b_offset], ldb, &dwork[ktau], nq, (ftnlen)4)
		;
	dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", nq, m, &c_b16, &
		dwork[ku], nq, &dwork[ktau], nq, (ftnlen)4, (ftnlen)5, (
		ftnlen)11, (ftnlen)7);
	if (nr > 0) {
	    s2norm = dlange_("Frobenius", p, m, &d__[d_offset], ldd, &dwork[1]
		    , (ftnlen)9);
	}
	d__1 = dlange_("Frobenius", nq, m, &dwork[ktau], nq, &dwork[1], (
		ftnlen)9) / scale;
	s2norm = dlapy2_(&s2norm, &d__1);
    }

    ret_val = s2norm;

    dwork[1] = wrkopt;

    return ret_val;
/* *** Last line of AB13BD *** */
} /* ab13bd_ */

