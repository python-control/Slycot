/* SB02MD.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b42 = 0.;
static doublereal c_b45 = 1.;
static doublereal c_b47 = .5;

/* Subroutine */ int sb02md_(char *dico, char *hinv, char *uplo, char *scal, 
	char *sort, integer *n, doublereal *a, integer *lda, doublereal *g, 
	integer *ldg, doublereal *q, integer *ldq, doublereal *rcond, 
	doublereal *wr, doublereal *wi, doublereal *s, integer *lds, 
	doublereal *u, integer *ldu, integer *iwork, doublereal *dwork, 
	integer *ldwork, logical *bwork, integer *info, ftnlen dico_len, 
	ftnlen hinv_len, ftnlen uplo_len, ftnlen scal_len, ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, s_dim1, 
	    s_offset, u_dim1, u_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, n2, np1, iscl, ierr, nrot;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgees_(char *, char *, L_fp, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern logical sb02mr_(), sb02ms_();
    extern /* Subroutine */ int sb02mu_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern logical sb02mv_(), sb02mw_();
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lhinv;
    static doublereal gnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical luplo;
    static doublereal qnorm, unorm;
    static logical lsort;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal rconda;
    extern /* Subroutine */ int dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
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

/*     To solve for X either the continuous-time algebraic Riccati */
/*     equation */
/*                              -1 */
/*        Q + A'*X + X*A - X*B*R  B'*X = 0                            (1) */

/*     or the discrete-time algebraic Riccati equation */
/*                                        -1 */
/*        X = A'*X*A - A'*X*B*(R + B'*X*B)  B'*X*A + Q                (2) */

/*     where A, B, Q and R are N-by-N, N-by-M, N-by-N and M-by-M matrices */
/*     respectively, with Q symmetric and R symmetric nonsingular; X is */
/*     an N-by-N symmetric matrix. */
/*                       -1 */
/*     The matrix G = B*R  B' must be provided on input, instead of B and */
/*     R, that is, for instance, the continuous-time equation */

/*        Q + A'*X + X*A - X*G*X = 0                                  (3) */

/*     is solved, where G is an N-by-N symmetric matrix. SLICOT Library */
/*     routine SB02MT should be used to compute G, given B and R. SB02MT */
/*     also enables to solve Riccati equations corresponding to optimal */
/*     problems with coupling terms. */

/*     The routine also returns the computed values of the closed-loop */
/*     spectrum of the optimal system, i.e., the stable eigenvalues */
/*     lambda(1),...,lambda(N) of the corresponding Hamiltonian or */
/*     symplectic matrix associated to the optimal problem. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of Riccati equation to be solved as */
/*             follows: */
/*             = 'C':  Equation (3), continuous-time case; */
/*             = 'D':  Equation (2), discrete-time case. */

/*     HINV    CHARACTER*1 */
/*             If DICO = 'D', specifies which symplectic matrix is to be */
/*             constructed, as follows: */
/*             = 'D':  The matrix H in (5) (see METHOD) is constructed; */
/*             = 'I':  The inverse of the matrix H in (5) is constructed. */
/*             HINV is not used if DICO = 'C'. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     SCAL    CHARACTER*1 */
/*             Specifies whether or not a scaling strategy should be */
/*             used, as follows: */
/*             = 'G':  General scaling should be used; */
/*             = 'N':  No scaling should be used. */

/*     SORT    CHARACTER*1 */
/*             Specifies which eigenvalues should be obtained in the top */
/*             of the Schur form, as follows: */
/*             = 'S':  Stable   eigenvalues come first; */
/*             = 'U':  Unstable eigenvalues come first. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, Q, G and X.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix A of the equation. */
/*             On exit, if DICO = 'D', and INFO = 0 or INFO > 1, the */
/*                                                                    -1 */
/*             leading N-by-N part of this array contains the matrix A  . */
/*             Otherwise, the array A is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             The leading N-by-N upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             must contain the upper triangular part or lower triangular */
/*             part, respectively, of the symmetric matrix G. The stricly */
/*             lower triangular part (if UPLO = 'U') or stricly upper */
/*             triangular part (if UPLO = 'L') is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U') or lower triangular part (if UPLO = 'L') of */
/*             this array must contain the upper triangular part or lower */
/*             triangular part, respectively, of the symmetric matrix Q. */
/*             The stricly lower triangular part (if UPLO = 'U') or */
/*             stricly upper triangular part (if UPLO = 'L') is not used. */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains the solution matrix X of the problem. */

/*     LDQ     INTEGER */
/*             The leading dimension of array N.  LDQ >= MAX(1,N). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             An estimate of the reciprocal of the condition number (in */
/*             the 1-norm) of the N-th order system of algebraic */
/*             equations from which the solution matrix X is obtained. */

/*     WR      (output) DOUBLE PRECISION array, dimension (2*N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (2*N) */
/*             If INFO = 0 or INFO = 5, these arrays contain the real and */
/*             imaginary parts, respectively, of the eigenvalues of the */
/*             2N-by-2N matrix S, ordered as specified by SORT (except */
/*             for the case HINV = 'D', when the order is opposite to */
/*             that specified by SORT). The leading N elements of these */
/*             arrays contain the closed-loop spectrum of the system */
/*                           -1 */
/*             matrix A - B*R  *B'*X, if DICO = 'C', or of the matrix */
/*                               -1 */
/*             A - B*(R + B'*X*B)  B'*X*A, if DICO = 'D'. Specifically, */
/*                lambda(k) = WR(k) + j*WI(k), for k = 1,2,...,N. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N) */
/*             If INFO = 0 or INFO = 5, the leading 2N-by-2N part of this */
/*             array contains the ordered real Schur form S of the */
/*             Hamiltonian or symplectic matrix H. That is, */

/*                    (S   S  ) */
/*                    ( 11  12) */
/*                S = (       ), */
/*                    (0   S  ) */
/*                    (     22) */

/*             where S  , S   and S   are N-by-N matrices. */
/*                    11   12      22 */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,2*N). */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N) */
/*             If INFO = 0 or INFO = 5, the leading 2N-by-2N part of this */
/*             array contains the transformation matrix U which reduces */
/*             the Hamiltonian or symplectic matrix H to the ordered real */
/*             Schur form S. That is, */

/*                    (U   U  ) */
/*                    ( 11  12) */
/*                U = (       ), */
/*                    (U   U  ) */
/*                    ( 21  22) */

/*             where U  , U  , U   and U   are N-by-N matrices. */
/*                    11   12   21      22 */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,2*N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) returns the scaling factor used */
/*             (set to 1 if SCAL = 'N'), also set if INFO = 5; */
/*             if DICO = 'D', DWORK(3) returns the reciprocal condition */
/*             number of the given matrix  A. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(2,6*N) if DICO = 'C'; */
/*             LDWORK >= MAX(3,6*N) if DICO = 'D'. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if matrix A is (numerically) singular in discrete- */
/*                   time case; */
/*             = 2:  if the Hamiltonian or symplectic matrix H cannot be */
/*                   reduced to real Schur form; */
/*             = 3:  if the real Schur form of the Hamiltonian or */
/*                   symplectic matrix H cannot be appropriately ordered; */
/*             = 4:  if the Hamiltonian or symplectic matrix H has less */
/*                   than N stable eigenvalues; */
/*             = 5:  if the N-th order system of linear algebraic */
/*                   equations, from which the solution matrix X would */
/*                   be obtained, is singular to working precision. */

/*     METHOD */

/*     The method used is the Schur vector approach proposed by Laub. */
/*     It is assumed that [A,B] is a stabilizable pair (where for (3) B */
/*     is any matrix such that B*B' = G with rank(B) = rank(G)), and */
/*     [E,A] is a detectable pair, where E is any matrix such that */
/*     E*E' = Q with rank(E) = rank(Q). Under these assumptions, any of */
/*     the algebraic Riccati equations (1)-(3) is known to have a unique */
/*     non-negative definite solution. See [2]. */
/*     Now consider the 2N-by-2N Hamiltonian or symplectic matrix */

/*                 ( A   -G ) */
/*            H =  (        ),                                    (4) */
/*                 (-Q   -A'), */

/*     for continuous-time equation, and */
/*                    -1        -1 */
/*                 ( A         A  *G   ) */
/*            H =  (   -1          -1  ),                         (5) */
/*                 (Q*A    A' + Q*A  *G) */
/*                                                            -1 */
/*     for discrete-time equation, respectively, where G = B*R  *B'. */
/*     The assumptions guarantee that H in (4) has no pure imaginary */
/*     eigenvalues, and H in (5) has no eigenvalues on the unit circle. */
/*     If Y is an N-by-N matrix then there exists an orthogonal matrix U */
/*     such that U'*Y*U is an upper quasi-triangular matrix. Moreover, U */
/*     can be chosen so that the 2-by-2 and 1-by-1 diagonal blocks */
/*     (corresponding to the complex conjugate eigenvalues and real */
/*     eigenvalues respectively) appear in any desired order. This is the */
/*     ordered real Schur form. Thus, we can find an orthogonal */
/*     similarity transformation U which puts (4) or (5) in ordered real */
/*     Schur form */

/*            U'*H*U = S = (S(1,1)  S(1,2)) */
/*                         (  0     S(2,2)) */

/*     where S(i,j) is an N-by-N matrix and the eigenvalues of S(1,1) */
/*     have negative real parts in case of (4), or moduli greater than */
/*     one in case of (5). If U is conformably partitioned into four */
/*     N-by-N blocks */

/*               U = (U(1,1)  U(1,2)) */
/*                   (U(2,1)  U(2,2)) */

/*     with respect to the assumptions we then have */
/*     (a) U(1,1) is invertible and X = U(2,1)*inv(U(1,1)) solves (1), */
/*         (2), or (3) with X = X' and non-negative definite; */
/*     (b) the eigenvalues of S(1,1) (if DICO = 'C') or S(2,2) (if */
/*         DICO = 'D') are equal to the eigenvalues of optimal system */
/*         (the 'closed-loop' spectrum). */

/*     [A,B] is stabilizable if there exists a matrix F such that (A-BF) */
/*     is stable. [E,A] is detectable if [A',E'] is stabilizable. */

/*     REFERENCES */

/*     [1] Laub, A.J. */
/*         A Schur Method for Solving Algebraic Riccati equations. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 913-921, 1979. */

/*     [2] Wonham, W.M. */
/*         On a matrix Riccati equation of stochastic control. */
/*         SIAM J. Contr., 6, pp. 681-697, 1968. */

/*     [3] Sima, V. */
/*         Algorithms for Linear-Quadratic Optimization. */
/*         Pure and Applied Mathematics: A Series of Monographs and */
/*         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     To obtain a stabilizing solution of the algebraic Riccati */
/*     equation for DICO = 'D', set SORT = 'U', if HINV = 'D', or set */
/*     SORT = 'S', if HINV = 'I'. */

/*     The routine can also compute the anti-stabilizing solutions of */
/*     the algebraic Riccati equations, by specifying */
/*         SORT = 'U' if DICO = 'D' and HINV = 'I', or DICO = 'C', or */
/*         SORT = 'S' if DICO = 'D' and HINV = 'D'. */

/*     Usually, the combinations HINV = 'D' and SORT = 'U', or HINV = 'I' */
/*     and SORT = 'U', will be faster then the other combinations [3]. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB02AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

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
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --wr;
    --wi;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    *info = 0;
    n2 = *n + *n;
    np1 = *n + 1;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    lscal = lsame_(scal, "G", (ftnlen)1, (ftnlen)1);
    lsort = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    if (discr) {
	lhinv = lsame_(hinv, "D", (ftnlen)1, (ftnlen)1);
    }

/*     Test the input scalar arguments. */

    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (discr) {
	if (! lhinv && ! lsame_(hinv, "I", (ftnlen)1, (ftnlen)1)) {
	    *info = -2;
	}
    }
    if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (! lscal && ! lsame_(scal, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -4;
    } else if (! lsort && ! lsame_(sort, "U", (ftnlen)1, (ftnlen)1)) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldg < max(1,*n)) {
	*info = -10;
    } else if (*ldq < max(1,*n)) {
	*info = -12;
    } else if (*lds < max(1,n2)) {
	*info = -17;
    } else if (*ldu < max(1,n2)) {
	*info = -19;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 2, i__2 = *n * 6;
/* Computing MAX */
	i__3 = 3, i__4 = *n * 6;
	if (! discr && *ldwork < max(i__1,i__2) || discr && *ldwork < max(
		i__3,i__4)) {
	    *info = -22;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB02MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*rcond = 1.;
	dwork[1] = 1.;
	dwork[2] = 1.;
	if (discr) {
	    dwork[3] = 1.;
	}
	return 0;
    }

    if (lscal) {

/*        Compute the norms of the matrices Q and G. */

	qnorm = dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[1], (
		ftnlen)6, (ftnlen)1);
	gnorm = dlansy_("1-norm", uplo, n, &g[g_offset], ldg, &dwork[1], (
		ftnlen)6, (ftnlen)1);
    }

/*     Initialise the Hamiltonian or symplectic matrix associated with */
/*     the problem. */
/*     Workspace:  need   1          if DICO = 'C'; */
/*                        max(2,4*N) if DICO = 'D'; */
/*                 prefer larger if DICO = 'D'. */

    sb02mu_(dico, hinv, uplo, n, &a[a_offset], lda, &g[g_offset], ldg, &q[
	    q_offset], ldq, &s[s_offset], lds, &iwork[1], &dwork[1], ldwork, 
	    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (*info != 0) {
	*info = 1;
	return 0;
    }

    wrkopt = dwork[1];
    if (discr) {
	rconda = dwork[2];
    }

    iscl = 0;
    if (lscal) {

/*        Scale the Hamiltonian or symplectic matrix. */

	if (qnorm > gnorm && gnorm > 0.) {
	    dlascl_("G", &c__0, &c__0, &qnorm, &gnorm, n, n, &s[np1 + s_dim1],
		     &n2, &ierr, (ftnlen)1);
	    dlascl_("G", &c__0, &c__0, &gnorm, &qnorm, n, n, &s[np1 * s_dim1 
		    + 1], &n2, &ierr, (ftnlen)1);
	    iscl = 1;
	}
    }

/*     Find the ordered Schur factorization of S,   S = U*H*U'. */
/*     Workspace:  need   6*N; */
/*                 prefer larger. */

    if (! discr) {
	if (lsort) {
	    dgees_("Vectors", "Sorted", (L_fp)sb02mv_, &n2, &s[s_offset], lds,
		     &nrot, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], 
		    ldwork, &bwork[1], info, (ftnlen)7, (ftnlen)6);
	} else {
	    dgees_("Vectors", "Sorted", (L_fp)sb02mr_, &n2, &s[s_offset], lds,
		     &nrot, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], 
		    ldwork, &bwork[1], info, (ftnlen)7, (ftnlen)6);
	}
    } else {
	if (lsort) {
	    dgees_("Vectors", "Sorted", (L_fp)sb02mw_, &n2, &s[s_offset], lds,
		     &nrot, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], 
		    ldwork, &bwork[1], info, (ftnlen)7, (ftnlen)6);
	} else {
	    dgees_("Vectors", "Sorted", (L_fp)sb02ms_, &n2, &s[s_offset], lds,
		     &nrot, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], 
		    ldwork, &bwork[1], info, (ftnlen)7, (ftnlen)6);
	}
	if (lhinv) {
	    dswap_(n, &wr[1], &c__1, &wr[np1], &c__1);
	    dswap_(n, &wi[1], &c__1, &wi[np1], &c__1);
	}
    }
    if (*info > n2) {
	*info = 3;
    } else if (*info > 0) {
	*info = 2;
    } else if (nrot != *n) {
	*info = 4;
    }
    if (*info != 0) {
	return 0;
    }

    wrkopt = max(wrkopt,dwork[1]);

/*     Check if U(1,1) is singular.  Use the (2,1) block of S as a */
/*     workspace for factoring U(1,1). */

    unorm = dlange_("1-norm", n, n, &u[u_offset], ldu, &dwork[1], (ftnlen)6);

    dlacpy_("Full", n, n, &u[u_offset], ldu, &s[np1 + s_dim1], lds, (ftnlen)4)
	    ;
    dgetrf_(n, n, &s[np1 + s_dim1], lds, &iwork[1], info);

    if (*info > 0) {

/*        Singular matrix.  Set INFO and RCOND for error return. */

	*info = 5;
	*rcond = 0.;
	goto L100;
    }

/*     Estimate the reciprocal condition of U(1,1). */
/*     Workspace: 6*N. */

    dgecon_("1-norm", n, &s[np1 + s_dim1], lds, &unorm, rcond, &dwork[1], &
	    iwork[np1], info, (ftnlen)6);

    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {

/*        Nearly singular matrix.  Set INFO for error return. */

	*info = 5;
	return 0;
    }

/*     Transpose U(2,1) in Q and compute the solution. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(n, &u[np1 + i__ * u_dim1], &c__1, &q[i__ + q_dim1], ldq);
/* L60: */
    }

    dgetrs_("Transpose", n, n, &s[np1 + s_dim1], lds, &iwork[1], &q[q_offset],
	     ldq, info, (ftnlen)9);

/*     Set S(2,1) to zero. */

    dlaset_("Full", n, n, &c_b42, &c_b42, &s[np1 + s_dim1], lds, (ftnlen)4);

/*     Make sure the solution matrix X is symmetric. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__;
	daxpy_(&i__2, &c_b45, &q[i__ + (i__ + 1) * q_dim1], ldq, &q[i__ + 1 + 
		i__ * q_dim1], &c__1);
	i__2 = *n - i__;
	dscal_(&i__2, &c_b47, &q[i__ + 1 + i__ * q_dim1], &c__1);
	i__2 = *n - i__;
	dcopy_(&i__2, &q[i__ + 1 + i__ * q_dim1], &c__1, &q[i__ + (i__ + 1) * 
		q_dim1], ldq);
/* L80: */
    }

    if (lscal) {

/*        Undo scaling for the solution matrix. */

	if (iscl == 1) {
	    dlascl_("G", &c__0, &c__0, &gnorm, &qnorm, n, n, &q[q_offset], 
		    ldq, &ierr, (ftnlen)1);
	}
    }

/*     Set the optimal workspace, the scaling factor, and reciprocal */
/*     condition number (if any). */

    dwork[1] = wrkopt;
L100:
    if (iscl == 1) {
	dwork[2] = qnorm / gnorm;
    } else {
	dwork[2] = 1.;
    }
    if (discr) {
	dwork[3] = rconda;
    }

    return 0;
/* *** Last line of SB02MD *** */
} /* sb02md_ */

