/* SG03AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b20 = 0.;
static doublereal c_b21 = 1.;
static integer c__1 = 1;

/* Subroutine */ int sg03ad_(char *dico, char *job, char *fact, char *trans, 
	char *uplo, integer *n, doublereal *a, integer *lda, doublereal *e, 
	integer *lde, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, doublereal *x, integer *ldx, doublereal *scale, doublereal *sep, 
	doublereal *ferr, doublereal *alphar, doublereal *alphai, doublereal *
	beta, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen 
	trans_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, x_dim1, 
	    x_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal eps, est;
    static integer kase, info1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb01rd_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dgegs_(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sg03ax_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), sg03ay_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), mb01rw_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal norma;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal norme;
    static logical wantx;
    static doublereal scale1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    static logical isfact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical isdisc, wantbh;
    static char etrans[1];
    static logical istran;
    static integer minwrk;
    static logical wantsp, isuppr;
    static integer optwrk;


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

/*     To solve for X either the generalized continuous-time Lyapunov */
/*     equation */

/*             T                T */
/*        op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1) */

/*     or the generalized discrete-time Lyapunov equation */

/*             T                T */
/*        op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2) */

/*     where op(M) is either M or M**T for M = A, E and the right hand */
/*     side Y is symmetric. A, E, Y, and the solution X are N-by-N */
/*     matrices. SCALE is an output scale factor, set to avoid overflow */
/*     in X. */

/*     Estimates of the separation and the relative forward error norm */
/*     are provided. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies which type of the equation is considered: */
/*             = 'C':  Continuous-time equation (1); */
/*             = 'D':  Discrete-time equation (2). */

/*     JOB     CHARACTER*1 */
/*             Specifies if the solution is to be computed and if the */
/*             separation is to be estimated: */
/*             = 'X':  Compute the solution only; */
/*             = 'S':  Estimate the separation only; */
/*             = 'B':  Compute the solution and estimate the separation. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether the generalized real Schur */
/*             factorization of the pencil A - lambda * E is supplied */
/*             on entry or not: */
/*             = 'N':  Factorization is not supplied; */
/*             = 'F':  Factorization is supplied. */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  op(A) = A,    op(E) = E; */
/*             = 'T':  op(A) = A**T, op(E) = E**T. */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the lower or the upper triangle of the */
/*             array X is needed on input: */
/*             = 'L':  Only the lower triangle is needed on input; */
/*             = 'U':  Only the upper triangle is needed on input. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N upper */
/*             Hessenberg part of this array must contain the */
/*             generalized Schur factor A_s of the matrix A (see */
/*             definition (3) in section METHOD). A_s must be an upper */
/*             quasitriangular matrix. The elements below the upper */
/*             Hessenberg part of the array A are not referenced. */
/*             If FACT = 'N', then the leading N-by-N part of this */
/*             array must contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the generalized Schur factor A_s of the matrix A. (A_s is */
/*             an upper quasitriangular matrix.) */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N upper */
/*             triangular part of this array must contain the */
/*             generalized Schur factor E_s of the matrix E (see */
/*             definition (4) in section METHOD). The elements below the */
/*             upper triangular part of the array E are not referenced. */
/*             If FACT = 'N', then the leading N-by-N part of this */
/*             array must contain the coefficient matrix E of the */
/*             equation. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the generalized Schur factor E_s of the matrix E. (E_s is */
/*             an upper triangular matrix.) */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N part of */
/*             this array must contain the orthogonal matrix Q from */
/*             the generalized Schur factorization (see definitions (3) */
/*             and (4) in section METHOD). */
/*             If FACT = 'N', Q need not be set on entry. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal matrix Q from the generalized Schur */
/*             factorization. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N part of */
/*             this array must contain the orthogonal matrix Z from */
/*             the generalized Schur factorization (see definitions (3) */
/*             and (4) in section METHOD). */
/*             If FACT = 'N', Z need not be set on entry. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal matrix Z from the generalized Schur */
/*             factorization. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if JOB = 'B' or 'X', then the leading N-by-N */
/*             part of this array must contain the right hand side matrix */
/*             Y of the equation. Either the lower or the upper */
/*             triangular part of this array is needed (see mode */
/*             parameter UPLO). */
/*             If JOB = 'S', X is not referenced. */
/*             On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then */
/*             the leading N-by-N part of this array contains the */
/*             solution matrix X of the equation. */
/*             If JOB = 'S', X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             (0 < SCALE <= 1) */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then */
/*             SEP contains an estimate of the separation of the */
/*             Lyapunov operator. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an */
/*             estimated forward error bound for the solution X. If XTRUE */
/*             is the true solution, FERR estimates the relative error */
/*             in the computed solution, measured in the Frobenius norm: */
/*             norm(X - XTRUE) / norm(XTRUE) */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             If FACT = 'N' and INFO = 0, 3, or 4, then */
/*             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the */
/*             eigenvalues of the matrix pencil A - lambda * E. */
/*             If FACT = 'F', ALPHAR, ALPHAI, and BETA are not */
/*             referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N**2) */
/*             IWORK is not referenced if JOB = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. The following table */
/*             contains the minimal work space requirements depending */
/*             on the choice of JOB and FACT. */

/*                    JOB        FACT    |  LDWORK */
/*                    -------------------+------------------- */
/*                    'X'        'F'     |  MAX(1,N) */
/*                    'X'        'N'     |  MAX(1,4*N) */
/*                    'B', 'S'   'F'     |  MAX(1,2*N**2) */
/*                    'B', 'S'   'N'     |  MAX(1,2*N**2,4*N) */

/*             For optimum performance, LDWORK should be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  FACT = 'F' and the matrix contained in the upper */
/*                   Hessenberg part of the array A is not in upper */
/*                   quasitriangular form; */
/*             = 2:  FACT = 'N' and the pencil A - lambda * E cannot be */
/*                   reduced to generalized Schur form: LAPACK routine */
/*                   DGEGS has failed to converge; */
/*             = 3:  DICO = 'D' and the pencil A - lambda * E has a */
/*                   pair of reciprocal eigenvalues. That is, lambda_i = */
/*                   1/lambda_j for some i and j, where lambda_i and */
/*                   lambda_j are eigenvalues of A - lambda * E. Hence, */
/*                   equation (2) is singular;  perturbed values were */
/*                   used to solve the equation (but the matrices A and */
/*                   E are unchanged); */
/*             = 4:  DICO = 'C' and the pencil A - lambda * E has a */
/*                   degenerate pair of eigenvalues. That is, lambda_i = */
/*                   -lambda_j for some i and j, where lambda_i and */
/*                   lambda_j are eigenvalues of A - lambda * E. Hence, */
/*                   equation (1) is singular;  perturbed values were */
/*                   used to solve the equation (but the matrices A and */
/*                   E are unchanged). */

/*     METHOD */

/*     A straightforward generalization [3] of the method proposed by */
/*     Bartels and Stewart [1] is utilized to solve (1) or (2). */

/*     First the pencil A - lambda * E is reduced to real generalized */
/*     Schur form A_s - lambda * E_s by means of orthogonal */
/*     transformations (QZ-algorithm): */

/*        A_s = Q**T * A * Z   (upper quasitriangular)                (3) */

/*        E_s = Q**T * E * Z   (upper triangular).                    (4) */

/*     If FACT = 'F', this step is omitted. Assuming SCALE = 1 and */
/*     defining */

/*              ( Z**T * Y * Z   :   TRANS = 'N' */
/*        Y_s = < */
/*              ( Q**T * Y * Q   :   TRANS = 'T' */


/*              ( Q**T * X * Q    if TRANS = 'N' */
/*        X_s = <                                                     (5) */
/*              ( Z**T * X * Z    if TRANS = 'T' */

/*     leads to the reduced Lyapunov equation */

/*               T                      T */
/*        op(A_s)  X_s op(E_s) + op(E_s)  X_s op(A_s) = Y_s,          (6) */

/*     or */
/*               T                      T */
/*        op(A_s)  X_s op(A_s) - op(E_s)  X_s op(E_s) = Y_s,          (7) */

/*     which are equivalent to (1) or (2), respectively. The solution X_s */
/*     of (6) or (7) is computed via block back substitution (if TRANS = */
/*     'N') or block forward substitution (if TRANS = 'T'), where the */
/*     block order is at most 2. (See [1] and [3] for details.) */
/*     Equation (5) yields the solution matrix X. */

/*     For fast computation the estimates of the separation and the */
/*     forward error are gained from (6) or (7) rather than (1) or */
/*     (2), respectively. We consider (6) and (7) as special cases of the */
/*     generalized Sylvester equation */

/*        R * X * S + U * X * V = Y,                                  (8) */

/*     whose separation is defined as follows */

/*        sep = sep(R,S,U,V) =   min   || R * X * S + U * X * V || . */
/*                            ||X|| = 1                           F */
/*                                 F */

/*     Equation (8) is equivalent to the system of linear equations */

/*        K * vec(X) = (kron(S**T,R) + kron(V**T,U)) * vec(X) = vec(Y), */

/*     where kron is the Kronecker product of two matrices and vec */
/*     is the mapping that stacks the columns of a matrix. If K is */
/*     nonsingular then */

/*        sep = 1 / ||K**(-1)|| . */
/*                             2 */

/*     We estimate ||K**(-1)|| by a method devised by Higham [2]. Note */
/*     that this method yields an estimation for the 1-norm but we use it */
/*     as an approximation for the 2-norm. Estimates for the forward */
/*     error norm are provided by */

/*        FERR = 2 * EPS * ||A_s||  * ||E_s||  / sep */
/*                                F          F */

/*     in the continuous-time case (1) and */

/*        FERR = EPS * ( ||A_s|| **2 + ||E_s|| **2 ) / sep */
/*                              F             F */

/*     in the discrete-time case (2). */
/*     The reciprocal condition number, RCOND, of the Lyapunov equation */
/*     can be estimated by FERR/EPS. */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or complex */
/*         matrix, with applications to condition estimation. */
/*         A.C.M. Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, 1988. */

/*     [3] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The number of flops required by the routine is given by the */
/*     following table. Note that we count a single floating point */
/*     arithmetic operation as one flop. c is an integer number of modest */
/*     size (say 4 or 5). */

/*                   |  FACT = 'F'            FACT = 'N' */
/*        -----------+------------------------------------------ */
/*        JOB = 'B'  |  (26+8*c)/3 * N**3     (224+8*c)/3 * N**3 */
/*        JOB = 'S'  |  8*c/3 * N**3          (198+8*c)/3 * N**3 */
/*        JOB = 'X'  |  26/3 * N**3           224/3 * N**3 */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if DICO = 'D' and the pencil A - lambda * E has a pair of almost */
/*     reciprocal eigenvalues, or DICO = 'C' and the pencil has an almost */
/*     degenerate pair of eigenvalues, then the Lyapunov equation will be */
/*     ill-conditioned. Perturbed values were used to solve the equation. */
/*     Ill-conditioning can be detected by a very small value of the */
/*     reciprocal condition number RCOND. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --alphar;
    --alphai;
    --beta;
    --iwork;
    --dwork;

    /* Function Body */
    isdisc = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    wantx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
    wantsp = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
    isfact = lsame_(fact, "F", (ftnlen)1, (ftnlen)1);
    istran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    isuppr = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (isdisc || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (wantx || wantsp || wantbh)) {
	*info = -2;
    } else if (! (isfact || lsame_(fact, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (istran || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (! (isuppr || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*lde < max(1,*n)) {
	*info = -10;
    } else if (*ldq < max(1,*n)) {
	*info = -12;
    } else if (*ldz < max(1,*n)) {
	*info = -14;
    } else if (*ldx < max(1,*n)) {
	*info = -16;
    } else {
	*info = 0;
    }
    if (*info == 0) {

/*        Compute minimal workspace. */

	if (wantx) {
	    if (isfact) {
		minwrk = max(*n,1);
	    } else {
/* Computing MAX */
		i__1 = *n << 2;
		minwrk = max(i__1,1);
	    }
	} else {
	    if (isfact) {
/* Computing MAX */
		i__1 = (*n << 1) * *n;
		minwrk = max(i__1,1);
	    } else {
/* Computing MAX */
		i__1 = (*n << 1) * *n, i__2 = *n << 2, i__1 = max(i__1,i__2);
		minwrk = max(i__1,1);
	    }
	}
	if (minwrk > *ldwork) {
	    *info = -25;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SG03AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*scale = 1.;
	if (! wantx) {
	    *sep = 0.;
	}
	if (wantbh) {
	    *ferr = 0.;
	}
	dwork[1] = 1.;
	return 0;
    }

    if (isfact) {

/*        Make sure the upper Hessenberg part of A is quasitriangular. */

	i__1 = *n - 2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (a[i__ + 1 + i__ * a_dim1] != 0. && a[i__ + 2 + (i__ + 1) * 
		    a_dim1] != 0.) {
		*info = 1;
		return 0;
	    }
/* L20: */
	}
    }

    if (! isfact) {

/*        Reduce A - lambda * E to generalized Schur form. */

/*           A := Q**T * A * Z   (upper quasitriangular) */
/*           E := Q**T * E * Z   (upper triangular) */

/*        ( Workspace: >= MAX(1,4*N) ) */

	dgegs_("Vectors", "Vectors", n, &a[a_offset], lda, &e[e_offset], lde, 
		&alphar[1], &alphai[1], &beta[1], &q[q_offset], ldq, &z__[
		z_offset], ldz, &dwork[1], ldwork, &info1, (ftnlen)7, (ftnlen)
		7);
	if (info1 != 0) {
	    *info = 2;
	    return 0;
	}
	optwrk = (integer) dwork[1];
    } else {
	optwrk = minwrk;
    }

    if (wantbh || wantx) {

/*        Transform right hand side. */

/*           X := Z**T * X * Z  or  X := Q**T * X * Q */

/*        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*        ( Workspace: >= N ) */

	if (*ldwork < *n * *n) {
	    if (istran) {
		mb01rw_(uplo, "Transpose", n, n, &x[x_offset], ldx, &q[
			q_offset], ldq, &dwork[1], &info1, (ftnlen)1, (ftnlen)
			9);
	    } else {
		mb01rw_(uplo, "Transpose", n, n, &x[x_offset], ldx, &z__[
			z_offset], ldz, &dwork[1], &info1, (ftnlen)1, (ftnlen)
			9);
	    }
	} else {
	    if (istran) {
		mb01rd_(uplo, "Transpose", n, n, &c_b20, &c_b21, &x[x_offset],
			 ldx, &q[q_offset], ldq, &x[x_offset], ldx, &dwork[1],
			 ldwork, info, (ftnlen)1, (ftnlen)9);
	    } else {
		mb01rd_(uplo, "Transpose", n, n, &c_b20, &c_b21, &x[x_offset],
			 ldx, &z__[z_offset], ldz, &x[x_offset], ldx, &dwork[
			1], ldwork, info, (ftnlen)1, (ftnlen)9);
	    }
	}
	if (! isuppr) {
	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - i__;
		dcopy_(&i__2, &x[i__ + 1 + i__ * x_dim1], &c__1, &x[i__ + (
			i__ + 1) * x_dim1], ldx);
/* L40: */
	    }
	}
/* Computing MAX */
	i__1 = optwrk, i__2 = *n * *n;
	optwrk = max(i__1,i__2);

/*        Solve reduced generalized Lyapunov equation. */

	if (isdisc) {
	    sg03ax_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &x[
		    x_offset], ldx, scale, &info1, (ftnlen)1);
	    if (info1 != 0) {
		*info = 3;
	    }
	} else {
	    sg03ay_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &x[
		    x_offset], ldx, scale, &info1, (ftnlen)1);
	    if (info1 != 0) {
		*info = 4;
	    }
	}

/*        Transform the solution matrix back. */

/*           X := Q * X * Q**T  or  X := Z * X * Z**T. */

/*        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*        ( Workspace: >= N ) */

	if (*ldwork < *n * *n) {
	    if (istran) {
		mb01rw_("Upper", "NoTranspose", n, n, &x[x_offset], ldx, &z__[
			z_offset], ldz, &dwork[1], &info1, (ftnlen)5, (ftnlen)
			11);
	    } else {
		mb01rw_("Upper", "NoTranspose", n, n, &x[x_offset], ldx, &q[
			q_offset], ldq, &dwork[1], &info1, (ftnlen)5, (ftnlen)
			11);
	    }
	} else {
	    if (istran) {
		mb01rd_("Upper", "NoTranspose", n, n, &c_b20, &c_b21, &x[
			x_offset], ldx, &z__[z_offset], ldz, &x[x_offset], 
			ldx, &dwork[1], ldwork, info, (ftnlen)5, (ftnlen)11);
	    } else {
		mb01rd_("Upper", "NoTranspose", n, n, &c_b20, &c_b21, &x[
			x_offset], ldx, &q[q_offset], ldq, &x[x_offset], ldx, 
			&dwork[1], ldwork, info, (ftnlen)5, (ftnlen)11);
	    }
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__;
	    dcopy_(&i__2, &x[i__ + (i__ + 1) * x_dim1], ldx, &x[i__ + 1 + i__ 
		    * x_dim1], &c__1);
/* L60: */
	}
    }

    if (wantbh || wantsp) {

/*        Estimate the 1-norm of the inverse Kronecker product matrix */
/*        belonging to the reduced generalized Lyapunov equation. */

/*        ( Workspace: 2*N*N ) */

	est = 0.;
	kase = 0;
L80:
	i__1 = *n * *n;
	dlacon_(&i__1, &dwork[*n * *n + 1], &dwork[1], &iwork[1], &est, &kase)
		;
	if (kase != 0) {
	    if (kase == 1 && ! istran || kase != 1 && istran) {
		*(unsigned char *)etrans = 'N';
	    } else {
		*(unsigned char *)etrans = 'T';
	    }
	    if (isdisc) {
		sg03ax_(etrans, n, &a[a_offset], lda, &e[e_offset], lde, &
			dwork[1], n, &scale1, &info1, (ftnlen)1);
		if (info1 != 0) {
		    *info = 3;
		}
	    } else {
		sg03ay_(etrans, n, &a[a_offset], lda, &e[e_offset], lde, &
			dwork[1], n, &scale1, &info1, (ftnlen)1);
		if (info1 != 0) {
		    *info = 4;
		}
	    }
	    goto L80;
	}
	*sep = scale1 / est;
    }

/*     Estimate the relative forward error. */

/*     ( Workspace: 2*N ) */

    if (wantbh) {
	eps = dlamch_("Precision", (ftnlen)9);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	    i__3 = i__ + 1;
	    i__2 = min(i__3,*n);
	    dwork[i__] = dnrm2_(&i__2, &a[i__ * a_dim1 + 1], &c__1);
	    dwork[*n + i__] = dnrm2_(&i__, &e[i__ * e_dim1 + 1], &c__1);
/* L100: */
	}
	norma = dnrm2_(n, &dwork[1], &c__1);
	norme = dnrm2_(n, &dwork[*n + 1], &c__1);
	if (isdisc) {
/* Computing 2nd power */
	    d__1 = norma;
/* Computing 2nd power */
	    d__2 = norme;
	    *ferr = (d__1 * d__1 + d__2 * d__2) * eps / *sep;
	} else {
	    *ferr = norma * 2. * norme * eps / *sep;
	}
    }

    dwork[1] = (doublereal) max(optwrk,minwrk);
    return 0;
/* *** Last line of SG03AD *** */
} /* sg03ad_ */

