/* SB02PD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static doublereal c_b57 = -.5;
static doublereal c_b66 = 0.;
static doublereal c_b67 = 1.;
static doublereal c_b81 = -1.;

/* Subroutine */ int sb02pd_(char *job, char *trana, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	q, integer *ldq, doublereal *x, integer *ldx, doublereal *rcond, 
	doublereal *ferr, doublereal *wr, doublereal *wi, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen trana_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, x_dim1, 
	    x_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, n2, ib, ic, ij, ji, ir, it, iu, ij1, ij2, iaf;
    static logical all;
    static integer ibr, ini, ifr;
    static doublereal eps, sep, tol;
    static integer isv, iscl, sdim, itau, iter;
    static doublereal conv, temp;
    static integer iwrk;
    static char loup[1];
    static integer info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02ed_(char *, integer *, doublereal *, integer *, ftnlen), 
	    dscal_(integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), sb02qd_(char *, char *, char *, char *
	    , char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char equed[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal hnorm;
    static logical bwork[1], lower;
    extern /* Subroutine */ int dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dgeqp3_(
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static doublereal gnorm2, qnorm2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern logical select_();
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static doublereal hinnrm;
    extern /* Subroutine */ int dgesvx_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, char 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static logical notrna;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dsytri_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);


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

/*     To solve the real continuous-time matrix algebraic Riccati */
/*     equation */

/*        op(A)'*X + X*op(A) + Q - X*G*X = 0, */

/*     where op(A) = A or A' = A**T and G, Q are symmetric (G = G**T, */
/*     Q = Q**T). The matrices A, G and Q are N-by-N and the solution X */
/*     is an N-by-N symmetric matrix. */

/*     An error bound on the solution and a condition estimate are also */
/*     optionally provided. */

/*     It is assumed that the matrices A, G and Q are such that the */
/*     corresponding Hamiltonian matrix has N eigenvalues with negative */
/*     real parts. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'X':  Compute the solution only; */
/*             = 'A':  Compute all: the solution, reciprocal condition */
/*                     number, and the error bound. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the option op(A): */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangles of G and Q are stored; */
/*             = 'L':  Lower triangles of G and Q are stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, Q, and X.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             coefficient matrix A of the equation. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix G. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix G. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= max(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix Q. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= max(1,N). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             If INFO = 0, INFO = 2, or INFO = 4, the leading N-by-N */
/*             part of this array contains the symmetric solution matrix */
/*             X of the algebraic Riccati equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'A', the estimate of the reciprocal condition */
/*             number of the Riccati equation. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'A', the estimated forward error bound for the */
/*             solution X. If XTRUE is the true solution, FERR bounds the */
/*             magnitude of the largest entry in (X - XTRUE) divided by */
/*             the magnitude of the largest entry in X. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             If JOB = 'A' and TRANA = 'N', WR and WI contain the real */
/*             and imaginary parts, respectively, of the eigenvalues of */
/*             the matrix A - G*X, i.e., the closed-loop system poles. */
/*             If JOB = 'A' and TRANA = 'T' or 'C', WR and WI contain the */
/*             real and imaginary parts, respectively, of the eigenvalues */
/*             of the matrix A - X*G, i.e., the closed-loop system poles. */
/*             If JOB = 'X', these arrays are not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK >= 2*N,          if JOB = 'X'; */
/*             LIWORK >= max(2*N,N*N), if JOB = 'A'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = 2, DWORK(1) contains the */
/*             optimal value of LDWORK. If JOB = 'A', then DWORK(2:N*N+1) */
/*             and DWORK(N*N+2:2*N*N+1) contain a real Schur form of the */
/*             closed-loop system matrix, Ac = A - G*X (if TRANA = 'N') */
/*             or Ac = A - X*G (if TRANA = 'T' or 'C'), and the */
/*             orthogonal matrix which reduced Ac to real Schur form, */
/*             respectively. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 4*N*N + 8*N + 1,               if JOB = 'X'; */
/*             LDWORK >= max( 4*N*N + 8*N, 6*N*N ) + 1, if JOB = 'A'. */
/*             For good performance, LDWORK should be larger, e.g., */
/*             LDWORK >= 4*N*N + 6*N +( 2*N+1 )*NB,     if JOB = 'X', */
/*             where NB is the optimal blocksize. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the Hamiltonian matrix has eigenvalues on the */
/*                   imaginary axis, so the solution and error bounds */
/*                   could not be computed; */
/*             = 2:  the iteration for the matrix sign function failed to */
/*                   converge after 50 iterations, but an approximate */
/*                   solution and error bounds (if JOB = 'A') have been */
/*                   computed; */
/*             = 3:  the system of linear equations for the solution is */
/*                   singular to working precision, so the solution and */
/*                   error bounds could not be computed; */
/*             = 4:  the matrix A-G*X (or A-X*G) cannot be reduced to */
/*                   Schur canonical form and condition number estimate */
/*                   and forward error estimate have not been computed. */

/*     METHOD */

/*     The Riccati equation is solved by the matrix sign function */
/*     approach [1], [2], implementing a scaling which enhances the */
/*     numerical stability [4]. */

/*     REFERENCES */

/*     [1] Bai, Z., Demmel, J., Dongarra, J., Petitet, A., Robinson, H., */
/*         and Stanley, K. */
/*         The spectral decomposition of nonsymmetric matrices on */
/*         distributed memory parallel computers. */
/*         SIAM J. Sci. Comput., vol. 18, pp. 1446-1461, 1997. */

/*     [2] Byers, R., He, C., and Mehrmann, V. */
/*         The matrix sign function method and the computation of */
/*         invariant subspaces. */
/*         SIAM J. Matrix Anal. Appl., vol. 18, pp. 615-632, 1997. */

/*     [3] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [4] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V., */
/*         DGRSVX and DMSRIC: Fortran 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Technical */
/*         University Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */

/*     The solution accuracy can be controlled by the output parameter */
/*     FERR. */

/*     FURTHER COMMENTS */

/*     The condition number of the Riccati equation is estimated as */

/*     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) + */
/*                 norm(Pi)*norm(G) ) / norm(X), */

/*     where Omega, Theta and Pi are linear operators defined by */

/*     Omega(W) = op(Ac)'*W + W*op(Ac), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))), */
/*        Pi(W) = inv(Omega(X*W*X)), */

/*     and the matrix Ac (the closed-loop system matrix) is given by */
/*        Ac = A - G*X, if TRANA = 'N', or */
/*        Ac = A - X*G, if TRANA = 'T' or 'C'. */

/*     The program estimates the quantities */

/*     sep(op(Ac),-op(Ac)') = 1 / norm(inv(Omega)), */

/*     norm(Theta) and norm(Pi) using 1-norm condition estimator. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [3]. */

/*     CONTRIBUTOR */

/*     P. Petkov, Tech. University of Sofia, March 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, continuous-time system, */
/*     optimal control, optimal regulator. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

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
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --wr;
    --wi;
    --iwork;
    --dwork;

    /* Function Body */
    all = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! all && ! lsame_(job, "X", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trana, 
	    "C", (ftnlen)1, (ftnlen)1) && ! notrna) {
	*info = -2;
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldg < max(1,*n)) {
	*info = -8;
    } else if (*ldq < max(1,*n)) {
	*info = -10;
    } else if (*ldx < max(1,*n)) {
	*info = -12;
    } else {

/*        Compute workspace. */

	if (all) {
/* Computing MAX */
	    i__1 = (*n << 2) * *n + (*n << 3) + 1, i__2 = *n * 6 * *n;
	    minwrk = max(i__1,i__2);
	} else {
	    minwrk = (*n << 2) * *n + (*n << 3) + 1;
	}
	if (*ldwork < minwrk) {
	    *info = -19;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB02PD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (all) {
	    *rcond = 1.;
	    *ferr = 0.;
	}
	dwork[1] = 1.;
	return 0;
    }

/*     Set tol. */

    eps = dlamch_("P", (ftnlen)1);
    tol = (doublereal) (*n) * 10. * eps;

/*     Compute the square-roots of the norms of the matrices Q and G . */

    qnorm2 = sqrt(dlansy_("1", uplo, n, &q[q_offset], ldq, &dwork[1], (ftnlen)
	    1, (ftnlen)1));
    gnorm2 = sqrt(dlansy_("1", uplo, n, &g[g_offset], ldg, &dwork[1], (ftnlen)
	    1, (ftnlen)1));

    n2 = *n << 1;

/*     Construct the lower (if UPLO = 'L') or upper (if UPLO = 'U') */
/*     triangle of the symmetric block-permuted Hamiltonian matrix. */
/*     During iteration, both the current iterate corresponding to the */
/*     Hamiltonian matrix, and its inverse are needed. To reduce the */
/*     workspace length, the transpose of the triangle specified by UPLO */
/*     of the current iterate H is saved in the opposite triangle, */
/*     suitably shifted with one column, and then the inverse of H */
/*     overwrites H. The triangles of the saved iterate and its inverse */
/*     are stored together in an 2*N-by-(2*N+1) matrix. For instance, if */
/*     UPLO = 'U', then the upper triangle is built starting from the */
/*     location 2*N+1 of the array DWORK, so that its transpose can be */
/*     stored in the lower triangle of DWORK. */
/*     Workspace: need   4*N*N,        if UPLO = 'L'; */
/*                       4*N*N + 2*N,  if UPLO = 'U'. */

    if (lower) {
	ini = 0;
	isv = n2;
	*(unsigned char *)loup = 'U';

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ij = (j - 1) * n2 + j;

	    i__2 = *n;
	    for (i__ = j; i__ <= i__2; ++i__) {
		dwork[ij] = -q[i__ + j * q_dim1];
		++ij;
/* L10: */
	    }

	    if (notrna) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ij] = -a[i__ + j * a_dim1];
		    ++ij;
/* L20: */
		}

	    } else {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ij] = -a[j + i__ * a_dim1];
		    ++ij;
/* L30: */
		}

	    }
/* L40: */
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ij = (*n + j - 1) * n2 + *n + j;

	    i__2 = *n;
	    for (i__ = j; i__ <= i__2; ++i__) {
		dwork[ij] = g[i__ + j * g_dim1];
		++ij;
/* L50: */
	    }

/* L60: */
	}

    } else {
	ini = n2;
	isv = 0;
	*(unsigned char *)loup = 'L';

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ij = j * n2 + 1;

	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[ij] = -q[i__ + j * q_dim1];
		++ij;
/* L70: */
	    }

/* L80: */
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ij = (*n + j) * n2 + 1;

	    if (notrna) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ij] = -a[j + i__ * a_dim1];
		    ++ij;
/* L90: */
		}

	    } else {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ij] = -a[i__ + j * a_dim1];
		    ++ij;
/* L100: */
		}

	    }

	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[ij] = g[i__ + j * g_dim1];
		++ij;
/* L110: */
	    }

/* L120: */
	}

    }

/*     Block-scaling. */

    iscl = 0;
    if (qnorm2 > gnorm2 && gnorm2 > 0.) {
	dlascl_(uplo, &c__0, &c__0, &qnorm2, &gnorm2, n, n, &dwork[ini + 1], &
		n2, &info2, (ftnlen)1);
	dlascl_(uplo, &c__0, &c__0, &gnorm2, &qnorm2, n, n, &dwork[n2 * *n + *
		n + ini + 1], &n2, &info2, (ftnlen)1);
	iscl = 1;
    }

/*     Workspace usage. */

    itau = n2 * n2;
    iwrk = itau + n2;

    lwamax = n2 * ilaenv_(&c__1, "DSYTRF", uplo, &n2, &c_n1, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);

/*     Compute the matrix sign function. */

    for (iter = 1; iter <= 50; ++iter) {

/*        Save the transpose of the corresponding triangle of the */
/*        current iterate in the free locations of the shifted opposite */
/*        triangle. */
/*        Workspace: need   4*N*N + 2*N. */

	if (lower) {

	    i__1 = n2;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(&i__, &dwork[i__], &n2, &dwork[i__ * n2 + 1], &c__1);
/* L130: */
	    }

	} else {

	    i__1 = n2;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(&i__, &dwork[i__ * n2 + 1], &c__1, &dwork[i__], &n2);
/* L140: */
	    }

	}

/*        Store the norm of the Hamiltonian matrix. */

	hnorm = dlansy_("F", uplo, &n2, &dwork[ini + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);

/*        Compute the inverse of the block-permuted Hamiltonian matrix. */
/*        Workspace: need   4*N*N + 2*N + 1; */
/*                   prefer 4*N*N + 2*N + 2*N*NB. */

	i__1 = *ldwork - iwrk;
	dsytrf_(uplo, &n2, &dwork[ini + 1], &n2, &iwork[1], &dwork[iwrk + 1], 
		&i__1, &info2, (ftnlen)1);
	if (info2 > 0) {
	    *info = 1;
	    return 0;
	}

/*        Workspace: need   4*N*N + 4*N. */

	dsytri_(uplo, &n2, &dwork[ini + 1], &n2, &iwork[1], &dwork[iwrk + 1], 
		&info2, (ftnlen)1);

/*        Block-permutation of the inverse matrix. */

	if (lower) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		ij2 = (*n + j - 1) * n2 + *n + j;

		i__2 = (j - 1) * n2 + *n;
		for (ij1 = (j - 1) * n2 + j; ij1 <= i__2; ++ij1) {
		    temp = dwork[ij1];
		    dwork[ij1] = -dwork[ij2];
		    dwork[ij2] = -temp;
		    ++ij2;
/* L150: */
		}

		i__2 = j - 1;
		dswap_(&i__2, &dwork[*n + j], &n2, &dwork[(j - 1) * n2 + *n + 
			1], &c__1);
/* L160: */
	    }

	} else {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		ij2 = (*n + j) * n2 + *n + 1;

		i__2 = j * n2 + j;
		for (ij1 = j * n2 + 1; ij1 <= i__2; ++ij1) {
		    temp = dwork[ij1];
		    dwork[ij1] = -dwork[ij2];
		    dwork[ij2] = -temp;
		    ++ij2;
/* L170: */
		}

		i__2 = j - 1;
		dswap_(&i__2, &dwork[(*n + 1) * n2 + j], &n2, &dwork[(*n + j) 
			* n2 + 1], &c__1);
/* L180: */
	    }

	}

/*        Scale the Hamiltonian matrix and its inverse and compute */
/*        the next iterate. */

	hinnrm = dlansy_("F", uplo, &n2, &dwork[ini + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);
	scale = sqrt(hinnrm / hnorm);

	if (lower) {

	    i__1 = n2;
	    for (j = 1; j <= i__1; ++j) {
		ji = (j - 1) * n2 + j;

		i__2 = j * n2;
		for (ij = ji; ij <= i__2; ++ij) {
		    ji += n2;
		    dwork[ij] = (dwork[ij] / scale + dwork[ji] * scale) / 2.;
		    dwork[ji] -= dwork[ij];
/* L190: */
		}

/* L200: */
	    }

	} else {

	    i__1 = n2;
	    for (j = 1; j <= i__1; ++j) {
		ji = j;

		i__2 = j * n2 + j;
		for (ij = j * n2 + 1; ij <= i__2; ++ij) {
		    dwork[ij] = (dwork[ij] / scale + dwork[ji] * scale) / 2.;
		    dwork[ji] -= dwork[ij];
		    ji += n2;
/* L210: */
		}

/* L220: */
	    }

	}

/*        Test for convergence. */

	conv = dlansy_("F", loup, &n2, &dwork[isv + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);
	if (conv <= tol * hnorm) {
	    goto L240;
	}
/* L230: */
    }

/*     No convergence after MAXIT iterations, but an approximate solution */
/*     has been found. */

    *info = 2;

L240:

/*     If UPLO = 'U', shift the upper triangle one column to the left. */

    if (! lower) {
	dlacpy_("U", &n2, &n2, &dwork[ini + 1], &n2, &dwork[1], &n2, (ftnlen)
		1);
    }

/*     Divide the triangle elements by -2 and then fill-in the other */
/*     triangle by symmetry. */

    if (lower) {

	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = n2 - i__ + 1;
	    dscal_(&i__2, &c_b57, &dwork[(i__ - 1) * n2 + i__], &c__1);
/* L250: */
	}

    } else {

	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dscal_(&i__, &c_b57, &dwork[(i__ - 1) * n2 + 1], &c__1);
/* L260: */
	}

    }
    ma02ed_(uplo, &n2, &dwork[1], &n2, (ftnlen)1);

/*     Back block-permutation. */

    i__1 = n2;
    for (j = 1; j <= i__1; ++j) {

	i__2 = (j - 1) * n2 + *n;
	for (i__ = (j - 1) * n2 + 1; i__ <= i__2; ++i__) {
	    temp = dwork[i__];
	    dwork[i__] = -dwork[i__ + *n];
	    dwork[i__ + *n] = temp;
/* L270: */
	}

/* L280: */
    }

/*     Compute the QR decomposition of the projector onto the stable */
/*     invariant subspace. */
/*     Workspace: need   4*N*N + 8*N + 1. */
/*                prefer 4*N*N + 6*N + ( 2*N+1 )*NB. */

    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
	dwork[(i__ - 1) * n2 + i__] += .5;
/* L290: */
    }

    i__1 = *ldwork - iwrk;
    dgeqp3_(&n2, &n2, &dwork[1], &n2, &iwork[1], &dwork[itau + 1], &dwork[
	    iwrk + 1], &i__1, &info2);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk + 1];
    lwamax = max(i__1,lwamax);

/*     Accumulate the orthogonal transformations. Note that only the */
/*     first N columns of the array DWORK, returned by DGEQP3, are */
/*     needed, so that the last N columns of DWORK are used to get the */
/*     orthogonal basis for the stable invariant subspace. */
/*     Workspace: need   4*N*N + 3*N. */
/*                prefer 4*N*N + 2*N + N*NB. */

    ib = *n * *n;
    iaf = n2 * *n;
    dlaset_("F", &n2, n, &c_b66, &c_b67, &dwork[iaf + 1], &n2, (ftnlen)1);
    i__1 = *ldwork - iwrk;
    dormqr_("L", "N", &n2, n, n, &dwork[1], &n2, &dwork[itau + 1], &dwork[iaf 
	    + 1], &n2, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk + 1];
    lwamax = iwrk + max(i__1,lwamax);

/*     Store the matrices V11 and V21' . */

    dlacpy_("F", n, n, &dwork[iaf + 1], &n2, &dwork[1], n, (ftnlen)1);
    ma02ad_("F", n, n, &dwork[iaf + *n + 1], &n2, &dwork[ib + 1], n, (ftnlen)
	    1);

    ir = iaf + ib;
    ic = ir + *n;
    ifr = ic + *n;
    ibr = ifr + *n;
    iwrk = ibr + *n;

/*     Compute the solution matrix X . */
/*     Workspace: need   3*N*N + 8*N. */

    dgesvx_("E", "T", n, n, &dwork[1], n, &dwork[iaf + 1], n, &iwork[1], 
	    equed, &dwork[ir + 1], &dwork[ic + 1], &dwork[ib + 1], n, &x[
	    x_offset], ldx, rcond, &dwork[ifr + 1], &dwork[ibr + 1], &dwork[
	    iwrk + 1], &iwork[*n + 1], &info2, (ftnlen)1, (ftnlen)1, (ftnlen)
	    1);
    if (info2 > 0) {
	*info = 3;
	return 0;
    }

/*     Symmetrize the solution. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    temp = (x[i__ + j * x_dim1] + x[j + i__ * x_dim1]) / 2.;
	    x[i__ + j * x_dim1] = temp;
	    x[j + i__ * x_dim1] = temp;
/* L300: */
	}

/* L310: */
    }

/*     Undo scaling for the solution matrix. */

    if (iscl == 1) {
	dlascl_("G", &c__0, &c__0, &gnorm2, &qnorm2, n, n, &x[x_offset], ldx, 
		&info2, (ftnlen)1);
    }

    if (all) {

/*        Compute the estimates of the reciprocal condition number and */
/*        error bound. */
/*        Workspace usage. */

	it = 1;
	iu = it + *n * *n;
	iwrk = iu + *n * *n;

	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[it + 1], n, (ftnlen)4)
		;
	if (notrna) {

/*           Compute Ac = A-G*X . */

	    dsymm_("L", uplo, n, n, &c_b81, &g[g_offset], ldg, &x[x_offset], 
		    ldx, &c_b67, &dwork[it + 1], n, (ftnlen)1, (ftnlen)1);
	} else {

/*           Compute Ac = A-X*G . */

	    dsymm_("R", uplo, n, n, &c_b81, &g[g_offset], ldg, &x[x_offset], 
		    ldx, &c_b67, &dwork[it + 1], n, (ftnlen)1, (ftnlen)1);
	}

/*        Compute the Schur factorization of Ac . */
/*        Workspace: need   2*N*N + 5*N + 1; */
/*                   prefer larger. */

	i__1 = *ldwork - iwrk;
	dgees_("V", "N", (L_fp)select_, n, &dwork[it + 1], n, &sdim, &wr[1], &
		wi[1], &dwork[iu + 1], n, &dwork[iwrk + 1], &i__1, bwork, &
		info2, (ftnlen)1, (ftnlen)1);
	if (info2 > 0) {
	    *info = 4;
	    return 0;
	}
/* Computing MAX */
	i__1 = (integer) dwork[iwrk + 1];
	lwamax = iwrk + max(i__1,lwamax);

/*        Estimate the reciprocal condition number and the forward error. */
/*        Workspace: need   6*N*N + 1; */
/*                   prefer larger. */

	i__1 = *ldwork - iwrk;
	sb02qd_("B", "F", trana, uplo, "O", n, &a[a_offset], lda, &dwork[it + 
		1], n, &dwork[iu + 1], n, &g[g_offset], ldg, &q[q_offset], 
		ldq, &x[x_offset], ldx, &sep, rcond, ferr, &iwork[1], &dwork[
		iwrk + 1], &i__1, &info2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk + 1];
	lwamax = iwrk + max(i__1,lwamax);
    }

    dwork[1] = (doublereal) lwamax;
    return 0;
/* *** Last line of SB02PD */
} /* sb02pd_ */

