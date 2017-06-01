/* SB03QD.f -- translated by f2c (version 20100827).
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

static doublereal c_b21 = 0.;
static doublereal c_b23 = -1.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b41 = 1.;

/* Subroutine */ int sb03qd_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *scale, doublereal *a, integer *
	lda, doublereal *t, integer *ldt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *x, integer *ldx, 
	doublereal *sep, doublereal *rcond, doublereal *ferr, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen fact_len, ftnlen trana_len, ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, jj, nn, jx, ldw;
    static doublereal eps;
    static logical jobb, jobc, jobe;
    static integer iabs, sdim;
    static char sjob[1];
    static integer ires, ixbs;
    static doublereal epsn, temp, tmax;
    static integer iwrk;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), mb01ud_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    extern /* Subroutine */ int mb01uw_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal anorm, cnorm;
    extern /* Subroutine */ int sb03qx_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), sb03qy_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical bwork[1];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lower;
    static doublereal xnorm;
    extern /* Subroutine */ int dsyr2k_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen), 
	    dlanhs_(char *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static logical update;
    static char tranat[1];
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical notrna;
    static doublereal xanorm, thnorm;
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

/*     To estimate the conditioning and compute an error bound on the */
/*     solution of the real continuous-time Lyapunov matrix equation */

/*         op(A)'*X + X*op(A) = scale*C */

/*     where op(A) = A or A' (A**T) and C is symmetric (C = C**T). The */
/*     matrix A is N-by-N, the right hand side C and the solution X are */
/*     N-by-N symmetric matrices, and scale is a given scale factor. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'C':  Compute the reciprocal condition number only; */
/*             = 'E':  Compute the error bound only; */
/*             = 'B':  Compute both the reciprocal condition number and */
/*                     the error bound. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix A is supplied on entry, as follows: */
/*             = 'F':  On entry, T and U (if LYAPUN = 'O') contain the */
/*                     factors from the real Schur factorization of the */
/*                     matrix A; */
/*             = 'N':  The Schur factorization of A will be computed */
/*                     and the factors will be stored in T and U (if */
/*                     LYAPUN = 'O'). */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the symmetric matrix C is to be */
/*             used, as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved in the iterative estimation process, */
/*             as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X and C.  N >= 0. */

/*     SCALE   (input) DOUBLE PRECISION */
/*             The scale factor, scale, set by a Lyapunov solver. */
/*             0 <= SCALE <= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If FACT = 'N' or LYAPUN = 'O', the leading N-by-N part of */
/*             this array must contain the original matrix A. */
/*             If FACT = 'F' and LYAPUN = 'R', A is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if FACT = 'N' or  LYAPUN = 'O'; */
/*             LDA >= 1,        if FACT = 'F' and LYAPUN = 'R'. */

/*     T       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDT,N) */
/*             If FACT = 'F', then on entry the leading N-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             quasi-triangular matrix T in Schur canonical form from a */
/*             Schur factorization of A. */
/*             If FACT = 'N', then this array need not be set on input. */
/*             On exit, (if INFO = 0 or INFO = N+1, for FACT = 'N') the */
/*             leading N-by-N upper Hessenberg part of this array */
/*             contains the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of A. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,N). */

/*     U       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDU,N) */
/*             If LYAPUN = 'O' and FACT = 'F', then U is an input */
/*             argument and on entry, the leading N-by-N part of this */
/*             array must contain the orthogonal matrix U from a real */
/*             Schur factorization of A. */
/*             If LYAPUN = 'O' and FACT = 'N', then U is an output */
/*             argument and on exit, if INFO = 0 or INFO = N+1, it */
/*             contains the orthogonal N-by-N matrix from a real Schur */
/*             factorization of A. */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix C of the original Lyapunov equation (with */
/*             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov */
/*             equation (with matrix T), if LYAPUN = 'R'. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix C of the original Lyapunov equation (with */
/*             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov */
/*             equation (with matrix T), if LYAPUN = 'R'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,N). */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             symmetric solution matrix X of the original Lyapunov */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             reduced Lyapunov equation (with matrix T), if */
/*             LYAPUN = 'R'. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', the estimated quantity */
/*             sep(op(A),-op(A)'). */
/*             If N = 0, or X = 0, or JOB = 'E', SEP is not referenced. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', an estimate of the reciprocal */
/*             condition number of the continuous-time Lyapunov equation. */
/*             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively. */
/*             If JOB = 'E', RCOND is not referenced. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'E' or JOB = 'B', an estimated forward error */
/*             bound for the solution X. If XTRUE is the true solution, */
/*             FERR bounds the magnitude of the largest entry in */
/*             (X - XTRUE) divided by the magnitude of the largest entry */
/*             in X. */
/*             If N = 0 or X = 0, FERR is set to 0. */
/*             If JOB = 'C', FERR is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             If JOB = 'C', then */
/*             LDWORK >= MAX(1,2*N*N),         if FACT = 'F'; */
/*             LDWORK >= MAX(1,2*N*N,5*N),     if FACT = 'N'. */
/*             If JOB = 'E', or JOB = 'B', and LYAPUN  = 'O', then */
/*             LDWORK >= MAX(1,3*N*N),         if FACT = 'F'; */
/*             LDWORK >= MAX(1,3*N*N,5*N),     if FACT = 'N'. */
/*             If JOB = 'E', or JOB = 'B', and LYAPUN  = 'R', then */
/*             LDWORK >= MAX(1,3*N*N+N-1),     if FACT = 'F'; */
/*             LDWORK >= MAX(1,3*N*N+N-1,5*N), if FACT = 'N'. */
/*             For optimum performance LDWORK should sometimes be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, i <= N, the QR algorithm failed to */
/*                   complete the reduction to Schur canonical form (see */
/*                   LAPACK Library routine DGEES); on exit, the matrix */
/*                   T(i+1:N,i+1:N) contains the partially converged */
/*                   Schur form, and DWORK(i+1:N) and DWORK(N+i+1:2*N) */
/*                   contain the real and imaginary parts, respectively, */
/*                   of the converged eigenvalues; this error is unlikely */
/*                   to appear; */
/*             = N+1:  if the matrices T and -T' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve Lyapunov equations, but the matrix T, if given */
/*                   (for FACT = 'F'), is unchanged. */

/*     METHOD */

/*     The condition number of the continuous-time Lyapunov equation is */
/*     estimated as */

/*     cond = (norm(Theta)*norm(A) + norm(inv(Omega))*norm(C))/norm(X), */

/*     where Omega and Theta are linear operators defined by */

/*     Omega(W) = op(A)'*W + W*op(A), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))). */

/*     The routine estimates the quantities */

/*     sep(op(A),-op(A)') = 1 / norm(inv(Omega)) */

/*     and norm(Theta) using 1-norm condition estimators. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [1]. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */
/*     The accuracy of the estimates obtained depends on the solution */
/*     accuracy and on the properties of the 1-norm estimator. */

/*     FURTHER COMMENTS */

/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */
/*     When SEP is computed and it is zero, the routine returns */
/*     immediately, with RCOND and FERR (if requested) set to 0 and 1, */
/*     respectively. In this case, the equation is singular. */

/*     CONTRIBUTORS */

/*     P. Petkov, Tech. University of Sofia, December 1998. */
/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 2003. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

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
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --iwork;
    --dwork;

    /* Function Body */
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
    jobb = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

    nn = *n * *n;
    if (jobc) {
	ldw = nn << 1;
    } else {
	ldw = nn * 3;
    }
    if (! (jobc || update)) {
	ldw = ldw + *n - 1;
    }

    *info = 0;
    if (! (jobb || jobc || jobe)) {
	*info = -1;
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*scale < 0. || *scale > 1.) {
	*info = -7;
    } else if (*lda < 1 || *lda < *n && (update || nofact)) {
	*info = -9;
    } else if (*ldt < max(1,*n)) {
	*info = -11;
    } else if (*ldu < 1 || *ldu < *n && update) {
	*info = -13;
    } else if (*ldc < max(1,*n)) {
	*info = -15;
    } else if (*ldx < max(1,*n)) {
	*info = -17;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = ldw, i__2 = *n * 5;
	if (*ldwork < 1 || *ldwork < ldw && ! nofact || *ldwork < max(i__1,
		i__2) && nofact) {
	    *info = -23;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB03QD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (! jobe) {
	    *rcond = 1.;
	}
	if (! jobc) {
	    *ferr = 0.;
	}
	dwork[1] = 1.;
	return 0;
    }

/*     Compute the 1-norm of the matrix X. */

    xnorm = dlansy_("1-norm", uplo, n, &x[x_offset], ldx, &dwork[1], (ftnlen)
	    6, (ftnlen)1);
    if (xnorm == 0.) {

/*        The solution is zero. */

	if (! jobe) {
	    *rcond = 0.;
	}
	if (! jobc) {
	    *ferr = 0.;
	}
	dwork[1] = (doublereal) (*n);
	return 0;
    }

/*     Compute the 1-norm of A or T. */

    if (nofact || update) {
	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);
    } else {
	anorm = dlanhs_("1-norm", n, &t[t_offset], ldt, &dwork[1], (ftnlen)6);
    }

/*     For the special case A = 0, set SEP and RCOND to 0. */
/*     For the special case A = I, set SEP to 2 and RCOND to 1. */
/*     A quick test is used in general. */

    if (anorm == 1.) {
	if (nofact || update) {
	    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
	} else {
	    dlacpy_("Full", n, n, &t[t_offset], ldt, &dwork[1], n, (ftnlen)4);
	    if (*n > 2) {
		i__1 = *n - 2;
		i__2 = *n - 2;
		dlaset_("Lower", &i__1, &i__2, &c_b21, &c_b21, &dwork[3], n, (
			ftnlen)5);
	    }
	}
	dwork[nn + 1] = 1.;
	i__1 = *n + 1;
	daxpy_(n, &c_b23, &dwork[nn + 1], &c__0, &dwork[1], &i__1);
	if (dlange_("Max", n, n, &dwork[1], n, &dwork[1], (ftnlen)3) == 0.) {
	    if (! jobe) {
		*sep = 2.;
		*rcond = 1.;
	    }
	    if (jobc) {
		dwork[1] = (doublereal) (nn + 1);
		return 0;
	    } else {

/*              Set FERR for the special case A = I. */

		dlacpy_(uplo, n, n, &x[x_offset], ldx, &dwork[1], n, (ftnlen)
			1);

		if (lower) {
		    i__1 = *n;
		    for (j = 1; j <= i__1; ++j) {
			i__2 = *n - j + 1;
			d__1 = -(*scale) / 2.;
			daxpy_(&i__2, &d__1, &c__[j + j * c_dim1], &c__1, &
				dwork[(j - 1) * *n + j], &c__1);
/* L10: */
		    }
		} else {
		    i__1 = *n;
		    for (j = 1; j <= i__1; ++j) {
			d__1 = -(*scale) / 2.;
			daxpy_(&j, &d__1, &c__[j * c_dim1 + 1], &c__1, &dwork[
				(j - 1) * *n + 1], &c__1);
/* L20: */
		    }
		}

/* Computing MIN */
		d__1 = 1., d__2 = dlansy_("1-norm", uplo, n, &dwork[1], n, &
			dwork[nn + 1], (ftnlen)6, (ftnlen)1) / xnorm;
		*ferr = min(d__1,d__2);
		dwork[1] = (doublereal) (nn + *n);
		return 0;
	    }
	}

    } else if (anorm == 0.) {
	if (! jobe) {
	    *sep = 0.;
	    *rcond = 0.;
	}
	if (! jobc) {
	    *ferr = 1.;
	}
	dwork[1] = (doublereal) (*n);
	return 0;
    }

/*     General case. */

    cnorm = dlansy_("1-norm", uplo, n, &c__[c_offset], ldc, &dwork[1], (
	    ftnlen)6, (ftnlen)1);

/*     Workspace usage. */

    iabs = 0;
    ixbs = iabs + nn;
    ires = ixbs + nn;
    iwrk = ires + nn;
    wrkopt = 0;

    if (nofact) {

/*        Compute the Schur factorization of A, A = U*T*U'. */
/*        Workspace:  need   5*N; */
/*                    prefer larger. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance.) */

	dlacpy_("Full", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)4)
		;
	if (update) {
	    *(unsigned char *)sjob = 'V';
	} else {
	    *(unsigned char *)sjob = 'N';
	}
	i__1 = *ldwork - (*n << 1);
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &dwork[1], &dwork[*n + 1], &u[u_offset], ldu, &dwork[(*
		n << 1) + 1], &i__1, bwork, info, (ftnlen)1, (ftnlen)11);
	if (*info > 0) {
	    return 0;
	}
	wrkopt = (integer) dwork[(*n << 1) + 1] + (*n << 1);
    }

    if (! jobe) {

/*        Estimate sep(op(A),-op(A)') = sep(op(T),-op(T)') and */
/*        norm(Theta). */
/*        Workspace 2*N*N. */

	sb03qy_("Both", trana, lyapun, n, &t[t_offset], ldt, &u[u_offset], 
		ldu, &x[x_offset], ldx, sep, &thnorm, &iwork[1], &dwork[1], 
		ldwork, info, (ftnlen)4, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
	i__1 = wrkopt, i__2 = nn << 1;
	wrkopt = max(i__1,i__2);

/*        Return if the equation is singular. */

	if (*sep == 0.) {
	    *rcond = 0.;
	    if (jobb) {
		*ferr = 1.;
	    }
	    dwork[1] = (doublereal) wrkopt;
	    return 0;
	}

/*        Estimate the reciprocal condition number. */

/* Computing MAX */
	d__1 = max(*sep,xnorm);
	tmax = max(d__1,anorm);
	if (tmax <= 1.) {
	    temp = *sep * xnorm;
	    denom = *scale * cnorm + *sep * anorm * thnorm;
	} else {
	    temp = *sep / tmax * (xnorm / tmax);
	    denom = *scale / tmax * (cnorm / tmax) + *sep / tmax * (anorm / 
		    tmax) * thnorm;
	}
	if (temp >= denom) {
	    *rcond = 1.;
	} else {
	    *rcond = temp / denom;
	}
    }

    if (! jobc) {

/*        Form a triangle of the residual matrix */
/*        R = op(A)'*X + X*op(A) - scale*C, or */
/*        R = op(T)'*X + X*op(T) - scale*C, */
/*        exploiting the symmetry. */
/*        Workspace 3*N*N. */

	if (notrna) {
	    *(unsigned char *)tranat = 'T';
	} else {
	    *(unsigned char *)tranat = 'N';
	}

	if (update) {

	    dlacpy_(uplo, n, n, &c__[c_offset], ldc, &dwork[ires + 1], n, (
		    ftnlen)1);
	    d__1 = -(*scale);
	    dsyr2k_(uplo, tranat, n, n, &c_b41, &a[a_offset], lda, &x[
		    x_offset], ldx, &d__1, &dwork[ires + 1], n, (ftnlen)1, (
		    ftnlen)1);
	} else {
	    mb01ud_("Right", trana, n, n, &c_b41, &t[t_offset], ldt, &x[
		    x_offset], ldx, &dwork[ires + 1], n, info, (ftnlen)5, (
		    ftnlen)1);
	    jj = ires + 1;
	    if (lower) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - j + 1;
		    daxpy_(&i__2, &c_b41, &dwork[jj], n, &dwork[jj], &c__1);
		    i__2 = *n - j + 1;
		    d__1 = -(*scale);
		    daxpy_(&i__2, &d__1, &c__[j + j * c_dim1], &c__1, &dwork[
			    jj], &c__1);
		    jj = jj + *n + 1;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    daxpy_(&j, &c_b41, &dwork[ires + j], n, &dwork[jj], &c__1)
			    ;
		    d__1 = -(*scale);
		    daxpy_(&j, &d__1, &c__[j * c_dim1 + 1], &c__1, &dwork[jj],
			     &c__1);
		    jj += *n;
/* L40: */
		}
	    }
	}

/* Computing MAX */
	i__1 = wrkopt, i__2 = nn * 3;
	wrkopt = max(i__1,i__2);

/*        Get the machine precision. */

	eps = dlamch_("Epsilon", (ftnlen)7);
	epsn = eps * (doublereal) (*n + 3);
	temp = eps * 3. * *scale;

/*        Add to abs(R) a term that takes account of rounding errors in */
/*        forming R: */
/*          abs(R) := abs(R) + EPS*(3*scale*abs(C) + */
/*                    (n+3)*(abs(op(A))'*abs(X) + abs(X)*abs(op(A)))), or */
/*          abs(R) := abs(R) + EPS*(3*scale*abs(C) + */
/*                    (n+3)*(abs(op(T))'*abs(X) + abs(X)*abs(op(T)))), */
/*        where EPS is the machine precision. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dwork[ixbs + (j - 1) * *n + i__] = (d__1 = x[i__ + j * x_dim1]
			, abs(d__1));
/* L50: */
	    }
/* L60: */
	}

	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    dwork[ires + (j - 1) * *n + i__] = temp * (d__1 = c__[i__ 
			    + j * c_dim1], abs(d__1)) + (d__2 = dwork[ires + (
			    j - 1) * *n + i__], abs(d__2));
/* L70: */
		}
/* L80: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ires + (j - 1) * *n + i__] = temp * (d__1 = c__[i__ 
			    + j * c_dim1], abs(d__1)) + (d__2 = dwork[ires + (
			    j - 1) * *n + i__], abs(d__2));
/* L90: */
		}
/* L100: */
	    }
	}

	if (update) {

/*           Workspace 3*N*N. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = a[i__ + j * 
			    a_dim1], abs(d__1));
/* L110: */
		}
/* L120: */
	    }

	    dsyr2k_(uplo, tranat, n, n, &epsn, &dwork[iabs + 1], n, &dwork[
		    ixbs + 1], n, &c_b41, &dwork[ires + 1], n, (ftnlen)1, (
		    ftnlen)1);
	} else {

/*           Workspace 3*N*N + N - 1. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__3 = j + 1;
		i__2 = min(i__3,*n);
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = t[i__ + j * 
			    t_dim1], abs(d__1));
/* L130: */
		}
/* L140: */
	    }

	    i__1 = *ldwork - iwrk;
	    mb01uw_("Left", tranat, n, n, &epsn, &dwork[iabs + 1], n, &dwork[
		    ixbs + 1], n, &dwork[iwrk + 1], &i__1, info, (ftnlen)4, (
		    ftnlen)1);
	    jj = ires + 1;
	    jx = ixbs + 1;
	    if (lower) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - j + 1;
		    daxpy_(&i__2, &c_b41, &dwork[jx], n, &dwork[jx], &c__1);
		    i__2 = *n - j + 1;
		    daxpy_(&i__2, &c_b41, &dwork[jx], &c__1, &dwork[jj], &
			    c__1);
		    jj = jj + *n + 1;
		    jx = jx + *n + 1;
/* L150: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    daxpy_(&j, &c_b41, &dwork[ixbs + j], n, &dwork[jx], &c__1)
			    ;
		    daxpy_(&j, &c_b41, &dwork[jx], &c__1, &dwork[jj], &c__1);
		    jj += *n;
		    jx += *n;
/* L160: */
		}
	    }

/* Computing MAX */
	    i__1 = wrkopt, i__2 = nn * 3 + *n - 1;
	    wrkopt = max(i__1,i__2);
	}

/*        Compute forward error bound, using matrix norm estimator. */
/*        Workspace 3*N*N. */

	xanorm = dlansy_("Max", uplo, n, &x[x_offset], ldx, &dwork[1], (
		ftnlen)3, (ftnlen)1);

	sb03qx_(trana, uplo, lyapun, n, &xanorm, &t[t_offset], ldt, &u[
		u_offset], ldu, &dwork[ires + 1], n, ferr, &iwork[1], &dwork[
		1], &ires, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    }

    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of SB03QD *** */
} /* sb03qd_ */

