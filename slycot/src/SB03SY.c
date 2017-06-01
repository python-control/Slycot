/* SB03SY.f -- translated by f2c (version 20100827).
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
static doublereal c_b22 = 1.;
static doublereal c_b23 = .5;

/* Subroutine */ int sb03sy_(char *job, char *trana, char *lyapun, integer *n,
	 doublereal *t, integer *ldt, doublereal *u, integer *ldu, doublereal 
	*xa, integer *ldxa, doublereal *sepd, doublereal *thnorm, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	job_len, ftnlen trana_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, u_dim1, u_offset, xa_dim1, xa_offset, i__1, 
	    i__2;

    /* Local variables */
    static integer nn;
    static doublereal est;
    static integer kase, itmp;
    static char uplo[1];
    static integer info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sb03mx_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen);
    static logical wants, wantt;
    extern /* Subroutine */ int dsyr2k_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical update;
    static char tranat[1];
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical notrna;


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

/*     To estimate the "separation" between the matrices op(A) and */
/*     op(A)', */

/*     sepd(op(A),op(A)') = min norm(op(A)'*X*op(A) - X)/norm(X) */
/*                        = 1 / norm(inv(Omega)) */

/*     and/or the 1-norm of Theta, where op(A) = A or A' (A**T), and */
/*     Omega and Theta are linear operators associated to the real */
/*     discrete-time Lyapunov matrix equation */

/*            op(A)'*X*op(A) - X = C, */

/*     defined by */

/*     Omega(W) = op(A)'*W*op(A) - W, */
/*     Theta(W) = inv(Omega(op(W)'*X*op(A) + op(A)'*X*op(W))). */

/*     The 1-norm condition estimators are used. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'S':  Compute the separation only; */
/*             = 'T':  Compute the norm of Theta only; */
/*             = 'B':  Compute both the separation and the norm of Theta. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved, as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and X.  N >= 0. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of A. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array must contain the */
/*             orthogonal matrix U from a real Schur factorization of A. */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     XA      (input) DOUBLE PRECISION array, dimension (LDXA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix product X*op(A), if LYAPUN = 'O', or U'*X*U*op(T), */
/*             if LYAPUN = 'R', in the Lyapunov equation. */
/*             If JOB = 'S', the array XA is not referenced. */

/*     LDXA    INTEGER */
/*             The leading dimension of array XA. */
/*             LDXA >= 1,        if JOB = 'S'; */
/*             LDXA >= MAX(1,N), if JOB = 'T' or 'B'. */

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO >= 0, SEPD contains */
/*             the estimated quantity sepd(op(A),op(A)'). */
/*             If JOB = 'T' or N = 0, SEPD is not referenced. */

/*     THNORM  (output) DOUBLE PRECISION */
/*             If JOB = 'T' or JOB = 'B', and INFO >= 0, THNORM contains */
/*             the estimated 1-norm of operator Theta. */
/*             If JOB = 'S' or N = 0, THNORM is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 0,            if N = 0; */
/*             LDWORK >= MAX(3,2*N*N), if N > 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = N+1:  if T has (almost) reciprocal eigenvalues; */
/*                   perturbed values were used to solve Lyapunov */
/*                   equations (but the matrix T is unchanged). */

/*     METHOD */

/*     SEPD is defined as */

/*            sepd( op(A), op(A)' ) = sigma_min( K ) */

/*     where sigma_min(K) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        K = kprod( op(A)', op(A)' ) - I(N**2). */

/*     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the */
/*     Kronecker product. The routine estimates sigma_min(K) by the */
/*     reciprocal of an estimate of the 1-norm of inverse(K), computed as */
/*     suggested in [1]. This involves the solution of several discrete- */
/*     time Lyapunov equations, either direct or transposed. The true */
/*     reciprocal 1-norm of inverse(K) cannot differ from sigma_min(K) by */
/*     more than a factor of N. */
/*     The 1-norm of Theta is estimated similarly. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or */
/*         complex matrix, with applications to condition estimation. */
/*         ACM Trans. Math. Softw., 14, pp. 381-396, 1988. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     When SEPD is zero, the routine returns immediately, with THNORM */
/*     (if requested) not set. In this case, the equation is singular. */
/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. Partly based on DDLSVX (and then SB03SD) by P. Petkov, */
/*     Tech. University of Sofia, March 1998 (and December 1998). */

/*     REVISIONS */

/*     February 6, 1999, V. Sima, Katholieke Univ. Leuven, Belgium. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    xa_dim1 = *ldxa;
    xa_offset = 1 + xa_dim1;
    xa -= xa_offset;
    --iwork;
    --dwork;

    /* Function Body */
    wants = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
    wantt = lsame_(job, "T", (ftnlen)1, (ftnlen)1);
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

    nn = *n * *n;
    *info = 0;
    if (! (wants || wantt || lsame_(job, "B", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldt < max(1,*n)) {
	*info = -6;
    } else if (*ldu < 1 || update && *ldu < *n) {
	*info = -8;
    } else if (*ldxa < 1 || ! wants && *ldxa < *n) {
	*info = -10;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 3, i__2 = nn << 1;
	if (*ldwork < 0 || *ldwork < max(i__1,i__2) && *n > 0) {
	    *info = -15;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB03SY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    itmp = nn + 1;

    if (notrna) {
	*(unsigned char *)tranat = 'T';
    } else {
	*(unsigned char *)tranat = 'N';
    }

    if (! wantt) {

/*        Estimate sepd(op(A),op(A)'). */
/*        Workspace:  max(3,2*N*N). */

	kase = 0;

/*        REPEAT */
L10:
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
		*(unsigned char *)uplo = 'U';
	    } else {
		*(unsigned char *)uplo = 'L';
	    }

	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
		i__1 = *n + 1;
		dscal_(n, &c_b23, &dwork[1], &i__1);
	    }
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

	    if (kase == 1) {

/*              Solve op(T)'*Y*op(T) - Y = scale*RHS. */

		sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
	    } else {

/*              Solve op(T)*W*op(T)' - W = scale*RHS. */

		sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
	    }

	    if (info2 > 0) {
		*info = *n + 1;
	    }

	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
		i__1 = *n + 1;
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
	    }

	    goto L10;
	}
/*        UNTIL KASE = 0 */

	if (est > scale) {
	    *sepd = scale / est;
	} else {
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
	    if (scale < est * bignum) {
		*sepd = scale / est;
	    } else {
		*sepd = bignum;
	    }
	}

/*        Return if the equation is singular. */

	if (*sepd == 0.) {
	    return 0;
	}
    }

    if (! wants) {

/*        Estimate norm(Theta). */
/*        Workspace:  max(3,2*N*N). */

	kase = 0;

/*        REPEAT */
L20:
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
		*(unsigned char *)uplo = 'U';
	    } else {
		*(unsigned char *)uplo = 'L';
	    }

/*           Fill in the remaining triangle of the symmetric matrix. */

	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

/*           Compute RHS = op(W)'*X*op(A) + op(A)'*X*op(W). */

	    dsyr2k_(uplo, tranat, n, n, &c_b22, &dwork[1], n, &xa[xa_offset], 
		    ldxa, &c_b21, &dwork[itmp], n, (ftnlen)1, (ftnlen)1);
	    dlacpy_(uplo, n, n, &dwork[itmp], n, &dwork[1], n, (ftnlen)1);

	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
		i__1 = *n + 1;
		dscal_(n, &c_b23, &dwork[1], &i__1);
	    }
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

	    if (kase == 1) {

/*              Solve op(T)'*Y*op(T) - Y = scale*RHS. */

		sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
	    } else {

/*              Solve op(T)*W*op(T)' - W = scale*RHS. */

		sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
	    }

	    if (info2 > 0) {
		*info = *n + 1;
	    }

	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
		i__1 = *n + 1;
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
	    }

	    goto L20;
	}
/*        UNTIL KASE = 0 */

	if (est < scale) {
	    *thnorm = est / scale;
	} else {
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
	    if (est < scale * bignum) {
		*thnorm = est / scale;
	    } else {
		*thnorm = bignum;
	    }
	}
    }

    return 0;
/* *** Last line of SB03SY *** */
} /* sb03sy_ */

