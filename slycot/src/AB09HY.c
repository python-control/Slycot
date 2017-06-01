/* AB09HY.f -- translated by f2c (version 20100827).
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
static logical c_true = TRUE_;
static doublereal c_b15 = 1.;
static doublereal c_b34 = -1.;
static doublereal c_b39 = 0.;

/* Subroutine */ int ab09hy_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *scalec, 
	doublereal *scaleo, doublereal *s, integer *lds, doublereal *r__, 
	integer *ldr, integer *iwork, doublereal *dwork, integer *ldwork, 
	logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, r_dim1, r_offset, s_dim1, s_offset, i__1, i__2, i__3, 
	    i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static integer i__, n2, kd, kg, kq, ks, ku, kw, lw, kbw, kcw, kdw, kwi, 
	    kwr, ierr, ktau;
    static doublereal rtol;
    extern /* Subroutine */ int sb02md_(char *, char *, char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, logical *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), dsyrk_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgerqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dorgrq_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
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

/*     To compute the Cholesky factors Su and Ru of the controllability */
/*     Grammian P = Su*Su' and observability Grammian Q = Ru'*Ru, */
/*     respectively, satisfying */

/*            A*P  + P*A' +  scalec^2*B*B'   = 0,       (1) */

/*            A'*Q + Q*A  +  scaleo^2*Cw'*Cw = 0,       (2) */

/*     where */
/*            Cw = Hw - Bw'*X, */
/*            Hw = inv(Dw)*C, */
/*            Bw = (B*D' + P*C')*inv(Dw'), */
/*            D*D' = Dw*Dw' (Dw upper triangular), */

/*     and, with Aw = A - Bw*Hw, X is the stabilizing solution of the */
/*     Riccati equation */

/*            Aw'*X + X*Aw + Hw'*Hw + X*Bw*Bw'*X = 0.   (3) */

/*     The P-by-M matrix D must have full row rank. Matrix A must be */
/*     stable and in a real Schur form. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of state-space representation, i.e., */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  M >= P >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             stable state dynamics matrix A in a real Schur canonical */
/*             form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B, corresponding to the Schur matrix A. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C, corresponding to the Schur */
/*             matrix A. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must */
/*             contain the full row rank input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     SCALEC  (output) DOUBLE PRECISION */
/*             Scaling factor for the controllability Grammian in (1). */

/*     SCALEO  (output) DOUBLE PRECISION */
/*             Scaling factor for the observability Grammian in (2). */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor Su of the cotrollability */
/*             Grammian P = Su*Su' satisfying (1). */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor Ru of the observability */
/*             Grammian Q = Ru'*Ru satisfying (2). */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*N */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) contains RCOND, the reciprocal */
/*             condition number of the U11 matrix from the expression */
/*             used to compute X = U21*inv(U11). A small value RCOND */
/*             indicates possible ill-conditioning of the Riccati */
/*             equation (3). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 2, N*(MAX(N,M,P)+5), */
/*                            2*N*P+MAX(P*(M+2),10*N*(N+1) ) ). */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension 2*N */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the state matrix A is not stable or is not in a */
/*                   real Schur form; */
/*             = 2:  the reduction of Hamiltonian matrix to real Schur */
/*                   form failed; */
/*             = 3:  the reordering of the real Schur form of the */
/*                   Hamiltonian matrix failed; */
/*             = 4:  the Hamiltonian matrix has less than N stable */
/*                   eigenvalues; */
/*             = 5:  the coefficient matrix U11 in the linear system */
/*                   X*U11 = U21, used to determine X, is singular to */
/*                   working precision; */
/*             = 6:  the feedthrough matrix D has not a full row rank P. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000. */
/*     D. Sima, University of Bucharest, May 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2000. */
/*     Based on the RASP routines SRGRO and SRGRO1, by A. Varga, 1992. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001. */

/*     KEYWORDS */

/*     Minimal realization, model reduction, multivariable system, */
/*     state-space model, state-space representation, */
/*     stochastic balancing. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    *info = 0;
/* Computing MAX */
/* Computing MAX */
    i__3 = max(*n,*m);
/* Computing MAX */
    i__4 = *p * (*m + 2), i__5 = *n * 10 * (*n + 1);
    i__1 = 2, i__2 = *n * (max(i__3,*p) + 5), i__1 = max(i__1,i__2), i__2 = (*
	    n << 1) * *p + max(i__4,i__5);
    lw = max(i__1,i__2);

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0 || *p > *m) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldc < max(1,*p)) {
	*info = -9;
    } else if (*ldd < max(1,*p)) {
	*info = -11;
    } else if (*lds < max(1,*n)) {
	*info = -15;
    } else if (*ldr < max(1,*n)) {
	*info = -17;
    } else if (*ldwork < lw) {
	*info = -20;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09HY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *scalec = 1.;
    *scaleo = 1.;
/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	dwork[1] = 2.;
	dwork[2] = 1.;
	return 0;
    }

/*     Solve for Su the Lyapunov equation */
/*                                      2 */
/*     A*(Su*Su') + (Su*Su')*A' + scalec *B*B' = 0 . */

/*     Workspace:  need   N*(MAX(N,M) + 5); */
/*                 prefer larger. */

    ku = 1;
    ktau = ku + *n * max(*n,*m);
    kw = ktau + *n;

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
    i__1 = *ldwork - kw + 1;
    sb03ou_(&c_false, &c_true, n, m, &a[a_offset], lda, &dwork[ku], n, &dwork[
	    ktau], &s[s_offset], lds, scalec, &dwork[kw], &i__1, &ierr);
    if (ierr != 0) {
	*info = 1;
	return 0;
    }
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Allocate workspace for Bw' (P*N), Cw (P*N), Q2 (P*M), */
/*     where Q2 = inv(Dw)*D. */
/*     Workspace:  need   2*N*P + P*M. */

    kbw = 1;
    kcw = kbw + *p * *n;
    kd = kcw + *p * *n;
    kdw = kd + *p * (*m - *p);
    ktau = kd + *p * *m;
    kw = ktau + *p;

/*     Compute an upper-triangular Dw such that D*D' = Dw*Dw', using */
/*     the RQ-decomposition of D: D = [0 Dw]*( Q1 ). */
/*                                           ( Q2 ) */
/*     Additional workspace:  need 2*P; prefer P + P*NB. */

    dlacpy_("F", p, m, &d__[d_offset], ldd, &dwork[kd], p, (ftnlen)1);
    i__1 = *ldwork - kw + 1;
    dgerqf_(p, m, &dwork[kd], p, &dwork[ktau], &dwork[kw], &i__1, &ierr);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__1,i__2);

/*     Check the full row rank of D. */

    rtol = (doublereal) (*m) * dlamch_("E", (ftnlen)1) * dlange_("1", p, m, &
	    d__[d_offset], ldd, &dwork[1], (ftnlen)1);
    i__1 = kdw + *p * *p - 1;
    i__2 = *p + 1;
    for (i__ = kdw; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if ((d__1 = dwork[i__], abs(d__1)) <= rtol) {
	    *info = 6;
	    return 0;
	}
/* L10: */
    }
/*                    -1 */
/*     Compute Hw = Dw  *C. */

    dlacpy_("F", p, n, &c__[c_offset], ldc, &dwork[kcw], p, (ftnlen)1);
    dtrsm_("Left", "Upper", "No-transpose", "Non-unit", p, n, &c_b15, &dwork[
	    kdw], p, &dwork[kcw], p, (ftnlen)4, (ftnlen)5, (ftnlen)12, (
	    ftnlen)8);

/*     Compute Bw' = inv(Dw)*(D*B' + C*Su*Su'). */

/*     Compute first Hw*Su*Su' in Bw'. */

    dlacpy_("F", p, n, &dwork[kcw], p, &dwork[kbw], p, (ftnlen)1);
    dtrmm_("Right", "Upper", "No-transpose", "Non-unit", p, n, &c_b15, &s[
	    s_offset], lds, &dwork[kbw], p, (ftnlen)5, (ftnlen)5, (ftnlen)12, 
	    (ftnlen)8);
    dtrmm_("Right", "Upper", "Transpose", "Non-unit", p, n, &c_b15, &s[
	    s_offset], lds, &dwork[kbw], p, (ftnlen)5, (ftnlen)5, (ftnlen)9, (
	    ftnlen)8);

/*     Compute Q2 = inv(Dw)*D, as the last P lines of the orthogonal */
/*     matrix ( Q1 ) from the RQ decomposition of D. */
/*            ( Q2 ) */
/*     Additional workspace:  need P; prefer P*NB. */

    i__2 = *ldwork - kw + 1;
    dorgrq_(p, m, p, &dwork[kd], p, &dwork[ktau], &dwork[kw], &i__2, &ierr);
/* Computing MAX */
    i__2 = wrkopt, i__1 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__2,i__1);

/*     Compute Bw' <- Bw' + Q2*B'. */

    dgemm_("No-transpose", "Transpose", p, n, m, &c_b15, &dwork[kd], p, &b[
	    b_offset], ldb, &c_b15, &dwork[kbw], p, (ftnlen)12, (ftnlen)9);

/*     Compute Aw = A - Bw*Hw in R. */

    dlacpy_("F", n, n, &a[a_offset], lda, &r__[r_offset], ldr, (ftnlen)1);
    dgemm_("Transpose", "No-transpose", n, n, p, &c_b34, &dwork[kbw], p, &
	    dwork[kcw], p, &c_b15, &r__[r_offset], ldr, (ftnlen)9, (ftnlen)12)
	    ;

/*     Allocate storage to solve the Riccati equation (3) for */
/*     G(N*N), Q(N*N), WR(2N), WI(2N), S(2N*2N), U(2N*2N). */

    n2 = *n + *n;
    kg = kd;
    kq = kg + *n * *n;
    kwr = kq + *n * *n;
    kwi = kwr + n2;
    ks = kwi + n2;
    ku = ks + n2 * n2;
    kw = ku + n2 * n2;

/*     Compute G = -Bw*Bw'. */

    dsyrk_("Upper", "Transpose", n, p, &c_b34, &dwork[kbw], p, &c_b39, &dwork[
	    kg], n, (ftnlen)5, (ftnlen)9);

/*     Compute Q = Hw'*Hw. */

    dsyrk_("Upper", "Transpose", n, p, &c_b15, &dwork[kcw], p, &c_b39, &dwork[
	    kq], n, (ftnlen)5, (ftnlen)9);

/*     Solve */

/*        Aw'*X + X*Aw + Q - X*G*X = 0, */

/*     with Q =  Hw'*Hw  and  G = -Bw*Bw'. */
/*     Additional workspace: need   6*N; */
/*                           prefer larger. */

    i__2 = *ldwork - kw + 1;
    sb02md_("Continuous", "None", "Upper", "General", "Stable", n, &r__[
	    r_offset], ldr, &dwork[kg], n, &dwork[kq], n, &rcond, &dwork[kwr],
	     &dwork[kwi], &dwork[ks], &n2, &dwork[ku], &n2, &iwork[1], &dwork[
	    kw], &i__2, &bwork[1], info, (ftnlen)10, (ftnlen)4, (ftnlen)5, (
	    ftnlen)7, (ftnlen)6);
    if (*info != 0) {
	return 0;
    }
/* Computing MAX */
    i__2 = wrkopt, i__1 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__2,i__1);

/*     Compute Cw = Hw - Bw'*X. */

    dgemm_("No-transpose", "No-transpose", p, n, n, &c_b34, &dwork[kbw], p, &
	    dwork[kq], n, &c_b15, &dwork[kcw], p, (ftnlen)12, (ftnlen)12);

/*     Solve for Ru the Lyapunov equation */
/*                                      2 */
/*     A'*(Ru'*Ru) + (Ru'*Ru)*A + scaleo  * Cw'*Cw = 0 . */

/*     Workspace:  need   N*(MAX(N,P) + 5); */
/*                 prefer larger. */

    ktau = kcw + *n * max(*n,*p);
    kw = ktau + *n;

    i__2 = *ldwork - kw + 1;
    sb03ou_(&c_false, &c_false, n, p, &a[a_offset], lda, &dwork[kcw], p, &
	    dwork[ktau], &r__[r_offset], ldr, scaleo, &dwork[kw], &i__2, &
	    ierr);
/* Computing MAX */
    i__2 = wrkopt, i__1 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__2,i__1);

/*     Save optimal workspace and RCOND. */

    dwork[1] = (doublereal) wrkopt;
    dwork[2] = rcond;

    return 0;
/* *** Last line of AB09HY *** */
} /* ab09hy_ */

