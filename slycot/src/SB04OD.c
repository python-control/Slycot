/* SB04OD.f -- translated by f2c (version 20100827).
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
static doublereal c_b52 = 1.;
static doublereal c_b53 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb04od_(char *reduce, char *trans, char *jobd, integer *
	m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *p, integer *ldp, doublereal *q, 
	integer *ldq, doublereal *u, integer *ldu, doublereal *v, integer *
	ldv, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen reduce_len, ftnlen trans_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, p_dim1, p_offset, 
	    q_dim1, q_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, mn, ijob;
    static doublereal anrm, bnrm, dnrm;
    static integer ierr;
    static doublereal enrm;
    static logical ljob1, ljob2;
    extern /* Subroutine */ int dgegs_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical ljobd, ljobf;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical ljobdf;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical ilascl, ilbscl, ildscl, ilescl, lredra, lredrb, lredua, 
	    lredub, lreduc;
    static doublereal bignum, safmax, safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lredur, ltrann;
    static doublereal anrmto, bnrmto, dnrmto, enrmto;
    extern /* Subroutine */ int dtgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static integer minwrk;
    static doublereal smlnum;
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

/*     To solve for R and L one of the generalized Sylvester equations */

/*        A * R - L * B = scale * C ) */
/*                                  )                                 (1) */
/*        D * R - L * E = scale * F ) */

/*     or */

/*        A' * R + D' * L = scale * C    ) */
/*                                       )                            (2) */
/*        R * B' + L * E' = scale * (-F) ) */

/*     where A and D are M-by-M matrices, B and E are N-by-N matrices and */
/*     C, F, R and L are M-by-N matrices. */

/*     The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an */
/*     output scaling factor chosen to avoid overflow. */

/*     The routine also optionally computes a Dif estimate, which */
/*     measures the separation of the spectrum of the matrix pair (A,D) */
/*     from the spectrum of the matrix pair (B,E), Dif[(A,D),(B,E)]. */

/*     ARGUMENTS */

/*     MODE PARAMETERS */

/*     REDUCE  CHARACTER*1 */
/*             Indicates whether the matrix pairs (A,D) and/or (B,E) are */
/*             to be reduced to generalized Schur form as follows: */
/*             = 'R':  The matrix pairs (A,D) and (B,E) are to be reduced */
/*                     to generalized (real) Schur canonical form; */
/*             = 'A':  The matrix pair (A,D) only is to be reduced */
/*                     to generalized (real) Schur canonical form, */
/*                     and the matrix pair (B,E) already is in this form; */
/*             = 'B':  The matrix pair (B,E) only is to be reduced */
/*                     to generalized (real) Schur canonical form, */
/*                     and the matrix pair (A,D) already is in this form; */
/*             = 'N':  The matrix pairs (A,D) and (B,E) are already in */
/*                     generalized (real) Schur canonical form, as */
/*                     produced by LAPACK routine DGEES. */

/*     TRANS   CHARACTER*1 */
/*             Indicates which of the equations, (1) or (2), is to be */
/*             solved as follows: */
/*             = 'N':  The generalized Sylvester equation (1) is to be */
/*                     solved; */
/*             = 'T':  The "transposed" generalized Sylvester equation */
/*                     (2) is to be solved. */

/*     JOBD    CHARACTER*1 */
/*             Indicates whether the Dif estimator is to be computed as */
/*             follows: */
/*             = '1':  Only the one-norm-based Dif estimate is computed */
/*                     and stored in DIF; */
/*             = '2':  Only the Frobenius norm-based Dif estimate is */
/*                     computed and stored in DIF; */
/*             = 'D':  The equation (1) is solved and the one-norm-based */
/*                     Dif estimate is computed and stored in DIF; */
/*             = 'F':  The equation (1) is solved and the Frobenius norm- */
/*                     based Dif estimate is computed and stored in DIF; */
/*             = 'N':  The Dif estimator is not required and hence DIF is */
/*                     not referenced. (Solve either (1) or (2) only.) */
/*             JOBD is not referenced if TRANS = 'T'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrices A and D and the number of rows */
/*             of the matrices C, F, R and L.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrices B and E and the number of */
/*             columns of the matrices C, F, R and L.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the coefficient matrix A of the equation; A must */
/*             be in upper quasi-triangular form if REDUCE = 'B' or 'N'. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the upper quasi-triangular form of A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix B of the equation; B must */
/*             be in upper quasi-triangular form if REDUCE = 'A' or 'N'. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper quasi-triangular form of B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand side matrix C of the first equation */
/*             in (1) or (2). */
/*             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N */
/*             part of this array contains the solution matrix R of the */
/*             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading */
/*             M-by-N part of this array contains the solution matrix R */
/*             achieved during the computation of the Dif estimate. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the coefficient matrix D of the equation; D must */
/*             be in upper triangular form if REDUCE = 'B' or 'N'. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the upper triangular form of D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix E of the equation; E must */
/*             be in upper triangular form if REDUCE = 'A' or 'N'. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper triangular form of E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand side matrix F of the second */
/*             equation in (1) or (2). */
/*             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N */
/*             part of this array contains the solution matrix L of the */
/*             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading */
/*             M-by-N part of this array contains the solution matrix L */
/*             achieved during the computation of the Dif estimate. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scaling factor in (1) or (2). If 0 < SCALE < 1, C and */
/*             F hold the solutions R and L, respectively, to a slightly */
/*             perturbed system (but the input or computed generalized */
/*             (real) Schur canonical form matrices A, B, D, and E */
/*             have not been changed). If SCALE = 0, C and F hold the */
/*             solutions R and L, respectively, to the homogeneous system */
/*             with C = F = 0. Normally, SCALE = 1. */

/*     DIF     (output) DOUBLE PRECISION */
/*             If TRANS = 'N' and JOBD <> 'N', then DIF contains the */
/*             value of the Dif estimator, which is an upper bound of */
/*                                                    -1 */
/*             Dif[(A,D),(B,E)] = sigma_min(Z) = 1/||Z  ||, in either the */
/*             one-norm, or Frobenius norm, respectively (see METHOD). */
/*             Otherwise, DIF is not referenced. */

/*     P       (output) DOUBLE PRECISION array, dimension (LDP,*) */
/*             If REDUCE = 'R' or 'A', then the leading M-by-M part of */
/*             this array contains the (left) transformation matrix used */
/*             to reduce (A,D) to generalized Schur form. */
/*             Otherwise, P is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDP = 1 and declare this */
/*             array to be P(1,1) in the calling program). */

/*     LDP     INTEGER */
/*             The leading dimension of array P. */
/*             LDP >= MAX(1,M) if REDUCE = 'R' or 'A', */
/*             LDP >= 1        if REDUCE = 'B' or 'N'. */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             If REDUCE = 'R' or 'A', then the leading M-by-M part of */
/*             this array contains the (right) transformation matrix used */
/*             to reduce (A,D) to generalized Schur form. */
/*             Otherwise, Q is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDQ = 1 and declare this */
/*             array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,M) if REDUCE = 'R' or 'A', */
/*             LDQ >= 1        if REDUCE = 'B' or 'N'. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             If REDUCE = 'R' or 'B', then the leading N-by-N part of */
/*             this array contains the (left) transformation matrix used */
/*             to reduce (B,E) to generalized Schur form. */
/*             Otherwise, U is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDU = 1 and declare this */
/*             array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= MAX(1,N) if REDUCE = 'R' or 'B', */
/*             LDU >= 1        if REDUCE = 'A' or 'N'. */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             If REDUCE = 'R' or 'B', then the leading N-by-N part of */
/*             this array contains the (right) transformation matrix used */
/*             to reduce (B,E) to generalized Schur form. */
/*             Otherwise, V is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDV = 1 and declare this */
/*             array to be V(1,1) in the calling program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. */
/*             LDV >= MAX(1,N) if REDUCE = 'R' or 'B', */
/*             LDV >= 1        if REDUCE = 'A' or 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+N+6) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If TRANS = 'N' and JOBD = 'D' or 'F', then */
/*                LDWORK = MAX(1,7*M,7*N,2*M*N) if REDUCE = 'R'; */
/*                LDWORK = MAX(1,7*M,2*M*N)     if REDUCE = 'A'; */
/*                LDWORK = MAX(1,7*N,2*M*N)     if REDUCE = 'B'; */
/*                LDWORK = MAX(1,2*M*N)         if REDUCE = 'N'. */
/*             Otherwise, the term 2*M*N above should be omitted. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if REDUCE <> 'N' and either (A,D) and/or (B,E) */
/*                   cannot be reduced to generalized Schur form; */
/*             = 2:  if REDUCE = 'N' and either A or B is not in */
/*                   upper quasi-triangular form; */
/*             = 3:  if a singular matrix was encountered during the */
/*                   computation of the solution matrices R and L, that */
/*                   is (A,D) and (B,E) have common or close eigenvalues. */

/*     METHOD */

/*     For the case TRANS = 'N', and REDUCE = 'R' or 'N', the algorithm */
/*     used by the routine consists of four steps (see [1] and [2]) as */
/*     follows: */

/*        (a) if REDUCE = 'R', then the matrix pairs (A,D) and (B,E) are */
/*            transformed to generalized Schur form, i.e. orthogonal */
/*            matrices P, Q, U and V are computed such that P' * A * Q */
/*            and U' * B * V are in upper quasi-triangular form and */
/*            P' * D * Q and U' * E * V are in upper triangular form; */
/*        (b) if REDUCE = 'R', then the matrices C and F are transformed */
/*            to give P' * C * V and P' * F * V respectively; */
/*        (c) if REDUCE = 'R', then the transformed system */

/*            P' * A * Q * R1 - L1 * U' * B * V = scale * P' * C * V */
/*            P' * D * Q * R1 - L1 * U' * E * V = scale * P' * F * V */

/*            is solved to give R1 and L1; otherwise, equation (1) is */
/*            solved to give R and L directly. The Dif estimator */
/*            is also computed if JOBD <> 'N'. */
/*        (d) if REDUCE = 'R', then the solution is transformed back */
/*            to give R = Q * R1 * V' and L = P * L1 * U'. */

/*     By using Kronecker products, equation (1) can also be written as */
/*     the system of linear equations Z * x = scale*y (see [1]), where */

/*            | I*A    I*D  | */
/*        Z = |             |. */
/*            |-B'*I  -E'*I | */

/*                                              -1 */
/*     If JOBD <> 'N', then a lower bound on ||Z  ||, in either the one- */
/*     norm or Frobenius norm, is computed, which in most cases is */
/*     a reliable estimate of the true value. Notice that since Z is a */
/*     matrix of order 2 * M * N, the exact value of Dif (i.e., in the */
/*     Frobenius norm case, the smallest singular value of Z) may be very */
/*     expensive to compute. */

/*     The case TRANS = 'N', and REDUCE = 'A' or 'B', is similar, but */
/*     only one of the matrix pairs should be reduced and the */
/*     calculations simplify. */

/*     For the case TRANS = 'T', and REDUCE = 'R' or 'N', the algorithm */
/*     is similar, but the steps (b), (c), and (d) are as follows: */

/*        (b) if REDUCE = 'R', then the matrices C and F are transformed */
/*            to give Q' * C * V and P' * F * U respectively; */
/*        (c) if REDUCE = 'R', then the transformed system */

/*            Q' * A' * P * R1 + Q' * D' * P * L1 =  scale * Q' * C * V */
/*            R1 * V' * B' * U + L1 * V' * E' * U = -scale * P' * F * U */

/*            is solved to give R1 and L1; otherwise, equation (2) is */
/*            solved to give R and L directly. */
/*        (d) if REDUCE = 'R', then the solution is transformed back */
/*            to give R = P * R1 * V' and L = P * L1 * V'. */

/*     REFERENCES */

/*     [1] Kagstrom, B. and Westin, L. */
/*         Generalized Schur Methods with Condition Estimators for */
/*         Solving the Generalized Sylvester Equation. */
/*         IEEE Trans. Auto. Contr., 34, pp. 745-751, 1989. */
/*     [2] Kagstrom, B. and Westin, L. */
/*         GSYLV - Fortran Routines for the Generalized Schur Method with */
/*         Dif Estimators for Solving the Generalized Sylvester */
/*         Equation. */
/*         Report UMINF-132.86, Institute of Information Processing, */
/*         Univ. of Umea, Sweden, July 1987. */
/*     [3] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur Method for the Problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */
/*     [4] Kagstrom, B. and Van Dooren, P. */
/*         Additive Decomposition of a Transfer Function with respect to */
/*         a Specified Region. */
/*         In: "Signal Processing, Scattering and Operator Theory, and */
/*         Numerical Methods" (Eds. M.A. Kaashoek et al.). */
/*         Proceedings of MTNS-89, Vol. 3, pp. 469-477, Birkhauser Boston */
/*         Inc., 1990. */
/*     [5] Kagstrom, B. and Van Dooren, P. */
/*         A Generalized State-space Approach for the Additive */
/*         Decomposition of a Transfer Matrix. */
/*         Report UMINF-91.12, Institute of Information Processing, Univ. */
/*         of Umea, Sweden, April 1991. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. A reliable estimate for the */
/*     condition number of Z in the Frobenius norm, is (see [1]) */

/*        K(Z) = SQRT(  ||A||**2 + ||B||**2 + ||C||**2 + ||D||**2 )/DIF. */

/*     If mu is an upper bound on the relative error of the elements of */
/*     the matrices A, B, C, D, E and F, then the relative error in the */
/*     actual solution is approximately mu * K(Z). */

/*     The relative error in the computed solution (due to rounding */
/*     errors) is approximately EPS * K(Z), where EPS is the machine */
/*     precision (see LAPACK Library routine DLAMCH). */

/*     FURTHER COMMENTS */

/*     For applications of the generalized Sylvester equation in control */
/*     theory, see [4] and [5]. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB04CD by Bo Kagstrom and Lars */
/*     Westin. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, Dec. 1999, */
/*     May 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, orthogonal transformation, real */
/*     Schur form, Sylvester equation. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    p_dim1 = *ldp;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    mn = max(*m,*n);
    lredur = lsame_(reduce, "R", (ftnlen)1, (ftnlen)1);
    lredua = lsame_(reduce, "A", (ftnlen)1, (ftnlen)1);
    lredub = lsame_(reduce, "B", (ftnlen)1, (ftnlen)1);
    lredra = lredur || lredua;
    lredrb = lredur || lredub;
    lreduc = lredra || lredub;
    if (lredur) {
/* Computing MAX */
	i__1 = 1, i__2 = mn * 7;
	minwrk = max(i__1,i__2);
    } else if (lredua) {
/* Computing MAX */
	i__1 = 1, i__2 = *m * 7;
	minwrk = max(i__1,i__2);
    } else if (lredub) {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 7;
	minwrk = max(i__1,i__2);
    } else {
	minwrk = 1;
    }
    ltrann = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    if (ltrann) {
	ljob1 = lsame_(jobd, "1", (ftnlen)1, (ftnlen)1);
	ljob2 = lsame_(jobd, "2", (ftnlen)1, (ftnlen)1);
	ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
	ljobf = lsame_(jobd, "F", (ftnlen)1, (ftnlen)1);
	ljobdf = ljob1 || ljob2 || ljobd || ljobf;
	if (ljobd || ljobf) {
/* Computing MAX */
	    i__1 = minwrk, i__2 = (*m << 1) * *n;
	    minwrk = max(i__1,i__2);
	}
    }

/*     Test the input scalar arguments. */

    if (! lreduc && ! lsame_(reduce, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ltrann && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (ltrann) {
	if (! ljobdf && ! lsame_(jobd, "N", (ftnlen)1, (ftnlen)1)) {
	    *info = -3;
	}
    }
    if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*m)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldc < max(1,*m)) {
	*info = -11;
    } else if (*ldd < max(1,*m)) {
	*info = -13;
    } else if (*lde < max(1,*n)) {
	*info = -15;
    } else if (*ldf < max(1,*m)) {
	*info = -17;
    } else if (! lredra && *ldp < 1 || lredra && *ldp < max(1,*m)) {
	*info = -21;
    } else if (! lredra && *ldq < 1 || lredra && *ldq < max(1,*m)) {
	*info = -23;
    } else if (! lredrb && *ldu < 1 || lredrb && *ldu < max(1,*n)) {
	*info = -25;
    } else if (! lredrb && *ldv < 1 || lredrb && *ldv < max(1,*n)) {
	*info = -27;
    } else if (*ldwork < minwrk) {
	*info = -30;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB04OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0) {
	*scale = 1.;
	dwork[1] = 1.;
	if (ltrann) {
	    if (ljobdf) {
		*dif = 1.;
	    }
	}
	return 0;
    }
    wrkopt = 1;
    sufwrk = *ldwork >= *m * *n;

/*     STEP 1: Reduce (A,D) and/or (B,E) to generalized Schur form. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    if (lreduc) {

/*        Get machine constants. */

	safmin = dlamch_("Safe minimum", (ftnlen)12);
	safmax = 1. / safmin;
	dlabad_(&safmin, &safmax);
	smlnum = sqrt(safmin) / dlamch_("Precision", (ftnlen)9);
	bignum = 1. / smlnum;

	if (! lredub) {

/*           Scale A if max element outside range [SMLNUM,BIGNUM]. */

	    anrm = dlange_("M", m, m, &a[a_offset], lda, &dwork[1], (ftnlen)1)
		    ;
	    ilascl = FALSE_;
	    if (anrm > 0. && anrm < smlnum) {
		anrmto = smlnum;
		ilascl = TRUE_;
	    } else if (anrm > bignum) {
		anrmto = bignum;
		ilascl = TRUE_;
	    }
	    if (ilascl) {
		dlascl_("G", &c__0, &c__0, &anrm, &anrmto, m, m, &a[a_offset],
			 lda, &ierr, (ftnlen)1);
	    }

/*           Scale D if max element outside range [SMLNUM,BIGNUM] */

	    dnrm = dlange_("M", m, m, &d__[d_offset], ldd, &dwork[1], (ftnlen)
		    1);
	    ildscl = FALSE_;
	    if (dnrm > 0. && dnrm < smlnum) {
		dnrmto = smlnum;
		ildscl = TRUE_;
	    } else if (dnrm > bignum) {
		dnrmto = bignum;
		ildscl = TRUE_;
	    }
	    if (ildscl) {
		dlascl_("G", &c__0, &c__0, &dnrm, &dnrmto, m, m, &d__[
			d_offset], ldd, &ierr, (ftnlen)1);
	    }

/*           Reduce (A,D) to generalized Schur form. */
/*           Workspace:  need   7*M; */
/*                       prefer 5*M + M*(NB+1). */

	    i__1 = *ldwork - *m * 3;
	    dgegs_("Vectors left", "Vectors right", m, &a[a_offset], lda, &
		    d__[d_offset], ldd, &dwork[1], &dwork[*m + 1], &dwork[(*m 
		    << 1) + 1], &p[p_offset], ldp, &q[q_offset], ldq, &dwork[*
		    m * 3 + 1], &i__1, info, (ftnlen)12, (ftnlen)13);

/*           Undo scaling */

	    if (ilascl) {
		dlascl_("H", &c__0, &c__0, &anrmto, &anrm, m, m, &a[a_offset],
			 lda, &ierr, (ftnlen)1);
	    }

	    if (ildscl) {
		dlascl_("U", &c__0, &c__0, &dnrmto, &dnrm, m, m, &d__[
			d_offset], ldd, &ierr, (ftnlen)1);
	    }

	    if (*info != 0) {
		*info = 1;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*m * 3 + 1] + *m * 3;
	    wrkopt = max(i__1,i__2);
	}
	if (! lredua) {

/*           Scale B if max element outside range [SMLNUM,BIGNUM] */

	    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &dwork[1], (ftnlen)1)
		    ;
	    ilbscl = FALSE_;
	    if (bnrm > 0. && bnrm < smlnum) {
		bnrmto = smlnum;
		ilbscl = TRUE_;
	    } else if (bnrm > bignum) {
		bnrmto = bignum;
		ilbscl = TRUE_;
	    }
	    if (ilbscl) {
		dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset],
			 ldb, &ierr, (ftnlen)1);
	    }

/*           Scale E if max element outside range [SMLNUM,BIGNUM] */

	    enrm = dlange_("M", n, n, &e[e_offset], lde, &dwork[1], (ftnlen)1)
		    ;
	    ilescl = FALSE_;
	    if (enrm > 0. && enrm < smlnum) {
		enrmto = smlnum;
		ilescl = TRUE_;
	    } else if (enrm > bignum) {
		enrmto = bignum;
		ilescl = TRUE_;
	    }
	    if (ilescl) {
		dlascl_("G", &c__0, &c__0, &enrm, &enrmto, n, n, &e[e_offset],
			 lde, &ierr, (ftnlen)1);
	    }

/*           Reduce (B,E) to generalized Schur form. */
/*           Workspace:  need   7*N; */
/*                       prefer 5*N + N*(NB+1). */

	    i__1 = *ldwork - *n * 3;
	    dgegs_("Vectors left", "Vectors right", n, &b[b_offset], ldb, &e[
		    e_offset], lde, &dwork[1], &dwork[*n + 1], &dwork[(*n << 
		    1) + 1], &u[u_offset], ldu, &v[v_offset], ldv, &dwork[*n *
		     3 + 1], &i__1, info, (ftnlen)12, (ftnlen)13);

/*           Undo scaling */

	    if (ilbscl) {
		dlascl_("H", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset],
			 ldb, &ierr, (ftnlen)1);
	    }

	    if (ilescl) {
		dlascl_("U", &c__0, &c__0, &enrmto, &enrm, n, n, &e[e_offset],
			 lde, &ierr, (ftnlen)1);
	    }

	    if (*info != 0) {
		*info = 1;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*n * 3 + 1] + *n * 3;
	    wrkopt = max(i__1,i__2);
	}
    }

    if (! lredur) {

/*        Set INFO = 2 if A and/or B are/is not in quasi-triangular form. */

	if (! lredua) {
	    i__ = 1;

L20:
	    if (i__ <= *m - 2) {
		if (a[i__ + 1 + i__ * a_dim1] != 0.) {
		    if (a[i__ + 2 + (i__ + 1) * a_dim1] != 0.) {
			*info = 2;
			return 0;
		    } else {
			++i__;
		    }
		}
		++i__;
		goto L20;
	    }
	}

	if (! lredub) {
	    i__ = 1;

L40:
	    if (i__ <= *n - 2) {
		if (b[i__ + 1 + i__ * b_dim1] != 0.) {
		    if (b[i__ + 2 + (i__ + 1) * b_dim1] != 0.) {
			*info = 2;
			return 0;
		    } else {
			++i__;
		    }
		}
		++i__;
		goto L40;
	    }
	}
    }

/*     STEP 2: Modify right hand sides (C,F). */

    if (lreduc) {
/* Computing MAX */
	i__1 = wrkopt, i__2 = *m * *n;
	wrkopt = max(i__1,i__2);
	if (sufwrk) {

/*           Enough workspace for a BLAS 3 calculation. */

	    if (ltrann) {

/*              Equation (1). */

		if (! lredub) {
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)9, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
		}
		if (! lredub) {
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &f[f_offset], ldf, &c_b53, &dwork[
			    1], m, (ftnlen)9, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
		}
	    } else {

/*              Equation (2). */

		if (! lredub) {
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &q[
			    q_offset], ldq, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)9, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
		}
		if (! lredub) {
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &f[f_offset], ldf, &c_b53, &dwork[
			    1], m, (ftnlen)9, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &u[u_offset], ldu, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
		}
	    }
	} else {

/*           Use a BLAS 2 calculation. */

	    if (ltrann) {

/*              Equation (1). */

		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				c__[i__ * c_dim1 + 1], &c__1, &c_b53, &dwork[
				1], &c__1, (ftnlen)9);
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
/* L60: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				c__[i__ + c_dim1], ldc, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
/* L80: */
		    }

		}
		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				f[i__ * f_dim1 + 1], &c__1, &c_b53, &dwork[1],
				 &c__1, (ftnlen)9);
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
/* L100: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				f[i__ + f_dim1], ldf, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
/* L120: */
		    }

		}
	    } else {

/*              Equation (2). */

		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", m, m, &c_b52, &q[q_offset], ldq, &
				c__[i__ * c_dim1 + 1], &c__1, &c_b53, &dwork[
				1], &c__1, (ftnlen)9);
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
/* L140: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				c__[i__ + c_dim1], ldc, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
/* L160: */
		    }

		}
		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				f[i__ * f_dim1 + 1], &c__1, &c_b53, &dwork[1],
				 &c__1, (ftnlen)9);
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
/* L180: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("Transpose", n, n, &c_b52, &u[u_offset], ldu, &
				f[i__ + f_dim1], ldf, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
/* L200: */
		    }

		}
	    }
	}
    }

/*     STEP 3: Solve the transformed system and compute the Dif */
/*             estimator. */

    if (ltrann) {
	if (ljobd) {
	    ijob = 1;
	} else if (ljobf) {
	    ijob = 2;
	} else if (ljob1) {
	    ijob = 3;
	} else if (ljob2) {
	    ijob = 4;
	} else {
	    ijob = 0;
	}
    } else {
	ijob = 0;
    }

/*     Workspace:  need 2*M*N if TRANS = 'N' and JOBD = 'D' or 'F'; */
/*                      1, otherwise. */

    dtgsyl_(trans, &ijob, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], lde, &f[
	    f_offset], ldf, scale, dif, &dwork[1], ldwork, &iwork[1], info, (
	    ftnlen)1);
    if (*info != 0) {
	*info = 3;
	return 0;
    }
    if (ltrann) {
	if (ljobd || ljobf) {
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (*m << 1) * *n;
	    wrkopt = max(i__1,i__2);
	}
    }

/*     STEP 4: Back transformation of the solution. */

    if (lreduc) {
	if (sufwrk) {

/*           Enough workspace for a BLAS 3 calculation. */

	    if (ltrann) {

/*              Equation (1). */

		if (! lredub) {
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    q[q_offset], ldq, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)9);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
		}
		if (! lredub) {
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &f[f_offset], ldf, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &u[u_offset], ldu, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)9);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
		}
	    } else {

/*              Equation (2). */

		if (! lredub) {
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)9);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
		}
		if (! lredub) {
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &f[f_offset], ldf, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
		} else {
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
		}
		if (! lredua) {
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)9);
		} else {
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
		}
	    }
	} else {

/*           Use a BLAS 2 calculation. */

	    if (ltrann) {

/*              Equation (1). */

		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", m, m, &c_b52, &q[q_offset], 
				ldq, &c__[i__ * c_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
/* L220: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &c__[i__ + c_dim1], ldc, &c_b53, &dwork[
				1], &c__1, (ftnlen)12);
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
/* L240: */
		    }

		}
		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &f[i__ * f_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
/* L260: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", n, n, &c_b52, &u[u_offset], 
				ldu, &f[i__ + f_dim1], ldf, &c_b53, &dwork[1],
				 &c__1, (ftnlen)12);
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
/* L280: */
		    }

		}
	    } else {

/*              Equation (2). */

		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &c__[i__ * c_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
/* L300: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &c__[i__ + c_dim1], ldc, &c_b53, &dwork[
				1], &c__1, (ftnlen)12);
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
/* L320: */
		    }

		}
		if (! lredub) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &f[i__ * f_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
/* L340: */
		    }

		}
		if (! lredua) {

		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &f[i__ + f_dim1], ldf, &c_b53, &dwork[1],
				 &c__1, (ftnlen)12);
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
/* L360: */
		    }

		}
	    }
	}
    }

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of SB04OD *** */
} /* sb04od_ */

