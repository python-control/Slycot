/* SG03BD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static integer c__2 = 2;
static doublereal c_b19 = 1.;
static integer c__1 = 1;
static doublereal c_b55 = -1.;

/* Subroutine */ int sg03bd_(char *dico, char *fact, char *trans, integer *n, 
	integer *m, doublereal *a, integer *lda, doublereal *e, integer *lde, 
	doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	doublereal *b, integer *ldb, doublereal *scale, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen dico_len, ftnlen fact_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, q_dim1, 
	    q_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal e1[4]	/* was [2][2] */, s1, s2, wi, wr1, wr2;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgegs_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), dgemm_(char *
	    , char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sg03bu_(char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen), sg03bv_(char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static integer minmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static logical isfact;
    extern /* Subroutine */ int dgerqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static doublereal safmin;
    static logical isdisc;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical istran;
    static integer minwrk, optwrk;


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

/*     To compute the Cholesky factor U of the matrix X, */

/*                 T */
/*        X = op(U)  * op(U), */

/*     which is the solution of either the generalized */
/*     c-stable continuous-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(E) + op(E)  * X * op(A) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (1) */

/*     or the generalized d-stable discrete-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(A) - op(E)  * X * op(E) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (2) */

/*     without first finding X and without the need to form the matrix */
/*     op(B)**T * op(B). */

/*     op(K) is either K or K**T for K = A, B, E, U. A and E are N-by-N */
/*     matrices, op(B) is an M-by-N matrix. The resulting matrix U is an */
/*     N-by-N upper triangular matrix with non-negative entries on its */
/*     main diagonal. SCALE is an output scale factor set to avoid */
/*     overflow in U. */

/*     In the continuous-time case (1) the pencil A - lambda * E must be */
/*     c-stable (that is, all eigenvalues must have negative real parts). */
/*     In the discrete-time case (2) the pencil A - lambda * E must be */
/*     d-stable (that is, the moduli of all eigenvalues must be smaller */
/*     than one). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies which type of the equation is considered: */
/*             = 'C':  Continuous-time equation (1); */
/*             = 'D':  Discrete-time equation (2). */

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

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of rows in the matrix op(B).  M >= 0. */

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

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N1) */
/*             On entry, if TRANS = 'T', the leading N-by-M part of this */
/*             array must contain the matrix B and N1 >= MAX(M,N). */
/*             If TRANS = 'N', the leading M-by-N part of this array */
/*             must contain the matrix B and N1 >= N. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the Cholesky factor U of the solution matrix X of the */
/*             problem, X = op(U)**T * op(U). */
/*             If M = 0 and N > 0, then U is set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             If TRANS = 'T', LDB >= MAX(1,N). */
/*             If TRANS = 'N', LDB >= MAX(1,M,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in U. */
/*             0 < SCALE <= 1. */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, 3, 5, 6, or 7, then */
/*             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the */
/*             eigenvalues of the matrix pencil A - lambda * E. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= MAX(1,4*N,6*N-6),  if FACT = 'N'; */
/*             LDWORK >= MAX(1,2*N,6*N-6),  if FACT = 'F'. */
/*             For good performance, LDWORK should be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the pencil A - lambda * E is (nearly) singular; */
/*                   perturbed values were used to solve the equation */
/*                   (but the reduced (quasi)triangular matrices A and E */
/*                   are unchanged); */
/*             = 2:  FACT = 'F' and the matrix contained in the upper */
/*                   Hessenberg part of the array A is not in upper */
/*                   quasitriangular form; */
/*             = 3:  FACT = 'F' and there is a 2-by-2 block on the main */
/*                   diagonal of the pencil A_s - lambda * E_s whose */
/*                   eigenvalues are not conjugate complex; */
/*             = 4:  FACT = 'N' and the pencil A - lambda * E cannot be */
/*                   reduced to generalized Schur form: LAPACK routine */
/*                   DGEGS has failed to converge; */
/*             = 5:  DICO = 'C' and the pencil A - lambda * E is not */
/*                   c-stable; */
/*             = 6:  DICO = 'D' and the pencil A - lambda * E is not */
/*                   d-stable; */
/*             = 7:  the LAPACK routine DSYEVX utilized to factorize M3 */
/*                   failed to converge in the discrete-time case (see */
/*                   section METHOD for SLICOT Library routine SG03BU). */
/*                   This error is unlikely to occur. */

/*     METHOD */

/*     An extension [2] of Hammarling's method [1] to generalized */
/*     Lyapunov equations is utilized to solve (1) or (2). */

/*     First the pencil A - lambda * E is reduced to real generalized */
/*     Schur form A_s - lambda * E_s by means of orthogonal */
/*     transformations (QZ-algorithm): */

/*        A_s = Q**T * A * Z   (upper quasitriangular)                (3) */

/*        E_s = Q**T * E * Z   (upper triangular).                    (4) */

/*     If the pencil A - lambda * E has already been factorized prior to */
/*     calling the routine however, then the factors A_s, E_s, Q and Z */
/*     may be supplied and the initial factorization omitted. */

/*     Depending on the parameters TRANS and M the N-by-N upper */
/*     triangular matrix B_s is defined as follows. In any case Q_B is */
/*     an M-by-M orthogonal matrix, which need not be accumulated. */

/*     1. If TRANS = 'N' and M < N, B_s is the upper triangular matrix */
/*        from the QR-factorization */

/*           ( Q_B  O )           ( B * Z ) */
/*           (        ) * B_s  =  (       ), */
/*           (  O   I )           (   O   ) */

/*        where the O's are zero matrices of proper size and I is the */
/*        identity matrix of order N-M. */

/*     2. If TRANS = 'N' and M >= N, B_s is the upper triangular matrix */
/*        from the (rectangular) QR-factorization */

/*                 ( B_s ) */
/*           Q_B * (     )  =  B * Z, */
/*                 (  O  ) */

/*        where O is the (M-N)-by-N zero matrix. */

/*     3. If TRANS = 'T' and M < N, B_s is the upper triangular matrix */
/*        from the RQ-factorization */

/*                       ( Q_B  O ) */
/*           (B_s  O ) * (        )  =  ( Q**T * B   O ). */
/*                       (  O   I ) */

/*     4. If TRANS = 'T' and M >= N, B_s is the upper triangular matrix */
/*        from the (rectangular) RQ-factorization */

/*           ( B_s   O ) * Q_B  =  Q**T * B, */

/*        where O is the N-by-(M-N) zero matrix. */

/*     Assuming SCALE = 1, the transformation of A, E and B described */
/*     above leads to the reduced continuous-time equation */

/*                 T        T */
/*          op(A_s)  op(U_s)  op(U_s) op(E_s) */

/*                 T        T */
/*        + op(E_s)  op(U_s)  op(U_s) op(A_s) */

/*                    T */
/*        =  - op(B_s)  op(B_s)                                       (5) */

/*     or to the reduced discrete-time equation */

/*                 T        T */
/*          op(A_s)  op(U_s)  op(U_s) op(A_s) */

/*                 T        T */
/*        - op(E_s)  op(U_s)  op(U_s) op(E_s) */

/*                    T */
/*        =  - op(B_s)  op(B_s).                                      (6) */

/*     For brevity we restrict ourself to equation (5) and the case */
/*     TRANS = 'N'. The other three cases can be treated in a similar */
/*     fashion. */

/*     We use the following partitioning for the matrices A_s, E_s, B_s */
/*     and U_s */

/*                 ( A11   A12 )          ( E11   E12 ) */
/*           A_s = (           ),   E_s = (           ), */
/*                 (   0   A22 )          (   0   E22 ) */

/*                 ( B11   B12 )          ( U11   U12 ) */
/*           B_s = (           ),   U_s = (           ).              (7) */
/*                 (   0   B22 )          (   0   U22 ) */

/*     The size of the (1,1)-blocks is 1-by-1 (iff A_s(2,1) = 0.0) or */
/*     2-by-2. */

/*     We compute U11 and U12**T in three steps. */

/*     Step I: */

/*        From (5) and (7) we get the 1-by-1 or 2-by-2 equation */

/*                T      T                   T      T */
/*             A11  * U11  * U11 * E11  + E11  * U11  * U11 * A11 */

/*                    T */
/*             = - B11  * B11. */

/*        For brevity, details are omitted here. See [2]. The technique */
/*        for computing U11 is similar to those applied to standard */
/*        Lyapunov equations in Hammarling's algorithm ([1], section 6). */

/*        Furthermore, the auxiliary matrices M1 and M2 defined as */
/*        follows */

/*                               -1      -1 */
/*           M1 = U11 * A11 * E11   * U11 */

/*                         -1      -1 */
/*           M2 = B11 * E11   * U11 */

/*        are computed in a numerically reliable way. */

/*     Step II: */

/*        The generalized Sylvester equation */

/*              T      T      T      T */
/*           A22  * U12  + E22  * U12  * M1  = */

/*                T           T      T      T      T */
/*           - B12  * M2 - A12  * U11  - E12  * U11  * M1 */

/*        is solved for U12**T. */

/*     Step III: */

/*        It can be shown that */

/*              T      T                  T      T */
/*           A22  * U22  * U22 * E22 + E22  * U22  * U22 * A22  = */

/*                T              T */
/*           - B22  * B22 - y * y                                     (8) */

/*        holds, where y is defined as */

/*                  T        T      T      T      T       T */
/*           y = B12  - ( E12  * U11  + E22  * U12  ) * M2 . */

/*        If B22_tilde is the square triangular matrix arising from the */
/*        (rectangular) QR-factorization */

/*                       ( B22_tilde )     ( B22  ) */
/*           Q_B_tilde * (           )  =  (      ), */
/*                       (     O     )     ( y**T ) */

/*        where Q_B_tilde is an orthogonal matrix of order N, then */

/*                T              T                T */
/*           - B22  * B22 - y * y   =  - B22_tilde  * B22_tilde. */

/*        Replacing the right hand side in (8) by the term */
/*        - B22_tilde**T * B22_tilde leads to a reduced generalized */
/*        Lyapunov equation of lower dimension compared to (5). */

/*     The recursive application of the steps I to III yields the */
/*     solution U_s of the equation (5). */

/*     It remains to compute the solution matrix U of the original */
/*     problem (1) or (2) from the matrix U_s. To this end we transform */
/*     the solution back (with respect to the transformation that led */
/*     from (1) to (5) (from (2) to (6)) and apply the QR-factorization */
/*     (RQ-factorization). The upper triangular solution matrix U is */
/*     obtained by */

/*        Q_U * U  =  U_s * Q**T     (if TRANS = 'N') */

/*     or */

/*        U * Q_U  =  Z * U_s        (if TRANS = 'T') */

/*     where Q_U is an N-by-N orthogonal matrix. Again, the orthogonal */
/*     matrix Q_U need not be accumulated. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-323, 1982. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The number of flops required by the routine is given by the */
/*     following table. Note that we count a single floating point */
/*     arithmetic operation as one flop. */

/*                 |           FACT = 'F'                  FACT = 'N' */
/*        ---------+-------------------------------------------------- */
/*         M <= N  |     (13*N**3+6*M*N**2         (211*N**3+6*M*N**2 */
/*                 |   +6*M**2*N-2*M**3)/3        +6*M**2*N-2*M**3)/3 */
/*                 | */
/*          M > N  | (11*N**3+12*M*N**2)/3     (209*N**3+12*M*N**2)/3 */

/*     FURTHER COMMENTS */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if DICO = 'D' and the pencil A - lambda * E has a pair of almost */
/*     reciprocal eigenvalues, or DICO = 'C' and the pencil has an almost */
/*     degenerate pair of eigenvalues, then the Lyapunov equation will be */
/*     ill-conditioned. Perturbed values were used to solve the equation. */
/*     A condition estimate can be obtained from the routine SG03AD. */
/*     When setting the error indicator INFO, the routine does not test */
/*     for near instability in the equation but only for exact */
/*     instability. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     May 1999 (V. Sima). */
/*     March 2002 (A. Varga). */
/*     Feb. 2004 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    --dwork;

    /* Function Body */
    isdisc = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    isfact = lsame_(fact, "F", (ftnlen)1, (ftnlen)1);
    istran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

/*     Compute minimal workspace. */

    if (isfact) {
/* Computing MAX */
	i__1 = 1, i__2 = *n << 1, i__1 = max(i__1,i__2), i__2 = *n * 6 - 6;
	minwrk = max(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = *n << 2, i__1 = max(i__1,i__2), i__2 = *n * 6 - 6;
	minwrk = max(i__1,i__2);
    }

/*     Check the scalar input parameters. */

    if (! (isdisc || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (isfact || lsame_(fact, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (istran || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*lde < max(1,*n)) {
	*info = -9;
    } else if (*ldq < max(1,*n)) {
	*info = -11;
    } else if (*ldz < max(1,*n)) {
	*info = -13;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (istran && *ldb < max(1,*n) || ! istran && *ldb < max(i__1,*n)) {
	    *info = -15;
	} else if (*ldwork < minwrk) {
	    *info = -21;
	} else {
	    *info = 0;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SG03BD", &i__1, (ftnlen)6);
	return 0;
    }

    *scale = 1.;

/*     Quick return if possible. */

    minmn = min(*m,*n);
    if (minmn == 0) {
	if (*n > 0) {
	    dlaset_("Full", n, n, &c_b10, &c_b10, &b[b_offset], ldb, (ftnlen)
		    4);
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
		*info = 2;
		return 0;
	    }
/* L20: */
	}
    }

    if (! isfact) {

/*        Reduce the pencil A - lambda * E to generalized Schur form. */

/*           A := Q**T * A * Z   (upper quasitriangular) */
/*           E := Q**T * E * Z   (upper triangular) */

/*        ( Workspace: >= MAX(1,4*N) ) */

	dgegs_("Vectors", "Vectors", n, &a[a_offset], lda, &e[e_offset], lde, 
		&alphar[1], &alphai[1], &beta[1], &q[q_offset], ldq, &z__[
		z_offset], ldz, &dwork[1], ldwork, &info1, (ftnlen)7, (ftnlen)
		7);
	if (info1 != 0) {
	    *info = 4;
	    return 0;
	}
	optwrk = (integer) dwork[1];
    } else {
	optwrk = minwrk;
    }

    if (isfact) {

/*        If the matrix pencil A - lambda * E has been in generalized */
/*        Schur form on entry, compute its eigenvalues. */

	safmin = dlamch_("Safe minimum", (ftnlen)12);
	e1[1] = 0.;
	i__ = 1;
/*        WHILE ( I .LE. N ) DO */
L30:
	if (i__ <= *n) {
/* Computing MIN */
	    i__1 = i__ + 1;
	    if (i__ == *n || a[min(i__1,*n) + i__ * a_dim1] == 0.) {
		alphar[i__] = a[i__ + i__ * a_dim1];
		alphai[i__] = 0.;
		beta[i__] = e[i__ + i__ * e_dim1];
		++i__;
	    } else {
		e1[0] = e[i__ + i__ * e_dim1];
		e1[2] = e[i__ + (i__ + 1) * e_dim1];
		e1[3] = e[i__ + 1 + (i__ + 1) * e_dim1];
		dlag2_(&a[i__ + i__ * a_dim1], lda, e1, &c__2, &safmin, &s1, &
			s2, &wr1, &wr2, &wi);
		if (wi == 0.) {
		    *info = 3;
		}
		alphar[i__] = wr1;
		alphai[i__] = wi;
		beta[i__] = s1;
		alphar[i__ + 1] = wr2;
		alphai[i__ + 1] = -wi;
		beta[i__ + 1] = s2;
		i__ += 2;
	    }
	    goto L30;
	}
/*        END WHILE 30 */
	if (*info != 0) {
	    return 0;
	}
    }

/*     Check on the stability of the matrix pencil A - lambda * E. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (isdisc) {
	    if (dlapy2_(&alphar[i__], &alphai[i__]) >= (d__1 = beta[i__], abs(
		    d__1))) {
		*info = 6;
		return 0;
	    }
	} else {
	    if (alphar[i__] == 0. || beta[i__] == 0. || d_sign(&c_b19, &
		    alphar[i__]) * d_sign(&c_b19, &beta[i__]) >= 0.) {
		*info = 5;
		return 0;
	    }
	}
/* L40: */
    }

/*     Transformation of the right hand side. */

/*        B := B * Z  or  B := Q**T * B */

/*     Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*     ( Workspace: max(1,N) ) */

    if (! istran) {
	if (*ldwork >= *n * *m) {
	    dgemm_("NoTranspose", "NoTranspose", m, n, n, &c_b19, &b[b_offset]
		    , ldb, &z__[z_offset], ldz, &c_b10, &dwork[1], m, (ftnlen)
		    11, (ftnlen)11);
	    dlacpy_("All", m, n, &dwork[1], m, &b[b_offset], ldb, (ftnlen)3);
	} else {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &b[i__ + b_dim1], ldb, &dwork[1], &c__1);
		dgemv_("Transpose", n, n, &c_b19, &z__[z_offset], ldz, &dwork[
			1], &c__1, &c_b10, &b[i__ + b_dim1], ldb, (ftnlen)9);
/* L60: */
	    }
	}
	if (*m < *n) {
	    i__1 = *n - *m;
	    dlaset_("All", &i__1, n, &c_b10, &c_b10, &b[*m + 1 + b_dim1], ldb,
		     (ftnlen)3);
	}
    } else {
	if (*ldwork >= *n * *m) {
	    dlacpy_("All", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)3);
	    dgemm_("Transpose", "NoTranspose", n, m, n, &c_b19, &q[q_offset], 
		    ldq, &dwork[1], n, &c_b10, &b[b_offset], ldb, (ftnlen)9, (
		    ftnlen)11);
	} else {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[1], &c__1);
		dgemv_("Transpose", n, n, &c_b19, &q[q_offset], ldq, &dwork[1]
			, &c__1, &c_b10, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)
			9);
/* L80: */
	    }
	}
	if (*m < *n) {
	    i__1 = *n - *m;
	    dlaset_("All", n, &i__1, &c_b10, &c_b10, &b[(*m + 1) * b_dim1 + 1]
		    , ldb, (ftnlen)3);
	}
    }
/* Computing MAX */
    i__1 = optwrk, i__2 = *n * *m;
    optwrk = max(i__1,i__2);

/*     Overwrite B with the triangular matrix of its QR-factorization */
/*     or its RQ-factorization. */
/*     (The entries on the main diagonal are non-negative.) */

/*     ( Workspace: >= max(1,2*N) ) */

    if (! istran) {
	if (*m >= 2) {
	    i__1 = *ldwork - *n;
	    dgeqrf_(m, n, &b[b_offset], ldb, &dwork[1], &dwork[*n + 1], &i__1,
		     &info1);
	    i__1 = max(*m,*n) - 1;
	    i__2 = min(*m,*n);
	    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[b_dim1 + 2], 
		    ldb, (ftnlen)5);
	}
	i__1 = minmn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (b[i__ + i__ * b_dim1] < 0.) {
		i__2 = *n + 1 - i__;
		dscal_(&i__2, &c_b55, &b[i__ + i__ * b_dim1], ldb);
	    }
/* L100: */
	}
    } else {
	if (*m >= 2) {
	    i__1 = *ldwork - *n;
	    dgerqf_(n, m, &b[b_offset], ldb, &dwork[1], &dwork[*n + 1], &i__1,
		     &info1);
	    if (*n >= *m) {
		i__1 = *m - 1;
		i__2 = *m - 1;
		dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[*n - *m + 2 
			+ b_dim1], ldb, (ftnlen)5);
		if (*n > *m) {
		    for (i__ = *m; i__ >= 1; --i__) {
			dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &b[(i__ + *n - 
				*m) * b_dim1 + 1], &c__1);
/* L120: */
		    }
		    i__1 = *n - *m;
		    dlaset_("All", n, &i__1, &c_b10, &c_b10, &b[b_dim1 + 1], 
			    ldb, (ftnlen)3);
		}
	    } else {
		if (*n > 1) {
		    i__1 = *n - 1;
		    i__2 = *n - 1;
		    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[(*m - *
			    n + 1) * b_dim1 + 2], ldb, (ftnlen)5);
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dcopy_(n, &b[(*m - *n + i__) * b_dim1 + 1], &c__1, &b[i__ 
			    * b_dim1 + 1], &c__1);
/* L140: */
		}
		i__1 = *m - *n;
		dlaset_("All", n, &i__1, &c_b10, &c_b10, &b[(*n + 1) * b_dim1 
			+ 1], ldb, (ftnlen)3);
	    }
	} else {
	    if (*n != 1) {
		dcopy_(n, &b[b_dim1 + 1], &c__1, &b[*n * b_dim1 + 1], &c__1);
		dlaset_("All", n, &c__1, &c_b10, &c_b10, &b[b_dim1 + 1], ldb, 
			(ftnlen)3);
	    }
	}
	i__1 = *n;
	for (i__ = *n - minmn + 1; i__ <= i__1; ++i__) {
	    if (b[i__ + i__ * b_dim1] < 0.) {
		dscal_(&i__, &c_b55, &b[i__ * b_dim1 + 1], &c__1);
	    }
/* L160: */
	}
    }
/* Computing MAX */
    i__1 = optwrk, i__2 = (integer) dwork[*n + 1] + *n;
    optwrk = max(i__1,i__2);

/*     Solve the reduced generalized Lyapunov equation. */

/*     ( Workspace: 6*N-6 ) */

    if (isdisc) {
	sg03bu_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &b[b_offset], 
		ldb, scale, &dwork[1], &info1, (ftnlen)1);
	if (info1 != 0) {
	    if (info1 == 1) {
		*info = 1;
	    }
	    if (info1 == 2) {
		*info = 3;
	    }
	    if (info1 == 3) {
		*info = 6;
	    }
	    if (info1 == 4) {
		*info = 7;
	    }
	    if (*info != 1) {
		return 0;
	    }
	}
    } else {
	sg03bv_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &b[b_offset], 
		ldb, scale, &dwork[1], &info1, (ftnlen)1);
	if (info1 != 0) {
	    if (info1 == 1) {
		*info = 1;
	    }
	    if (info1 >= 2) {
		*info = 3;
	    }
	    if (info1 == 3) {
		*info = 5;
	    }
	    if (*info != 1) {
		return 0;
	    }
	}
    }

/*     Transform the solution matrix back. */

/*        U := U * Q**T   or   U := Z * U */

/*     Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*     ( Workspace: max(1,N) ) */

    if (! istran) {
	if (*ldwork >= *n * *n) {
	    dlacpy_("All", n, n, &q[q_offset], ldq, &dwork[1], n, (ftnlen)3);
	    dtrmm_("Right", "Upper", "Transpose", "NonUnit", n, n, &c_b19, &b[
		    b_offset], ldb, &dwork[1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)7);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &dwork[*n * (i__ - 1) + 1], &c__1, &b[i__ + b_dim1],
			 ldb);
/* L170: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - i__ + 1;
		dcopy_(&i__2, &b[i__ + i__ * b_dim1], ldb, &dwork[1], &c__1);
		i__2 = *n - i__ + 1;
		dgemv_("NoTranspose", n, &i__2, &c_b19, &q[i__ * q_dim1 + 1], 
			ldq, &dwork[1], &c__1, &c_b10, &b[i__ + b_dim1], ldb, 
			(ftnlen)11);
/* L180: */
	    }
	}
    } else {
	if (*ldwork >= *n * *n) {
	    dlacpy_("All", n, n, &z__[z_offset], ldz, &dwork[1], n, (ftnlen)3)
		    ;
	    dtrmm_("Right", "Upper", "NoTranspose", "NonUnit", n, n, &c_b19, &
		    b[b_offset], ldb, &dwork[1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)11, (ftnlen)7);
	    dlacpy_("All", n, n, &dwork[1], n, &b[b_offset], ldb, (ftnlen)3);
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(&i__, &b[i__ * b_dim1 + 1], &c__1, &dwork[1], &c__1);
		dgemv_("NoTranspose", n, &i__, &c_b19, &z__[z_offset], ldz, &
			dwork[1], &c__1, &c_b10, &b[i__ * b_dim1 + 1], &c__1, 
			(ftnlen)11);
/* L200: */
	    }
	}
    }
/* Computing MAX */
    i__1 = optwrk, i__2 = *n * *n;
    optwrk = max(i__1,i__2);

/*     Overwrite U with the triangular matrix of its QR-factorization */
/*     or its RQ-factorization. */
/*     (The entries on the main diagonal are non-negative.) */

/*     ( Workspace: >= max(1,2*N) ) */

    if (! istran) {
	i__1 = *ldwork - *n;
	dgeqrf_(n, n, &b[b_offset], ldb, &dwork[1], &dwork[*n + 1], &i__1, &
		info1);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[b_dim1 + 2], 
		    ldb, (ftnlen)5);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (b[i__ + i__ * b_dim1] < 0.) {
		i__2 = *n + 1 - i__;
		dscal_(&i__2, &c_b55, &b[i__ + i__ * b_dim1], ldb);
	    }
/* L220: */
	}
    } else {
	i__1 = *ldwork - *n;
	dgerqf_(n, n, &b[b_offset], ldb, &dwork[1], &dwork[*n + 1], &i__1, &
		info1);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[b_dim1 + 2], 
		    ldb, (ftnlen)5);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (b[i__ + i__ * b_dim1] < 0.) {
		dscal_(&i__, &c_b55, &b[i__ * b_dim1 + 1], &c__1);
	    }
/* L240: */
	}
    }
/* Computing MAX */
    i__1 = optwrk, i__2 = (integer) dwork[*n + 1] + *n;
    optwrk = max(i__1,i__2);

    dwork[1] = (doublereal) max(optwrk,minwrk);
    return 0;
/* *** Last line of SG03BD *** */
} /* sg03bd_ */

