/* SB03OD.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;
static doublereal c_b22 = 1.;

/* Subroutine */ int sb03od_(char *dico, char *fact, char *trans, integer *n, 
	integer *m, doublereal *a, integer *lda, doublereal *q, integer *ldq, 
	doublereal *b, integer *ldb, doublereal *scale, doublereal *wr, 
	doublereal *wi, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen fact_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, i__1, i__2;
    logical L__1;

    /* Local variables */
    static integer i__, j, k, l, ne;
    static doublereal emax;
    static integer sdim, itau;
    static logical cont;
    static doublereal temp;
    static integer ifail;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sb03ou_(logical *,
	     logical *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer minmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical bwork[1];
    static integer jwork;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgerqf_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, integer *);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static integer inform__;
    static logical ltrans;
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

/*     To solve for X = op(U)'*op(U) either the stable non-negative */
/*     definite continuous-time Lyapunov equation */
/*                                   2 */
/*        op(A)'*X + X*op(A) = -scale *op(B)'*op(B)                   (1) */

/*     or the convergent non-negative definite discrete-time Lyapunov */
/*     equation */
/*                                   2 */
/*        op(A)'*X*op(A) - X = -scale *op(B)'*op(B)                   (2) */

/*     where op(K) = K or K' (i.e., the transpose of the matrix K), A is */
/*     an N-by-N matrix, op(B) is an M-by-N matrix, U is an upper */
/*     triangular matrix containing the Cholesky factor of the solution */
/*     matrix X, X = op(U)'*op(U), and scale is an output scale factor, */
/*     set less than or equal to 1 to avoid overflow in X. If matrix B */
/*     has full rank then the solution matrix X will be positive-definite */
/*     and hence the Cholesky factor U will be nonsingular, but if B is */
/*     rank deficient then X may be only positive semi-definite and U */
/*     will be singular. */

/*     In the case of equation (1) the matrix A must be stable (that */
/*     is, all the eigenvalues of A must have negative real parts), */
/*     and for equation (2) the matrix A must be convergent (that is, */
/*     all the eigenvalues of A must lie inside the unit circle). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of Lyapunov equation to be solved as */
/*             follows: */
/*             = 'C':  Equation (1), continuous-time case; */
/*             = 'D':  Equation (2), discrete-time case. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix A is supplied on entry, as follows: */
/*             = 'F':  On entry, A and Q contain the factors from the */
/*                     real Schur factorization of the matrix A; */
/*             = 'N':  The Schur factorization of A will be computed */
/*                     and the factors will be stored in A and Q. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = 'N':  op(K) = K    (No transpose); */
/*             = 'T':  op(K) = K**T (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and the number of columns in */
/*             matrix op(B).  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of rows in matrix op(B).  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. If FACT = 'F', then A contains */
/*             an upper quasi-triangular matrix S in Schur canonical */
/*             form; the elements below the upper Hessenberg part of the */
/*             array A are not referenced. */
/*             On exit, the leading N-by-N upper Hessenberg part of this */
/*             array contains the upper quasi-triangular matrix S in */
/*             Schur canonical form from the Shur factorization of A. */
/*             The contents of array A is not modified if FACT = 'F'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     Q       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDQ,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N part of */
/*             this array must contain the orthogonal matrix Q of the */
/*             Schur factorization of A. */
/*             Otherwise, Q need not be set on entry. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal matrix Q of the Schur factorization of A. */
/*             The contents of array Q is not modified if FACT = 'F'. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             if TRANS = 'N', and dimension (LDB,max(M,N)), if */
/*             TRANS = 'T'. */
/*             On entry, if TRANS = 'N', the leading M-by-N part of this */
/*             array must contain the coefficient matrix B of the */
/*             equation. */
/*             On entry, if TRANS = 'T', the leading N-by-M part of this */
/*             array must contain the coefficient matrix B of the */
/*             equation. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper triangular Cholesky factor U of the solution */
/*             matrix X of the problem, X = op(U)'*op(U). */
/*             If M = 0 and N > 0, then U is set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N,M), if TRANS = 'N'; */
/*             LDB >= MAX(1,N),   if TRANS = 'T'. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             If FACT = 'N', and INFO >= 0 and INFO <= 2, WR and WI */
/*             contain the real and imaginary parts, respectively, of */
/*             the eigenvalues of A. */
/*             If FACT = 'F', WR and WI are not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, or INFO = 1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If M > 0, LDWORK >= MAX(1,4*N + MIN(M,N)); */
/*             If M = 0, LDWORK >= 1. */
/*             For optimum performance LDWORK should sometimes be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the Lyapunov equation is (nearly) singular */
/*                   (warning indicator); */
/*                   if DICO = 'C' this means that while the matrix A */
/*                   (or the factor S) has computed eigenvalues with */
/*                   negative real parts, it is only just stable in the */
/*                   sense that small perturbations in A can make one or */
/*                   more of the eigenvalues have a non-negative real */
/*                   part; */
/*                   if DICO = 'D' this means that while the matrix A */
/*                   (or the factor S) has computed eigenvalues inside */
/*                   the unit circle, it is nevertheless only just */
/*                   convergent, in the sense that small perturbations */
/*                   in A can make one or more of the eigenvalues lie */
/*                   outside the unit circle; */
/*                   perturbed values were used to solve the equation; */
/*             = 2:  if FACT = 'N' and DICO = 'C', but the matrix A is */
/*                   not stable (that is, one or more of the eigenvalues */
/*                   of A has a non-negative real part), or DICO = 'D', */
/*                   but the matrix A is not convergent (that is, one or */
/*                   more of the eigenvalues of A lies outside the unit */
/*                   circle); however, A will still have been factored */
/*                   and the eigenvalues of A returned in WR and WI. */
/*             = 3:  if FACT = 'F' and DICO = 'C', but the Schur factor S */
/*                   supplied in the array A is not stable (that is, one */
/*                   or more of the eigenvalues of S has a non-negative */
/*                   real part), or DICO = 'D', but the Schur factor S */
/*                   supplied in the array A is not convergent (that is, */
/*                   one or more of the eigenvalues of S lies outside the */
/*                   unit circle); */
/*             = 4:  if FACT = 'F' and the Schur factor S supplied in */
/*                   the array A has two or more consecutive non-zero */
/*                   elements on the first sub-diagonal, so that there is */
/*                   a block larger than 2-by-2 on the diagonal; */
/*             = 5:  if FACT = 'F' and the Schur factor S supplied in */
/*                   the array A has a 2-by-2 diagonal block with real */
/*                   eigenvalues instead of a complex conjugate pair; */
/*             = 6:  if FACT = 'N' and the LAPACK Library routine DGEES */
/*                   has failed to converge. This failure is not likely */
/*                   to occur. The matrix B will be unaltered but A will */
/*                   be destroyed. */

/*     METHOD */

/*     The method used by the routine is based on the Bartels and Stewart */
/*     method [1], except that it finds the upper triangular matrix U */
/*     directly without first finding X and without the need to form the */
/*     normal matrix op(B)'*op(B). */

/*     The Schur factorization of a square matrix A is given by */

/*        A = QSQ', */

/*     where Q is orthogonal and S is an N-by-N block upper triangular */
/*     matrix with 1-by-1 and 2-by-2 blocks on its diagonal (which */
/*     correspond to the eigenvalues of A). If A has already been */
/*     factored prior to calling the routine however, then the factors */
/*     Q and S may be supplied and the initial factorization omitted. */

/*     If TRANS = 'N', the matrix B is factored as (QR factorization) */
/*            _   _                   _   _  _ */
/*        B = P ( R ),  M >= N,   B = P ( R  Z ),  M < N, */
/*              ( 0 ) */
/*           _                                    _ */
/*     where P is an M-by-M orthogonal matrix and R is a square upper */
/*                                         _   _      _     _  _ */
/*     triangular matrix. Then, the matrix B = RQ, or B = ( R  Z )Q (if */
/*     M < N) is factored as */
/*        _                       _ */
/*        B = P ( R ),  M >= N,   B = P ( R  Z ),  M < N. */

/*     If TRANS = 'T', the matrix B is factored as (RQ factorization) */
/*                                         _ */
/*                 _   _                 ( Z ) _ */
/*        B = ( 0  R ) P,  M >= N,   B = ( _ ) P,  M < N, */
/*                                       ( R ) */
/*           _                                    _ */
/*     where P is an M-by-M orthogonal matrix and R is a square upper */
/*                                         _     _     _       _   _ */
/*     triangular matrix. Then, the matrix B = Q'R, or B = Q'( Z'  R' )' */
/*     (if M < N) is factored as */
/*        _                       _ */
/*        B = ( R ) P,  M >= N,   B = ( Z ) P,  M < N. */
/*                                    ( R ) */

/*     These factorizations are utilised to either transform the */
/*     continuous-time Lyapunov equation to the canonical form */
/*                                                        2 */
/*       op(S)'*op(V)'*op(V) + op(V)'*op(V)*op(S) = -scale *op(F)'*op(F), */

/*     or the discrete-time Lyapunov equation to the canonical form */
/*                                                        2 */
/*       op(S)'*op(V)'*op(V)*op(S) - op(V)'*op(V) = -scale *op(F)'*op(F), */

/*     where V and F are upper triangular, and */

/*        F = R,  M >= N,   F = ( R  Z ),  M < N, if TRANS = 'N'; */
/*                              ( 0  0 ) */

/*        F = R,  M >= N,   F = ( 0  Z ),  M < N, if TRANS = 'T'. */
/*                              ( 0  R ) */

/*     The transformed equation is then solved for V, from which U is */
/*     obtained via the QR factorization of V*Q', if TRANS = 'N', or */
/*     via the RQ factorization of Q*V, if TRANS = 'T'. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W. */
/*         Solution of the matrix equation  A'X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if A is only just stable (or convergent) then the Lyapunov */
/*     equation will be ill-conditioned.  A symptom of ill-conditioning */
/*     is "large" elements in U relative to those of A and B, or a */
/*     "small" value for scale. A condition estimate can be computed */
/*     using SLICOT Library routine SB03MD. */

/*     SB03OD routine can be also used for solving "unstable" Lyapunov */
/*     equations, i.e., when matrix A has all eigenvalues with positive */
/*     real parts, if DICO = 'C', or with moduli greater than one, */
/*     if DICO = 'D'. Specifically, one may solve for X = op(U)'*op(U) */
/*     either the continuous-time Lyapunov equation */
/*                                  2 */
/*        op(A)'*X + X*op(A) = scale *op(B)'*op(B),                   (3) */

/*     or the discrete-time Lyapunov equation */
/*                                  2 */
/*        op(A)'*X*op(A) - X = scale *op(B)'*op(B),                   (4) */

/*     provided, for equation (3), the given matrix A is replaced by -A, */
/*     or, for equation (4), the given matrices A and B are replaced by */
/*     inv(A) and B*inv(A), if TRANS = 'N' (or inv(A)*B, if TRANS = 'T'), */
/*     respectively. Although the inversion generally can rise numerical */
/*     problems, in case of equation (4) it is expected that the matrix A */
/*     is enough well-conditioned, having only eigenvalues with moduli */
/*     greater than 1. However, if A is ill-conditioned, it could be */
/*     preferable to use the more general SLICOT Lyapunov solver SB03MD. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB03CD by Sven Hammarling, */
/*     NAG Ltd, United Kingdom. */

/*     REVISIONS */

/*     Dec. 1997, April 1998, May 1998, May 1999, Oct. 2001 (V. Sima). */
/*     March 2002 (A. Varga). */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

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

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    cont = lsame_(dico, "C", (ftnlen)1, (ftnlen)1);
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    minmn = min(*m,*n);

    *info = 0;
    if (! cont && ! lsame_(dico, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldq < max(1,*n)) {
	*info = -9;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*n);
	if (*ldb < max(1,*n) || *ldb < max(i__1,*m) && ! ltrans) {
	    *info = -11;
	} else if (*ldwork < 1 || *m > 0 && *ldwork < (*n << 2) + minmn) {
	    *info = -16;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB03OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (minmn == 0) {
	if (*m == 0) {
	    dlaset_("Full", n, n, &c_b10, &c_b10, &b[b_offset], ldb, (ftnlen)
		    4);
	}
	*scale = 1.;
	dwork[1] = 1.;
	return 0;
    }

/*     Start the solution. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    if (nofact) {

/*        Find the Schur factorization of A,   A = Q*S*Q'. */
/*        Workspace:  need   3*N; */
/*                    prefer larger. */

	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, 
		&sdim, &wr[1], &wi[1], &q[q_offset], ldq, &dwork[1], ldwork, 
		bwork, &inform__, (ftnlen)7, (ftnlen)11);
	if (inform__ != 0) {
	    *info = 6;
	    return 0;
	}
	wrkopt = (integer) dwork[1];

/*        Check the eigenvalues for stability. */

	if (cont) {
	    emax = wr[1];

	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		if (wr[j] > emax) {
		    emax = wr[j];
		}
/* L20: */
	    }

	} else {
	    emax = dlapy2_(&wr[1], &wi[1]);

	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		temp = dlapy2_(&wr[j], &wi[j]);
		if (temp > emax) {
		    emax = temp;
		}
/* L40: */
	    }

	}

	if (cont && emax >= 0. || ! cont && emax >= 1.) {
	    *info = 2;
	    return 0;
	}
    } else {
	wrkopt = 0;
    }

/*     Perform the QR or RQ factorization of B, */
/*            _   _           _   _  _ */
/*        B = P ( R ), or B = P ( R  Z ), if TRANS = 'N', or */
/*              ( 0 ) */
/*                                 _ */
/*                 _   _         ( Z ) _ */
/*        B = ( 0  R ) P, or B = ( _ ) P, if TRANS = 'T'. */
/*                               ( R ) */
/*     Workspace:  need   MIN(M,N) + N; */
/*                 prefer MIN(M,N) + N*NB. */

    itau = 1;
    jwork = itau + minmn;
    if (ltrans) {
	i__1 = *ldwork - jwork + 1;
	dgerqf_(n, m, &b[b_offset], ldb, &dwork[itau], &dwork[jwork], &i__1, &
		ifail);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1, i__1 = max(
		i__1,i__2), i__2 = minmn * *n;
	wrkopt = max(i__1,i__2);
	jwork = itau;

/*        Form in B */
/*        _      _              _         _   _                    _ */
/*        B := Q'R,   m >= n,   B := Q'*( Z'  R' )',   m < n, with B an */
/*        n-by-min(m,n) matrix. */
/*        Use a BLAS 3 operation if enough workspace, and BLAS 2, */
/*                   _ */
/*        otherwise: B is formed column by column. */

	if (*ldwork >= jwork + minmn * *n - 1) {
	    k = jwork;

	    i__1 = minmn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &q[*n - minmn + i__ + q_dim1], ldq, &dwork[k], &
			c__1);
		k += *n;
/* L60: */
	    }

	    dtrmm_("Right", "Upper", "No transpose", "Non-unit", n, &minmn, &
		    c_b22, &b[*n - minmn + 1 + (*m - minmn + 1) * b_dim1], 
		    ldb, &dwork[jwork], n, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
		    ftnlen)8);
	    if (*m < *n) {
		i__1 = *n - *m;
		dgemm_("Transpose", "No transpose", n, m, &i__1, &c_b22, &q[
			q_offset], ldq, &b[b_offset], ldb, &c_b22, &dwork[
			jwork], n, (ftnlen)9, (ftnlen)12);
	    }
	    dlacpy_("Full", n, &minmn, &dwork[jwork], n, &b[b_offset], ldb, (
		    ftnlen)4);
	} else {
	    ne = *n - minmn;

	    i__1 = minmn;
	    for (j = 1; j <= i__1; ++j) {
		++ne;
		dcopy_(&ne, &b[(*m - minmn + j) * b_dim1 + 1], &c__1, &dwork[
			jwork], &c__1);
		dgemv_("Transpose", &ne, n, &c_b22, &q[q_offset], ldq, &dwork[
			jwork], &c__1, &c_b10, &b[j * b_dim1 + 1], &c__1, (
			ftnlen)9);
/* L80: */
	    }

	}
    } else {
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(m, n, &b[b_offset], ldb, &dwork[itau], &dwork[jwork], &i__1, &
		ifail);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1, i__1 = max(
		i__1,i__2), i__2 = minmn * *n;
	wrkopt = max(i__1,i__2);
	jwork = itau;

/*        Form in B */
/*        _    _               _      _  _                    _ */
/*        B := RQ,   m >= n,   B := ( R  Z )*Q,   m < n, with B an */
/*        min(m,n)-by-n matrix. */
/*        Use a BLAS 3 operation if enough workspace, and BLAS 2, */
/*                   _ */
/*        otherwise: B is formed row by row. */

	if (*ldwork >= jwork + minmn * *n - 1) {
	    dlacpy_("Full", &minmn, n, &q[q_offset], ldq, &dwork[jwork], &
		    minmn, (ftnlen)4);
	    dtrmm_("Left", "Upper", "No transpose", "Non-unit", &minmn, n, &
		    c_b22, &b[b_offset], ldb, &dwork[jwork], &minmn, (ftnlen)
		    4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
	    if (*m < *n) {
		i__1 = *n - *m;
		dgemm_("No transpose", "No transpose", m, n, &i__1, &c_b22, &
			b[(*m + 1) * b_dim1 + 1], ldb, &q[*m + 1 + q_dim1], 
			ldq, &c_b22, &dwork[jwork], &minmn, (ftnlen)12, (
			ftnlen)12);
	    }
	    dlacpy_("Full", &minmn, n, &dwork[jwork], &minmn, &b[b_offset], 
		    ldb, (ftnlen)4);
	} else {
/* Computing MAX */
	    i__1 = 0, i__2 = *n - *m;
	    ne = minmn + max(i__1,i__2);

	    i__1 = minmn;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(&ne, &b[j + j * b_dim1], ldb, &dwork[jwork], &c__1);
		dgemv_("Transpose", &ne, n, &c_b22, &q[j + q_dim1], ldq, &
			dwork[jwork], &c__1, &c_b10, &b[j + b_dim1], ldb, (
			ftnlen)9);
		--ne;
/* L100: */
	    }

	}
    }
    jwork = itau + minmn;

/*     Solve for U the transformed Lyapunov equation */
/*                                                      2    _      _ */
/*     op(S)'*op(U)'*op(U) + op(U)'*op(U)*op(S) = -scale *op(B)'*op(B), */

/*     or */
/*                                                      2    _      _ */
/*     op(S)'*op(U)'*op(U)*op(S) - op(U)'*op(U) = -scale *op(B)'*op(B) */

/*     Workspace:  need   MIN(M,N) + 4*N; */
/*                 prefer larger. */

    L__1 = ! cont;
    i__1 = *ldwork - jwork + 1;
    sb03ou_(&L__1, &ltrans, n, &minmn, &a[a_offset], lda, &b[b_offset], ldb, &
	    dwork[itau], &b[b_offset], ldb, scale, &dwork[jwork], &i__1, info)
	    ;
    if (*info > 1) {
	++(*info);
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);
    jwork = itau;

/*     Form   U :=  U*Q' or U := Q*U in the array B. */
/*     Use a BLAS 3 operation if enough workspace, and BLAS 2, otherwise. */
/*     Workspace:  need   N; */
/*                 prefer N*N; */

    if (*ldwork >= jwork + *n * *n - 1) {
	if (ltrans) {
	    dlacpy_("Full", n, n, &q[q_offset], ldq, &dwork[jwork], n, (
		    ftnlen)4);
	    dtrmm_("Right", "Upper", "No transpose", "Non-unit", n, n, &c_b22,
		     &b[b_offset], ldb, &dwork[jwork], n, (ftnlen)5, (ftnlen)
		    5, (ftnlen)12, (ftnlen)8);
	} else {
	    k = jwork;

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &q[i__ * q_dim1 + 1], &c__1, &dwork[k], n);
		++k;
/* L120: */
	    }

	    dtrmm_("Left", "Upper", "No transpose", "Non-unit", n, n, &c_b22, 
		    &b[b_offset], ldb, &dwork[jwork], n, (ftnlen)4, (ftnlen)5,
		     (ftnlen)12, (ftnlen)8);
	}
	dlacpy_("Full", n, n, &dwork[jwork], n, &b[b_offset], ldb, (ftnlen)4);
/* Computing MAX */
	i__1 = wrkopt, i__2 = jwork + *n * *n - 1;
	wrkopt = max(i__1,i__2);
    } else {
	if (ltrans) {

/*           U is formed column by column ( U := Q*U ). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(&i__, &b[i__ * b_dim1 + 1], &c__1, &dwork[jwork], &
			c__1);
		dgemv_("No transpose", n, &i__, &c_b22, &q[q_offset], ldq, &
			dwork[jwork], &c__1, &c_b10, &b[i__ * b_dim1 + 1], &
			c__1, (ftnlen)12);
/* L140: */
	    }
	} else {

/*           U is formed row by row ( U' := Q*U' ). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - i__ + 1;
		dcopy_(&i__2, &b[i__ + i__ * b_dim1], ldb, &dwork[jwork], &
			c__1);
		i__2 = *n - i__ + 1;
		dgemv_("No transpose", n, &i__2, &c_b22, &q[i__ * q_dim1 + 1],
			 ldq, &dwork[jwork], &c__1, &c_b10, &b[i__ + b_dim1], 
			ldb, (ftnlen)12);
/* L160: */
	    }
	}
    }

/*     Lastly find the QR or RQ factorization of U, overwriting on B, */
/*     to give the required Cholesky factor. */
/*     Workspace:  need   2*N; */
/*                 prefer N + N*NB; */

    jwork = itau + *n;
    if (ltrans) {
	i__1 = *ldwork - jwork + 1;
	dgerqf_(n, n, &b[b_offset], ldb, &dwork[itau], &dwork[jwork], &i__1, &
		ifail);
    } else {
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(n, n, &b[b_offset], ldb, &dwork[itau], &dwork[jwork], &i__1, &
		ifail);
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Make the diagonal elements of U non-negative. */

    if (ltrans) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (b[j + j * b_dim1] < 0.) {

		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] = -b[i__ + j * b_dim1];
/* L180: */
		}

	    }
/* L200: */
	}

    } else {
	k = jwork;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dwork[k] = b[j + j * b_dim1];
	    l = jwork;

	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (dwork[l] < 0.) {
		    b[i__ + j * b_dim1] = -b[i__ + j * b_dim1];
		}
		++l;
/* L220: */
	    }

	    ++k;
/* L240: */
	}
    }

    if (*n > 1) {
	i__1 = *n - 1;
	i__2 = *n - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &b[b_dim1 + 2], ldb, (
		ftnlen)5);
    }

/*     Set the optimal workspace. */

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of SB03OD *** */
} /* sb03od_ */

