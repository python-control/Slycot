/* SG02AD.f -- translated by f2c (version 20100827).
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
static doublereal c_b30 = 1.;
static integer c__1 = 1;
static doublereal c_b71 = 0.;
static doublereal c_b107 = .5;

/* Subroutine */ int sg02ad_(char *dico, char *jobb, char *fact, char *uplo, 
	char *jobl, char *scal, char *sort, char *acc, integer *n, integer *m,
	 integer *p, doublereal *a, integer *lda, doublereal *e, integer *lde,
	 doublereal *b, integer *ldb, doublereal *q, integer *ldq, doublereal 
	*r__, integer *ldr, doublereal *l, integer *ldl, doublereal *rcondu, 
	doublereal *x, integer *ldx, doublereal *alfar, doublereal *alfai, 
	doublereal *beta, doublereal *s, integer *lds, doublereal *t, integer 
	*ldt, doublereal *u, integer *ldu, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, logical *bwork, integer *iwarn, 
	integer *info, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, 
	ftnlen uplo_len, ftnlen jobl_len, ftnlen scal_len, ftnlen sort_len, 
	ftnlen acc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, l_dim1, 
	    l_offset, q_dim1, q_offset, r_dim1, r_offset, s_dim1, s_offset, 
	    t_dim1, t_offset, u_dim1, u_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn, mp, np, iw, np1, iwb, iwc, iwf, ldw, nnm;
    static doublereal eps, u12m;
    static integer iwr, ndim;
    static doublereal asym, seps;
    static integer info1;
    static logical lfacb, lfacn;
    extern /* Subroutine */ int mb02pd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, char 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical lfacq, lfacr, ljobb;
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgges_(char *, char *, char *, L_fp, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, logical *, integer *, ftnlen,
	     ftnlen, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *), mb01sd_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, ftnlen);
    static logical lscal;
    extern /* Subroutine */ int mb02vd_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobl;
    static char equed[1];
    static logical discr;
    extern logical sb02ou_(), sb02ov_(), sb02ow_();
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), sb02oy_(char *, char *, char *, char 
	    *, char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen)
	    , daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical luplo;
    static doublereal rnorm, unorm;
    static char qtype[1], rtype[1];
    static logical lsort;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical refine;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical ljobln;
    static doublereal rcondl;
    static logical colequ;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical rowequ;
    static integer wrkopt;
    static doublereal pivotu;


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
/*                                   -1 */
/*        Q + A'XE + E'XA - (L+E'XB)R  (L+E'XB)' = 0 ,              (1) */

/*     or the discrete-time algebraic Riccati equation */
/*                                        -1 */
/*        E'XE = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q ,         (2) */

/*     where A, E, B, Q, R, and L are N-by-N, N-by-N, N-by-M, N-by-N, */
/*     M-by-M and N-by-M matrices, respectively, such that Q = C'C, */
/*     R = D'D and L = C'D; X is an N-by-N symmetric matrix. */
/*     The routine also returns the computed values of the closed-loop */
/*     spectrum of the system, i.e., the stable eigenvalues */
/*     lambda(1),...,lambda(N) of the pencil (A - BF,E), where F is */
/*     the optimal gain matrix, */
/*             -1 */
/*        F = R  (L+E'XB)' ,        for (1), */

/*     and */
/*                    -1 */
/*        F = (R+B'XB)  (L+A'XB)' , for (2). */
/*                              -1 */
/*     Optionally, matrix G = BR  B' may be given instead of B and R. */
/*     Other options include the case with Q and/or R given in a */
/*     factored form, Q = C'C, R = D'D, and with L a zero matrix. */

/*     The routine uses the method of deflating subspaces, based on */
/*     reordering the eigenvalues in a generalized Schur matrix pair. */

/*     It is assumed that E is nonsingular, but this condition is not */
/*     checked. Note that the definition (1) of the continuous-time */
/*     algebraic Riccati equation, and the formula for the corresponding */
/*     optimal gain matrix, require R to be nonsingular, but the */
/*     associated linear quadratic optimal problem could have a unique */
/*     solution even when matrix R is singular, under mild assumptions */
/*     (see METHOD). The routine SG02AD works accordingly in this case. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of Riccati equation to be solved as */
/*             follows: */
/*             = 'C':  Equation (1), continuous-time case; */
/*             = 'D':  Equation (2), discrete-time case. */

/*     JOBB    CHARACTER*1 */
/*             Specifies whether or not the matrix G is given, instead */
/*             of the matrices B and R, as follows: */
/*             = 'B':  B and R are given; */
/*             = 'G':  G is given. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the matrices Q and/or R (if */
/*             JOBB = 'B') are factored, as follows: */
/*             = 'N':  Not factored, Q and R are given; */
/*             = 'C':  C is given, and Q = C'C; */
/*             = 'D':  D is given, and R = D'D; */
/*             = 'B':  Both factors C and D are given, Q = C'C, R = D'D. */

/*     UPLO    CHARACTER*1 */
/*             If JOBB = 'G', or FACT = 'N', specifies which triangle of */
/*             the matrices G, or Q and R, is stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */
/*             JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed. */
/*             SLICOT Library routine SB02MT should be called just before */
/*             SG02AD, for obtaining the results when JOBB = 'G' and */
/*             JOBL = 'N'. */

/*     SCAL    CHARACTER*1 */
/*             If JOBB = 'B', specifies whether or not a scaling strategy */
/*             should be used to scale Q, R, and L, as follows: */
/*             = 'G':  General scaling should be used; */
/*             = 'N':  No scaling should be used. */
/*             SCAL is not used if JOBB = 'G'. */

/*     SORT    CHARACTER*1 */
/*             Specifies which eigenvalues should be obtained in the top */
/*             of the generalized Schur form, as follows: */
/*             = 'S':  Stable   eigenvalues come first; */
/*             = 'U':  Unstable eigenvalues come first. */

/*     ACC     CHARACTER*1 */
/*             Specifies whether or not iterative refinement should be */
/*             used to solve the system of algebraic equations giving */
/*             the solution matrix X, as follows: */
/*             = 'R':  Use iterative refinement; */
/*             = 'N':  Do not use iterative refinement. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices A, E, Q, and X, and the number of rows of the */
/*             matrices B and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs. If JOBB = 'B', M is the */
/*             order of the matrix R, and the number of columns of the */
/*             matrix B.  M >= 0. */
/*             M is not used if JOBB = 'G'. */

/*     P       (input) INTEGER */
/*             The number of system outputs. If FACT = 'C' or 'D' or 'B', */
/*             P is the number of rows of the matrices C and/or D. */
/*             P >= 0. */
/*             Otherwise, P is not used. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the descriptor system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix E of the descriptor system. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,*) */
/*             If JOBB = 'B', the leading N-by-M part of this array must */
/*             contain the input matrix B of the system. */
/*             If JOBB = 'G', the leading N-by-N upper triangular part */
/*             (if UPLO = 'U') or lower triangular part (if UPLO = 'L') */
/*             of this array must contain the upper triangular part or */
/*             lower triangular part, respectively, of the matrix */
/*                   -1 */
/*             G = BR  B'. The stricly lower triangular part (if */
/*             UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If FACT = 'N' or 'D', the leading N-by-N upper triangular */
/*             part (if UPLO = 'U') or lower triangular part (if UPLO = */
/*             'L') of this array must contain the upper triangular part */
/*             or lower triangular part, respectively, of the symmetric */
/*             state weighting matrix Q. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             If FACT = 'C' or 'B', the leading P-by-N part of this */
/*             array must contain the output matrix C of the system. */
/*             If JOBB = 'B' and SCAL = 'G', then Q is modified */
/*             internally, but is restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,N) if FACT = 'N' or 'D'; */
/*             LDQ >= MAX(1,P) if FACT = 'C' or 'B'. */

/*     R       (input) DOUBLE PRECISION array, dimension (LDR,*) */
/*             If FACT = 'N' or 'C', the leading M-by-M upper triangular */
/*             part (if UPLO = 'U') or lower triangular part (if UPLO = */
/*             'L') of this array must contain the upper triangular part */
/*             or lower triangular part, respectively, of the symmetric */
/*             input weighting matrix R. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             If FACT = 'D' or 'B', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D of the */
/*             system. */
/*             If JOBB = 'B' and SCAL = 'G', then R is modified */
/*             internally, but is restored on exit. */
/*             If JOBB = 'G', this array is not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R. */
/*             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C'; */
/*             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B'; */
/*             LDR >= 1        if JOBB = 'G'. */

/*     L       (input) DOUBLE PRECISION array, dimension (LDL,*) */
/*             If JOBL = 'N' and JOBB = 'B', the leading N-by-M part of */
/*             this array must contain the cross weighting matrix L. */
/*             If JOBB = 'B' and SCAL = 'G', then L is modified */
/*             internally, but is restored on exit. */
/*             If JOBL = 'Z' or JOBB = 'G', this array is not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N' and JOBB = 'B'; */
/*             LDL >= 1        if JOBL = 'Z' or  JOBB = 'G'. */

/*     RCONDU  (output) DOUBLE PRECISION */
/*             If N > 0 and INFO = 0 or INFO = 7, an estimate of the */
/*             reciprocal of the condition number (in the 1-norm) of */
/*             the N-th order system of algebraic equations from which */
/*             the solution matrix X is obtained. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             If INFO = 0, the leading N-by-N part of this array */
/*             contains the solution matrix X of the problem. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     ALFAR   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     ALFAI   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (2*N) */
/*             The generalized eigenvalues of the 2N-by-2N matrix pair, */
/*             ordered as specified by SORT (if INFO = 0, or INFO >= 5). */
/*             For instance, if SORT = 'S', the leading N elements of */
/*             these arrays contain the closed-loop spectrum of the */
/*             system. Specifically, */
/*                lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for */
/*             k = 1,2,...,N. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,*) */
/*             The leading 2N-by-2N part of this array contains the */
/*             ordered real Schur form S of the first matrix in the */
/*             reduced matrix pencil associated to the optimal problem, */
/*             corresponding to the scaled Q, R, and L, if JOBB = 'B' */
/*             and SCAL = 'G'. That is, */

/*                    (S   S  ) */
/*                    ( 11  12) */
/*                S = (       ), */
/*                    (0   S  ) */
/*                    (     22) */

/*             where S  , S   and S   are N-by-N matrices. */
/*                    11   12      22 */
/*             Array S must have 2*N+M columns if JOBB = 'B', and 2*N */
/*             columns, otherwise. */

/*     LDS     INTEGER */
/*             The leading dimension of array S. */
/*             LDS >= MAX(1,2*N+M) if JOBB = 'B'; */
/*             LDS >= MAX(1,2*N)   if JOBB = 'G'. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,2*N) */
/*             The leading 2N-by-2N part of this array contains the */
/*             ordered upper triangular form T of the second matrix in */
/*             the reduced matrix pencil associated to the optimal */
/*             problem, corresponding to the scaled Q, R, and L, if */
/*             JOBB = 'B' and SCAL = 'G'. That is, */

/*                    (T   T  ) */
/*                    ( 11  12) */
/*                T = (       ), */
/*                    (0   T  ) */
/*                    (     22) */

/*             where T  , T   and T   are N-by-N matrices. */
/*                    11   12      22 */

/*     LDT     INTEGER */
/*             The leading dimension of array T. */
/*             LDT >= MAX(1,2*N+M) if JOBB = 'B'; */
/*             LDT >= MAX(1,2*N)   if JOBB = 'G'. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N) */
/*             The leading 2N-by-2N part of this array contains the right */
/*             transformation matrix U which reduces the 2N-by-2N matrix */
/*             pencil to the ordered generalized real Schur form (S,T). */
/*             That is, */

/*                    (U   U  ) */
/*                    ( 11  12) */
/*                U = (       ), */
/*                    (U   U  ) */
/*                    ( 21  22) */

/*             where U  , U  , U   and U   are N-by-N matrices. */
/*                    11   12   21      22 */
/*             If JOBB = 'B' and SCAL = 'G', then U corresponds to the */
/*             scaled pencil. If a basis for the stable deflating */
/*             subspace of the original problem is needed, then the */
/*             submatrix U   must be multiplied by the scaling factor */
/*                        21 */
/*             contained in DWORK(4). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,2*N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the original matrix pencil, specifically of the triangular */
/*             M-by-M factor obtained during the reduction process. If */
/*             the user sets TOL > 0, then the given value of TOL is used */
/*             as a lower bound for the reciprocal condition number of */
/*             that matrix; a matrix whose estimated condition number is */
/*             less than 1/TOL is considered to be nonsingular. If the */
/*             user sets TOL <= 0, then a default tolerance, defined by */
/*             TOLDEF = EPS, is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not referenced if JOBB = 'G'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,M,2*N) if JOBB = 'B'; */
/*             LIWORK >= MAX(1,2*N)   if JOBB = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. If JOBB = 'B' and N > 0, DWORK(2) returns the */
/*             reciprocal of the condition number of the M-by-M bottom */
/*             right lower triangular matrix obtained while compressing */
/*             the matrix pencil of order 2N+M to obtain a pencil of */
/*             order 2N. If ACC = 'R', and INFO = 0 or INFO = 7, DWORK(3) */
/*             returns the reciprocal pivot growth factor (see SLICOT */
/*             Library routine MB02PD) for the LU factorization of the */
/*             coefficient matrix of the system of algebraic equations */
/*             giving the solution matrix X; if DWORK(3) is much */
/*             less than 1, then the computed X and RCONDU could be */
/*             unreliable. If INFO = 0 or INFO = 7, DWORK(4) returns the */
/*             scaling factor used to scale Q, R, and L. DWORK(4) is set */
/*             to 1 if JOBB = 'G' or SCAL = 'N'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G'; */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the computed solution may be inaccurate due to poor */
/*                   scaling or eigenvalues too close to the boundary of */
/*                   the stability domain (the imaginary axis, if */
/*                   DICO = 'C', or the unit circle, if DICO = 'D'). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the computed extended matrix pencil is singular, */
/*                   possibly due to rounding errors; */
/*             = 2:  if the QZ algorithm failed; */
/*             = 3:  if reordering of the generalized eigenvalues failed; */
/*             = 4:  if after reordering, roundoff changed values of */
/*                   some complex eigenvalues so that leading eigenvalues */
/*                   in the generalized Schur form no longer satisfy the */
/*                   stability condition; this could also be caused due */
/*                   to scaling; */
/*             = 5:  if the computed dimension of the solution does not */
/*                   equal N; */
/*             = 6:  if the spectrum is too close to the boundary of */
/*                   the stability domain; */
/*             = 7:  if a singular matrix was encountered during the */
/*                   computation of the solution matrix X. */

/*     METHOD */

/*     The routine uses a variant of the method of deflating subspaces */
/*     proposed by van Dooren [1]. See also [2], [3], [4]. */
/*     It is assumed that E is nonsingular, the triple (E,A,B) is */
/*     strongly stabilizable and detectable (see [3]); if, in addition, */

/*        -    [ Q   L ] */
/*        R := [       ] >= 0 , */
/*             [ L'  R ] */

/*     then the pencils */

/*           discrete-time                   continuous-time */

/*     |A   0   B|     |E   0   0|    |A   0   B|     |E   0   0| */
/*     |Q  -E'  L| - z |0  -A'  0| ,  |Q   A'  L| - s |0  -E'  0| ,   (3) */
/*     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0| */

/*     are dichotomic, i.e., they have no eigenvalues on the boundary of */
/*     the stability domain. The above conditions are sufficient for */
/*     regularity of these pencils. A necessary condition is that */
/*     rank([ B'  L'  R']') = m. */

/*     Under these assumptions the algebraic Riccati equation is known to */
/*     have a unique non-negative definite solution. */
/*     The first step in the method of deflating subspaces is to form the */
/*     extended matrices in (3), of order 2N + M. Next, these pencils are */
/*     compressed to a form of order 2N (see [1]) */

/*        lambda x A  - B . */
/*                  f    f */

/*     This generalized eigenvalue problem is then solved using the QZ */
/*     algorithm and the stable deflating subspace Ys is determined. */
/*     If [Y1'|Y2']' is a basis for Ys, then the required solution is */
/*                       -1 */
/*            X = Y2 x Y1  . */

/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         A Generalized Eigenvalue Approach for Solving Riccati */
/*         Equations. */
/*         SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981. */

/*     [2] Arnold, III, W.F. and Laub, A.J. */
/*         Generalized Eigenproblem Algorithms and Software for */
/*         Algebraic Riccati Equations. */
/*         Proc. IEEE, 72, 1746-1754, 1984. */

/*     [3] Mehrmann, V. */
/*         The Autonomous Linear Quadratic Control Problem. Theory and */
/*         Numerical Solution. */
/*         Lect. Notes in Control and Information Sciences, vol. 163, */
/*         Springer-Verlag, Berlin, 1991. */

/*     [4] Sima, V. */
/*         Algorithms for Linear-Quadratic Optimization. */
/*         Pure and Applied Mathematics: A Series of Monographs and */
/*         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */

/*     This routine is particularly suited for systems where the matrix R */
/*     is ill-conditioned, or even singular. */

/*     FURTHER COMMENTS */

/*     To obtain a stabilizing solution of the algebraic Riccati */
/*     equations set SORT = 'S'. */

/*     The routine can also compute the anti-stabilizing solutions of */
/*     the algebraic Riccati equations, by specifying SORT = 'U'. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 2002. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, September 2002, */
/*     December 2002. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    l_dim1 = *ldl;
    l_offset = 1 + l_dim1;
    l -= l_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --alfar;
    --alfai;
    --beta;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    *iwarn = 0;
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    ljobb = lsame_(jobb, "B", (ftnlen)1, (ftnlen)1);
    lfacn = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
    lfacq = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
    lfacr = lsame_(fact, "D", (ftnlen)1, (ftnlen)1);
    lfacb = lsame_(fact, "B", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    lsort = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
    refine = lsame_(acc, "R", (ftnlen)1, (ftnlen)1);
    nn = *n << 1;
    if (ljobb) {
	ljobl = lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1);
	ljobln = lsame_(jobl, "N", (ftnlen)1, (ftnlen)1);
	lscal = lsame_(scal, "G", (ftnlen)1, (ftnlen)1);
	nnm = nn + *m;
/* Computing MAX */
	i__1 = nnm, i__2 = *m * 3;
	ldw = max(i__1,i__2);
    } else {
	lscal = FALSE_;
	nnm = nn;
	ldw = 1;
    }
    np1 = *n + 1;

/*     Test the input scalar arguments. */

    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ljobb && ! lsame_(jobb, "G", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! lfacq && ! lfacr && ! lfacb && ! lfacn) {
	*info = -3;
    } else if (! ljobb || lfacn) {
	if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	    *info = -4;
	}
    }
    if (*info == 0 && ljobb) {
	if (! ljobl && ! ljobln) {
	    *info = -5;
	} else if (! lscal && ! lsame_(scal, "N", (ftnlen)1, (ftnlen)1)) {
	    *info = -6;
	}
    }
    if (*info == 0) {
	if (! lsort && ! lsame_(sort, "U", (ftnlen)1, (ftnlen)1)) {
	    *info = -7;
	} else if (! refine && ! lsame_(acc, "N", (ftnlen)1, (ftnlen)1)) {
	    *info = -8;
	} else if (*n < 0) {
	    *info = -9;
	} else if (ljobb) {
	    if (*m < 0) {
		*info = -10;
	    }
	}
    }
    if (*info == 0 && ! lfacn) {
	if (*p < 0) {
	    *info = -11;
	}
    }
    if (*info == 0) {
	if (*lda < max(1,*n)) {
	    *info = -13;
	} else if (*lde < max(1,*n)) {
	    *info = -15;
	} else if (*ldb < max(1,*n)) {
	    *info = -17;
	} else if ((lfacn || lfacr) && *ldq < max(1,*n) || (lfacq || lfacb) &&
		 *ldq < max(1,*p)) {
	    *info = -19;
	} else if (ljobb) {
	    if ((lfacn || lfacq) && *ldr < max(1,*m) || (lfacr || lfacb) && *
		    ldr < max(1,*p)) {
		*info = -21;
	    } else if (ljobln && *ldl < max(1,*n) || ljobl && *ldl < 1) {
		*info = -23;
	    }
	} else {
	    if (*ldr < 1) {
		*info = -21;
	    } else if (*ldl < 1) {
		*info = -23;
	    }
	}
    }
    if (*info == 0) {
	if (*ldx < max(1,*n)) {
	    *info = -26;
	} else if (*lds < max(1,nnm)) {
	    *info = -31;
	} else if (*ldt < max(1,nnm)) {
	    *info = -33;
	} else if (*ldu < max(1,nn)) {
	    *info = -35;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2);
	    if (*ldwork < max(i__1,ldw)) {
		*info = -39;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SG02AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 4.;
	dwork[4] = 1.;
	return 0;
    }

/*     Start computations. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    lscal = lscal && ljobb;
    if (lscal) {

/*        Scale the matrices Q, R (or G), and L so that */
/*           norm(Q) + norm(R) + norm(L) = 1, */
/*        using the 1-norm. If Q and/or R are factored, the norms of */
/*        the factors are used. */
/*        Workspace: need   max(N,M), if FACT = 'N'; */
/*                          N,        if FACT = 'D'; */
/*                          M,        if FACT = 'C'. */

	if (lfacn || lfacr) {
	    scale = dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6, (ftnlen)1);
	    *(unsigned char *)qtype = *(unsigned char *)uplo;
	    np = *n;
	} else {
	    scale = dlange_("1-norm", p, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6);
	    *(unsigned char *)qtype = 'G';
	    np = *p;
	}

	if (lfacn || lfacq) {
	    rnorm = dlansy_("1-norm", uplo, m, &r__[r_offset], ldr, &dwork[1],
		     (ftnlen)6, (ftnlen)1);
	    *(unsigned char *)rtype = *(unsigned char *)uplo;
	    mp = *m;
	} else {
	    rnorm = dlange_("1-norm", p, m, &r__[r_offset], ldr, &dwork[1], (
		    ftnlen)6);
	    *(unsigned char *)rtype = 'G';
	    mp = *p;
	}
	scale += rnorm;

	if (ljobln) {
	    scale += dlange_("1-norm", n, m, &l[l_offset], ldl, &dwork[1], (
		    ftnlen)6);
	}
	if (scale == 0.) {
	    scale = 1.;
	}

	dlascl_(qtype, &c__0, &c__0, &scale, &c_b30, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
	dlascl_(rtype, &c__0, &c__0, &scale, &c_b30, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
	if (ljobln) {
	    dlascl_("G", &c__0, &c__0, &scale, &c_b30, n, m, &l[l_offset], 
		    ldl, &info1, (ftnlen)1);
	}
    } else {
	scale = 1.;
    }

/*     Construct the extended matrix pair. */
/*     Workspace: need   1,                if JOBB = 'G', */
/*                       max(1,2*N+M,3*M), if JOBB = 'B'; */
/*                prefer larger. */

    sb02oy_("Optimal control", dico, jobb, fact, uplo, jobl, "Not identity E",
	     n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq,
	     &r__[r_offset], ldr, &l[l_offset], ldl, &e[e_offset], lde, &s[
	    s_offset], lds, &t[t_offset], ldt, tol, &iwork[1], &dwork[1], 
	    ldwork, info, (ftnlen)15, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)14);

    if (lscal) {

/*        Undo scaling of the data arrays. */

	dlascl_(qtype, &c__0, &c__0, &c_b30, &scale, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
	dlascl_(rtype, &c__0, &c__0, &c_b30, &scale, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
	if (ljobln) {
	    dlascl_("G", &c__0, &c__0, &c_b30, &scale, n, m, &l[l_offset], 
		    ldl, &info1, (ftnlen)1);
	}
    }

    if (*info != 0) {
	return 0;
    }
    wrkopt = (integer) dwork[1];
    if (ljobb) {
	rcondl = dwork[2];
    }

/*     Workspace: need   max(7*(2*N+1)+16,16*N); */
/*                prefer larger. */

    if (discr) {
	if (lsort) {

/*           The natural tendency of the QZ algorithm to get the largest */
/*           eigenvalues in the leading part of the matrix pair is */
/*           exploited, by computing the unstable eigenvalues of the */
/*           permuted matrix pair. */

	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ov_, &nn, &t[
		    t_offset], ldt, &s[s_offset], lds, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
	    dswap_(n, &alfar[np1], &c__1, &alfar[1], &c__1);
	    dswap_(n, &alfai[np1], &c__1, &alfai[1], &c__1);
	    dswap_(n, &beta[np1], &c__1, &beta[1], &c__1);
	} else {
	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ov_, &nn, &s[
		    s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
	}
    } else {
	if (lsort) {
	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ow_, &nn, &s[
		    s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
	} else {
	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ou_, &nn, &s[
		    s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
	}
    }
    if (info1 > 0 && info1 <= nn + 1) {
	*info = 2;
    } else if (info1 == nn + 2) {
	*info = 4;
    } else if (info1 == nn + 3) {
	*info = 3;
    } else if (ndim != *n) {
	*info = 5;
    }
    if (*info != 0) {
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[1];
    wrkopt = max(i__1,i__2);

/*     Take the non-identity matrix E into account and orthogonalize the */
/*     basis. Use the array X as workspace. */
/*     Workspace: need   N; */
/*                prefer N*NB. */

    dgemm_("No transpose", "No transpose", n, n, n, &c_b30, &e[e_offset], lde,
	     &u[u_offset], ldu, &c_b71, &x[x_offset], ldx, (ftnlen)12, (
	    ftnlen)12);
    dlacpy_("Full", n, n, &x[x_offset], ldx, &u[u_offset], ldu, (ftnlen)4);
    dgeqrf_(&nn, n, &u[u_offset], ldu, &x[x_offset], &dwork[1], ldwork, &
	    info1);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[1];
    wrkopt = max(i__1,i__2);
    dorgqr_(&nn, n, n, &u[u_offset], ldu, &x[x_offset], &dwork[1], ldwork, &
	    info1);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[1];
    wrkopt = max(i__1,i__2);

/*     Check for the symmetry of the solution. The array X is again used */
/*     as workspace. */

    dgemm_("Transpose", "No transpose", n, n, n, &c_b30, &u[u_offset], ldu, &
	    u[np1 + u_dim1], ldu, &c_b71, &x[x_offset], ldx, (ftnlen)9, (
	    ftnlen)12);
    u12m = 0.;
    asym = 0.;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = u12m, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
	    u12m = max(d__2,d__3);
/* Computing MAX */
	    d__2 = asym, d__3 = (d__1 = x[i__ + j * x_dim1] - x[j + i__ * 
		    x_dim1], abs(d__1));
	    asym = max(d__2,d__3);
/* L10: */
	}

/* L20: */
    }

    eps = dlamch_("Epsilon", (ftnlen)7);
    seps = sqrt(eps);
    asym -= seps;
    if (asym > u12m * .1) {
	*info = 6;
	return 0;
    } else if (asym > seps) {
	*iwarn = 1;
    }

/*     Compute the solution of X*U(1,1) = U(2,1). Use the (2,1) block */
/*     of S as a workspace for factoring U(1,1). */

    if (refine) {

/*        Use LU factorization and iterative refinement for finding X. */
/*        Workspace:  need   8*N. */

/*        First transpose U(2,1) in-situ. */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__;
	    dswap_(&i__2, &u[*n + i__ + (i__ + 1) * u_dim1], ldu, &u[*n + i__ 
		    + 1 + i__ * u_dim1], &c__1);
/* L30: */
	}

	iwr = 1;
	iwc = iwr + *n;
	iwf = iwc + *n;
	iwb = iwf + *n;
	iw = iwb + *n;

	mb02pd_("Equilibrate", "Transpose", n, n, &u[u_offset], ldu, &s[np1 + 
		s_dim1], lds, &iwork[1], equed, &dwork[iwr], &dwork[iwc], &u[
		np1 + u_dim1], ldu, &x[x_offset], ldx, rcondu, &dwork[iwf], &
		dwork[iwb], &iwork[np1], &dwork[iw], &info1, (ftnlen)11, (
		ftnlen)9, (ftnlen)1);

/*        Transpose U(2,1) back in-situ. */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__;
	    dswap_(&i__2, &u[*n + i__ + (i__ + 1) * u_dim1], ldu, &u[*n + i__ 
		    + 1 + i__ * u_dim1], &c__1);
/* L40: */
	}

	if (! lsame_(equed, "N", (ftnlen)1, (ftnlen)1)) {

/*           Undo the equilibration of U(1,1) and U(2,1). */

	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);

	    if (rowequ) {

		i__1 = *n - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
		    dwork[iwr + i__] = 1. / dwork[iwr + i__];
/* L50: */
		}

		mb01sd_("Row scaling", n, n, &u[u_offset], ldu, &dwork[iwr], &
			dwork[iwc], (ftnlen)11);
	    }

	    if (colequ) {

		i__1 = *n - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
		    dwork[iwc + i__] = 1. / dwork[iwc + i__];
/* L60: */
		}

		mb01sd_("Column scaling", &nn, n, &u[u_offset], ldu, &dwork[
			iwr], &dwork[iwc], (ftnlen)14);
	    }
	}

	pivotu = dwork[iw];

	if (info1 > 0) {

/*           Singular matrix. Set INFO and DWORK for error return. */

	    *info = 7;
	    goto L80;
	}

    } else {

/*        Use LU factorization and a standard solution algorithm. */

	dlacpy_("Full", n, n, &u[u_offset], ldu, &s[np1 + s_dim1], lds, (
		ftnlen)4);
	dlacpy_("Full", n, n, &u[np1 + u_dim1], ldu, &x[x_offset], ldx, (
		ftnlen)4);

/*        Solve the system X*U(1,1) = U(2,1). */

	mb02vd_("No Transpose", n, n, &s[np1 + s_dim1], lds, &iwork[1], &x[
		x_offset], ldx, &info1, (ftnlen)12);

	if (info1 != 0) {
	    *info = 7;
	    *rcondu = 0.;
	    goto L80;
	} else {

/*           Compute the norm of U(1,1). */

	    unorm = dlange_("1-norm", n, n, &u[u_offset], ldu, &dwork[1], (
		    ftnlen)6);

/*           Estimate the reciprocal condition of U(1,1). */
/*           Workspace: need 4*N. */

	    dgecon_("1-norm", n, &s[np1 + s_dim1], lds, &unorm, rcondu, &
		    dwork[1], &iwork[np1], info, (ftnlen)6);

	    if (*rcondu < eps) {

/*              Nearly singular matrix. Set IWARN for warning indication. */

		*iwarn = 1;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = *n << 2;
	    wrkopt = max(i__1,i__2);
	}
    }

/*     Set S(2,1) to zero. */

    dlaset_("Full", n, n, &c_b71, &c_b71, &s[np1 + s_dim1], lds, (ftnlen)4);

/*     Make sure the solution matrix X is symmetric. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__;
	daxpy_(&i__2, &c_b30, &x[i__ + (i__ + 1) * x_dim1], ldx, &x[i__ + 1 + 
		i__ * x_dim1], &c__1);
	i__2 = *n - i__;
	dscal_(&i__2, &c_b107, &x[i__ + 1 + i__ * x_dim1], &c__1);
	i__2 = *n - i__;
	dcopy_(&i__2, &x[i__ + 1 + i__ * x_dim1], &c__1, &x[i__ + (i__ + 1) * 
		x_dim1], ldx);
/* L70: */
    }

    if (lscal) {

/*        Undo scaling for the solution X. */

	dlascl_("G", &c__0, &c__0, &c_b30, &scale, n, n, &x[x_offset], ldx, &
		info1, (ftnlen)1);
    }

    dwork[1] = (doublereal) wrkopt;

L80:
    if (ljobb) {
	dwork[2] = rcondl;
    }
    if (refine) {
	dwork[3] = pivotu;
    }
    dwork[4] = scale;

    return 0;
/* *** Last line of SG02AD *** */
} /* sg02ad_ */

