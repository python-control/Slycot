/* SB02OD.f -- translated by f2c (version 20100827).
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
static doublereal c_b26 = 1.;
static integer c__1 = 1;
static doublereal c_b66 = -1.;
static doublereal c_b104 = 0.;

/* Subroutine */ int sb02od_(char *dico, char *jobb, char *fact, char *uplo, 
	char *jobl, char *sort, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, doublereal *l, 
	integer *ldl, doublereal *rcond, doublereal *x, integer *ldx, 
	doublereal *alfar, doublereal *alfai, doublereal *beta, doublereal *s,
	 integer *lds, doublereal *t, integer *ldt, doublereal *u, integer *
	ldu, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *info, ftnlen dico_len, ftnlen 
	jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen 
	sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, l_dim1, l_offset, q_dim1, 
	    q_offset, r_dim1, r_offset, s_dim1, s_offset, t_dim1, t_offset, 
	    u_dim1, u_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn, mp, np, np1, ldw, nnm;
    static doublereal dum[1];
    static integer ndim;
    static logical lscl;
    static integer info1;
    static logical lfacb, lfacn, lfacq, lfacr, ljobb;
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgges_(char *, char *, char *, L_fp, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *, ftnlen, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobl;
    static doublereal qscal;
    static logical discr;
    static doublereal rscal;
    extern logical sb02mr_(), sb02mv_(), sb02ou_(), sb02ov_(), sb02ow_();
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
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobln;
    static doublereal rcondl;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
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

/*     To solve for X either the continuous-time algebraic Riccati */
/*     equation */
/*                              -1 */
/*        Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1) */

/*     or the discrete-time algebraic Riccati equation */
/*                                     -1 */
/*        X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2) */

/*     where A, B, Q, R, and L are N-by-N, N-by-M, N-by-N, M-by-M and */
/*     N-by-M matrices, respectively, such that Q = C'C, R = D'D and */
/*     L = C'D; X is an N-by-N symmetric matrix. */
/*     The routine also returns the computed values of the closed-loop */
/*     spectrum of the system, i.e., the stable eigenvalues lambda(1), */
/*     ..., lambda(N) of the corresponding Hamiltonian or symplectic */
/*     pencil, in the continuous-time case or discrete-time case, */
/*     respectively. */
/*                              -1 */
/*     Optionally, matrix G = BR  B' may be given instead of B and R. */
/*     Other options include the case with Q and/or R given in a */
/*     factored form, Q = C'C, R = D'D, and with L a zero matrix. */

/*     The routine uses the method of deflating subspaces, based on */
/*     reordering the eigenvalues in a generalized Schur matrix pair. */
/*     A standard eigenproblem is solved in the continuous-time case */
/*     if G is given. */

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
/*             the matrices G and Q (if FACT = 'N'), or Q and R (if */
/*             JOBB = 'B'), is stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */
/*             JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed. */
/*             SLICOT Library routine SB02MT should be called just before */
/*             SB02OD, for obtaining the results when JOBB = 'G' and */
/*             JOBL = 'N'. */

/*     SORT    CHARACTER*1 */
/*             Specifies which eigenvalues should be obtained in the top */
/*             of the generalized Schur form, as follows: */
/*             = 'S':  Stable   eigenvalues come first; */
/*             = 'U':  Unstable eigenvalues come first. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the matrices */
/*             A, Q, and X, and the number of rows of the matrices B */
/*             and L.  N >= 0. */

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
/*             state matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

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
/*             If JOBB = 'B', the triangular part of this array defined */
/*             by UPLO is modified internally, but is restored on exit. */
/*             If FACT = 'C' or 'B', the leading P-by-N part of this */
/*             array must contain the output matrix C of the system. */
/*             If JOBB = 'B', this part is modified internally, but is */
/*             restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,N) if FACT = 'N' or 'D', */
/*             LDQ >= MAX(1,P) if FACT = 'C' or 'B'. */

/*     R       (input) DOUBLE PRECISION array, dimension (LDR,M) */
/*             If FACT = 'N' or 'C', the leading M-by-M upper triangular */
/*             part (if UPLO = 'U') or lower triangular part (if UPLO = */
/*             'L') of this array must contain the upper triangular part */
/*             or lower triangular part, respectively, of the symmetric */
/*             input weighting matrix R. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             The triangular part of this array defined by UPLO is */
/*             modified internally, but is restored on exit. */
/*             If FACT = 'D' or 'B', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D of the */
/*             system. This part is modified internally, but is restored */
/*             on exit. */
/*             If JOBB = 'G', this array is not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R. */
/*             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C'; */
/*             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B'; */
/*             LDR >= 1        if JOBB = 'G'. */

/*     L       (input) DOUBLE PRECISION array, dimension (LDL,M) */
/*             If JOBL = 'N' (and JOBB = 'B'), the leading N-by-M part of */
/*             this array must contain the cross weighting matrix L. */
/*             This part is modified internally, but is restored on exit. */
/*             If JOBL = 'Z' or JOBB = 'G', this array is not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N' and JOBB = 'B'; */
/*             LDL >= 1        if JOBL = 'Z' or  JOBB = 'G'. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             An estimate of the reciprocal of the condition number (in */
/*             the 1-norm) of the N-th order system of algebraic */
/*             equations from which the solution matrix X is obtained. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the */
/*             solution matrix X of the problem. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     ALFAR   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     ALFAI   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (2*N) */
/*             The generalized eigenvalues of the 2N-by-2N matrix pair, */
/*             ordered as specified by SORT (if INFO = 0). For instance, */
/*             if SORT = 'S', the leading N elements of these arrays */
/*             contain the closed-loop spectrum of the system matrix */
/*             A - BF, where F is the optimal feedback matrix computed */
/*             based on the solution matrix X. Specifically, */
/*                lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for */
/*             k = 1,2,...,N. */
/*             If DICO = 'C' and JOBB = 'G', the elements of BETA are */
/*             set to 1. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,*) */
/*             The leading 2N-by-2N part of this array contains the */
/*             ordered real Schur form S of the first matrix in the */
/*             reduced matrix pencil associated to the optimal problem, */
/*             or of the corresponding Hamiltonian matrix, if DICO = 'C' */
/*             and JOBB = 'G'. That is, */

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
/*             LDS >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDS >= MAX(1,2*N)   if JOBB = 'G'. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,2*N) */
/*             If DICO = 'D' or JOBB = 'B', the leading 2N-by-2N part of */
/*             this array contains the ordered upper triangular form T of */
/*             the second matrix in the reduced matrix pencil associated */
/*             to the optimal problem. That is, */

/*                    (T   T  ) */
/*                    ( 11  12) */
/*                T = (       ), */
/*                    (0   T  ) */
/*                    (     22) */

/*             where T  , T   and T   are N-by-N matrices. */
/*                    11   12      22 */
/*             If DICO = 'C' and JOBB = 'G' this array is not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of array T. */
/*             LDT >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDT >= MAX(1,2*N)   if JOBB = 'G' and DICO = 'D', */
/*             LDT >= 1            if JOBB = 'G' and DICO = 'C'. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N) */
/*             The leading 2N-by-2N part of this array contains the right */
/*             transformation matrix U which reduces the 2N-by-2N matrix */
/*             pencil to the ordered generalized real Schur form (S,T), */
/*             or the Hamiltonian matrix to the ordered real Schur */
/*             form S, if DICO = 'C' and JOBB = 'G'. That is, */

/*                    (U   U  ) */
/*                    ( 11  12) */
/*                U = (       ), */
/*                    (U   U  ) */
/*                    ( 21  22) */

/*             where U  , U  , U   and U   are N-by-N matrices. */
/*                    11   12   21      22 */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,2*N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the original matrix pencil, specifically of the triangular */
/*             factor obtained during the reduction process. If the user */
/*             sets TOL > 0, then the given value of TOL is used as a */
/*             lower bound for the reciprocal condition number of that */
/*             matrix; a matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. If the user */
/*             sets TOL <= 0, then a default tolerance, defined by */
/*             TOLDEF = EPS, is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not referenced if JOBB = 'G'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,M,2*N) if JOBB = 'B', */
/*             LIWORK >= MAX(1,2*N)   if JOBB = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. If JOBB = 'B' and N > 0, DWORK(2) returns the */
/*             reciprocal of the condition number of the M-by-M lower */
/*             triangular matrix obtained after compressing the matrix */
/*             pencil of order 2N+M to obtain a pencil of order 2N. */
/*             If INFO = 0 or INFO = 6, DWORK(3) returns the scaling */
/*             factor used internally, which should multiply the */
/*             submatrix Y2 to recover X from the first N columns of U */
/*             (see METHOD). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(3,6*N),                       if JOBB = 'G', */
/*                                                            DICO = 'C'; */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G', */
/*                                                            DICO = 'D'; */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the computed extended matrix pencil is singular, */
/*                   possibly due to rounding errors; */
/*             = 2:  if the QZ (or QR) algorithm failed; */
/*             = 3:  if reordering of the (generalized) eigenvalues */
/*                   failed; */
/*             = 4:  if after reordering, roundoff changed values of */
/*                   some complex eigenvalues so that leading eigenvalues */
/*                   in the (generalized) Schur form no longer satisfy */
/*                   the stability condition; this could also be caused */
/*                   due to scaling; */
/*             = 5:  if the computed dimension of the solution does not */
/*                   equal N; */
/*             = 6:  if a singular matrix was encountered during the */
/*                   computation of the solution matrix X. */

/*     METHOD */

/*     The routine uses a variant of the method of deflating subspaces */
/*     proposed by van Dooren [1]. See also [2], [3]. */
/*     It is assumed that (A,B) is stabilizable and (C,A) is detectable. */
/*     Under these assumptions the algebraic Riccati equation is known to */
/*     have a unique non-negative definite solution. */
/*     The first step in the method of deflating subspaces is to form the */
/*     extended Hamiltonian matrices, dimension 2N + M given by */

/*           discrete-time                   continuous-time */

/*     |A   0   B|     |I   0   0|    |A   0   B|     |I   0   0| */
/*     |Q  -I   L| - z |0  -A'  0|,   |Q   A'  L| - s |0  -I   0|. */
/*     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0| */

/*     Next, these pencils are compressed to a form (see [1]) */

/*        lambda x A  - B . */
/*                  f    f */

/*     This generalized eigenvalue problem is then solved using the QZ */
/*     algorithm and the stable deflating subspace Ys is determined. */
/*     If [Y1'|Y2']' is a basis for Ys, then the required solution is */
/*                       -1 */
/*            X = Y2 x Y1  . */
/*     A standard eigenvalue problem is solved using the QR algorithm in */
/*     the continuous-time case when G is given (DICO = 'C', JOBB = 'G'). */

/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         A Generalized Eigenvalue Approach for Solving Riccati */
/*         Equations. */
/*         SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981. */

/*     [2] Mehrmann, V. */
/*         The Autonomous Linear Quadratic Control Problem. Theory and */
/*         Numerical Solution. */
/*         Lect. Notes in Control and Information Sciences, vol. 163, */
/*         Springer-Verlag, Berlin, 1991. */

/*     [3] Sima, V. */
/*         Algorithms for Linear-Quadratic Optimization. */
/*         Pure and Applied Mathematics: A Series of Monographs and */
/*         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */

/*     This routine is particularly suited for systems where the matrix R */
/*     is ill-conditioned. Internal scaling is used. */

/*     FURTHER COMMENTS */

/*     To obtain a stabilizing solution of the algebraic Riccati */
/*     equations set SORT = 'S'. */

/*     The routine can also compute the anti-stabilizing solutions of */
/*     the algebraic Riccati equations, by specifying SORT = 'U'. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB02CD by T.G.J. Beelen, Philips, */
/*     Eindhoven, Holland. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, June 2002, */
/*     December 2002, January 2005. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

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

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    ljobb = lsame_(jobb, "B", (ftnlen)1, (ftnlen)1);
    lfacn = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
    lfacq = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
    lfacr = lsame_(fact, "D", (ftnlen)1, (ftnlen)1);
    lfacb = lsame_(fact, "B", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    lsort = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

    nn = *n << 1;
    if (ljobb) {
	ljobl = lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1);
	ljobln = lsame_(jobl, "N", (ftnlen)1, (ftnlen)1);
	nnm = nn + *m;
/* Computing MAX */
	i__1 = nnm, i__2 = *m * 3;
	ldw = max(i__1,i__2);
    } else {
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
	}
    }
    if (*info == 0) {
	if (! lsort && ! lsame_(sort, "U", (ftnlen)1, (ftnlen)1)) {
	    *info = -6;
	} else if (*n < 0) {
	    *info = -7;
	} else if (ljobb) {
	    if (*m < 0) {
		*info = -8;
	    }
	}
    }
    if (*info == 0 && ! lfacn) {
	if (*p < 0) {
	    *info = -9;
	}
    }
    if (*info == 0) {
	if (*lda < max(1,*n)) {
	    *info = -11;
	} else if (*ldb < max(1,*n)) {
	    *info = -13;
	} else if ((lfacn || lfacr) && *ldq < max(1,*n) || (lfacq || lfacb) &&
		 *ldq < max(1,*p)) {
	    *info = -15;
	} else if (*ldr < 1) {
	    *info = -17;
	} else if (*ldl < 1) {
	    *info = -19;
	} else if (ljobb) {
	    if ((lfacn || lfacq) && *ldr < *m || (lfacr || lfacb) && *ldr < *
		    p) {
		*info = -17;
	    } else if (ljobln && *ldl < *n) {
		*info = -19;
	    }
	}
    }
    if (*info == 0) {
	if (*ldx < max(1,*n)) {
	    *info = -22;
	} else if (*lds < max(1,nnm)) {
	    *info = -27;
	} else if (*ldt < 1) {
	    *info = -29;
	} else if (*ldu < max(1,nn)) {
	    *info = -31;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 3, i__2 = *n * 6;
	    if (*ldwork < max(i__1,i__2)) {
		*info = -35;
	    } else if (discr || ljobb) {
		if (*ldt < nnm) {
		    *info = -29;
		} else /* if(complicated condition) */ {
/* Computing MAX */
		    i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2)
			    ;
		    if (*ldwork < max(i__1,ldw)) {
			*info = -35;
		    }
		}
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB02OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*rcond = 1.;
	dwork[1] = 3.;
	dwork[3] = 1.;
	return 0;
    }

/*     Always scale the matrix pencil. */

    lscal = TRUE_;

/*     Start computations. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    if (lscal && ljobb) {

/*        Scale the matrices Q, R, and L so that */
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

	if (lfacn || lfacr) {
	    qscal = scale;
	} else {
	    qscal = sqrt(scale);
	}

	if (lfacn || lfacq) {
	    rscal = scale;
	} else {
	    rscal = sqrt(scale);
	}

	dlascl_(qtype, &c__0, &c__0, &qscal, &c_b26, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
	dlascl_(rtype, &c__0, &c__0, &rscal, &c_b26, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
	if (ljobln) {
	    dlascl_("G", &c__0, &c__0, &scale, &c_b26, n, m, &l[l_offset], 
		    ldl, &info1, (ftnlen)1);
	}
    }

/*     Construct the extended matrix pair. */

/*     Workspace: need   1,                if JOBB = 'G', */
/*                       max(1,2*N+M,3*M), if JOBB = 'B'; */
/*                prefer larger. */

    sb02oy_("Optimal control", dico, jobb, fact, uplo, jobl, "Identity E", n, 
	    m, p, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
	    r__[r_offset], ldr, &l[l_offset], ldl, &u[u_offset], &c__1, &s[
	    s_offset], lds, &t[t_offset], ldt, tol, &iwork[1], &dwork[1], 
	    ldwork, info, (ftnlen)15, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)10);

    if (lscal && ljobb) {

/*        Undo scaling of the data arrays. */

	dlascl_(qtype, &c__0, &c__0, &c_b26, &qscal, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
	dlascl_(rtype, &c__0, &c__0, &c_b26, &rscal, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
	if (ljobln) {
	    dlascl_("G", &c__0, &c__0, &c_b26, &scale, n, m, &l[l_offset], 
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

    if (lscal && ! ljobb) {

/*        This part of the code is used when G is given (JOBB = 'G'). */
/*        A standard eigenproblem is solved in the continuous-time case. */
/*        Scale the Hamiltonian matrix S, if DICO = 'C', or the */
/*        symplectic pencil (S,T), if DICO = 'D', using the square roots */
/*        of the norms of the matrices Q and G. */
/*        Workspace: need   N. */

	if (lfacn || lfacr) {
	    scale = sqrt(dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[
		    1], (ftnlen)6, (ftnlen)1));
	} else {
	    scale = dlange_("1-norm", p, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6);
	}
	rnorm = sqrt(dlansy_("1-norm", uplo, n, &b[b_offset], ldb, &dwork[1], 
		(ftnlen)6, (ftnlen)1));

	lscl = min(scale,rnorm) > 0. && scale != rnorm;

	if (lscl) {
	    if (discr) {
		dlascl_("G", &c__0, &c__0, &scale, &rnorm, n, n, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
		dlascl_("G", &c__0, &c__0, &rnorm, &scale, n, n, &t[np1 * 
			t_dim1 + 1], ldt, &info1, (ftnlen)1);
	    } else {
		d__1 = -rnorm;
		dlascl_("G", &c__0, &c__0, &scale, &d__1, n, n, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
		dlascl_("G", &c__0, &c__0, &rnorm, &scale, n, n, &s[np1 * 
			s_dim1 + 1], lds, &info1, (ftnlen)1);
		dlascl_("G", &c__0, &c__0, &c_b26, &c_b66, n, n, &s[np1 + np1 
			* s_dim1], lds, &info1, (ftnlen)1);
	    }
	} else {
	    if (! discr) {
		dlascl_("G", &c__0, &c__0, &c_b26, &c_b66, n, &nn, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
	    }
	}
    } else {
	lscl = FALSE_;
    }

/*     Workspace: need   max(7*(2*N+1)+16,16*N), */
/*                                          if JOBB = 'B' or  DICO = 'D'; */
/*                       6*N,               if JOBB = 'G' and DICO = 'C'; */
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
	if (ljobb) {
	    if (lsort) {
		dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ow_, &nn, &
			s[s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1],
			 &alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset],
			 ldu, &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)
			10, (ftnlen)7, (ftnlen)4);
	    } else {
		dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ou_, &nn, &
			s[s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1],
			 &alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset],
			 ldu, &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)
			10, (ftnlen)7, (ftnlen)4);
	    }
	} else {
	    if (lsort) {
		dgees_("Vectors", "Sort", (L_fp)sb02mv_, &nn, &s[s_offset], 
			lds, &ndim, &alfar[1], &alfai[1], &u[u_offset], ldu, &
			dwork[1], ldwork, &bwork[1], &info1, (ftnlen)7, (
			ftnlen)4);
	    } else {
		dgees_("Vectors", "Sort", (L_fp)sb02mr_, &nn, &s[s_offset], 
			lds, &ndim, &alfar[1], &alfai[1], &u[u_offset], ldu, &
			dwork[1], ldwork, &bwork[1], &info1, (ftnlen)7, (
			ftnlen)4);
	    }
	    dum[0] = 1.;
	    dcopy_(&nn, dum, &c__0, &beta[1], &c__1);
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

/*     Select submatrices U1 and U2 out of the array U which define the */
/*     solution X = U2 x inv(U1). */
/*     Since X = X' we may obtain X as the solution of the system of */
/*     linear equations U1' x X = U2', where */
/*        U1 = U(1:n, 1:n), */
/*        U2 = U(n+1:2n, 1:n). */
/*     Use the (2,1) block of S as a workspace for factoring U1. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dcopy_(n, &u[np1 + j * u_dim1], &c__1, &x[j + x_dim1], ldx);
/* L20: */
    }

    dlacpy_("Full", n, n, &u[u_offset], ldu, &s[np1 + s_dim1], lds, (ftnlen)4)
	    ;

/*     Check if U1 is singular. */

    unorm = dlange_("1-norm", n, n, &s[np1 + s_dim1], lds, &dwork[1], (ftnlen)
	    6);

/*     Solve the system U1' x X = U2'. */

    dgetrf_(n, n, &s[np1 + s_dim1], lds, &iwork[1], &info1);
    if (info1 != 0) {
	*info = 6;
	dwork[3] = 1.;
	if (lscal) {
	    if (ljobb) {
		dwork[3] = scale;
	    } else if (lscl) {
		dwork[3] = scale / rnorm;
	    }
	}
	return 0;
    } else {

/*        Estimate the reciprocal condition of U1. */
/*        Workspace: need 3*N. */

	dgecon_("1-norm", n, &s[np1 + s_dim1], lds, &unorm, rcond, &dwork[1], 
		&iwork[np1], info, (ftnlen)6);

	if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {

/*           Nearly singular matrix.  Set INFO for error return. */

	    *info = 6;
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = *n * 3;
	wrkopt = max(i__1,i__2);
	dgetrs_("Transpose", n, n, &s[np1 + s_dim1], lds, &iwork[1], &x[
		x_offset], ldx, &info1, (ftnlen)9);

/*        Set S(2,1) to zero. */

	dlaset_("Full", n, n, &c_b104, &c_b104, &s[np1 + s_dim1], lds, (
		ftnlen)4);

	if (lscal) {

/*           Prepare to undo scaling for the solution X. */

	    if (! ljobb) {
		if (lscl) {
		    scale /= rnorm;
		} else {
		    scale = 1.;
		}
	    }
	    dwork[3] = scale;
	    scale *= .5;
	} else {
	    dwork[3] = 1.;
	    scale = .5;
	}

/*        Make sure the solution matrix X is symmetric. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__ + 1;
	    daxpy_(&i__2, &c_b26, &x[i__ + i__ * x_dim1], ldx, &x[i__ + i__ * 
		    x_dim1], &c__1);
	    i__2 = *n - i__ + 1;
	    dscal_(&i__2, &scale, &x[i__ + i__ * x_dim1], &c__1);
	    i__2 = *n - i__ + 1;
	    dcopy_(&i__2, &x[i__ + i__ * x_dim1], &c__1, &x[i__ + i__ * 
		    x_dim1], ldx);
/* L40: */
	}
    }

    dwork[1] = (doublereal) wrkopt;
    if (ljobb) {
	dwork[2] = rcondl;
    }

    return 0;
/* *** Last line of SB02OD *** */
} /* sb02od_ */

