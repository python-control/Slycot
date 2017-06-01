/* SB02OY.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b26 = 1.;
static doublereal c_b27 = 0.;
static doublereal c_b46 = -1.;

/* Subroutine */ int sb02oy_(char *type__, char *dico, char *jobb, char *fact,
	 char *uplo, char *jobl, char *jobe, integer *n, integer *m, integer *
	p, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *q, integer *ldq, doublereal *r__, integer *ldr, 
	doublereal *l, integer *ldl, doublereal *e, integer *lde, doublereal *
	af, integer *ldaf, doublereal *bf, integer *ldbf, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen type_len, ftnlen dico_len, ftnlen jobb_len, ftnlen fact_len, 
	ftnlen uplo_len, ftnlen jobl_len, ftnlen jobe_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, bf_dim1, 
	    bf_offset, e_dim1, e_offset, l_dim1, l_offset, q_dim1, q_offset, 
	    r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, nm, np1, n2p1, nnm, itau;
    static logical optc, lfacb, lfacn, lfacq, lfacr, ljobb, ljobe;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobl, discr;
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical luplo;
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqlf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrcon_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), dormql_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
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

/*     To construct the extended matrix pairs for the computation of the */
/*     solution of the algebraic matrix Riccati equations arising in the */
/*     problems of optimal control, both discrete and continuous-time, */
/*     and of spectral factorization, both discrete and continuous-time. */
/*     These matrix pairs, of dimension 2N + M, are given by */

/*           discrete-time                   continuous-time */

/*     |A   0   B|     |E   0   0|    |A   0   B|     |E   0   0| */
/*     |Q  -E'  L| - z |0  -A'  0|,   |Q   A'  L| - s |0  -E'  0|.   (1) */
/*     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0| */

/*     After construction, these pencils are compressed to a form */
/*     (see [1]) */

/*        lambda x A  - B , */
/*                  f    f */

/*     where A  and B  are 2N-by-2N matrices. */
/*            f      f */
/*                              -1 */
/*     Optionally, matrix G = BR  B' may be given instead of B and R; */
/*     then, for L = 0, 2N-by-2N matrix pairs are directly constructed as */

/*         discrete-time            continuous-time */

/*     |A   0 |     |E   G |    |A  -G |     |E   0 | */
/*     |      | - z |      |,   |      | - s |      |.               (2) */
/*     |Q  -E'|     |0  -A'|    |Q   A'|     |0  -E'| */

/*     Similar pairs are obtained for non-zero L, if SLICOT Library */
/*     routine SB02MT is called before SB02OY. */
/*     Other options include the case with E identity matrix, L a zero */
/*     matrix, or Q and/or R given in a factored form, Q = C'C, R = D'D. */
/*     For spectral factorization problems, there are minor differences */
/*     (e.g., B is replaced by C'). */
/*     The second matrix in (2) is not constructed in the continuous-time */
/*     case if E is specified as being an identity matrix. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             Specifies the type of problem to be addressed as follows: */
/*             = 'O':  Optimal control problem; */
/*             = 'S':  Spectral factorization problem. */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of linear system considered as follows: */
/*             = 'C':  Continuous-time system; */
/*             = 'D':  Discrete-time system. */

/*     JOBB    CHARACTER*1 */
/*             Specifies whether or not the matrix G is given, instead */
/*             of the matrices B and R, as follows: */
/*             = 'B':  B and R are given; */
/*             = 'G':  G is given. */
/*             For JOBB = 'G', a 2N-by-2N matrix pair is directly */
/*             obtained assuming L = 0 (see the description of JOBL). */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the matrices Q and/or R (if */
/*             JOBB = 'B') are factored, as follows: */
/*             = 'N':  Not factored, Q and R are given; */
/*             = 'C':  C is given, and Q = C'C; */
/*             = 'D':  D is given, and R = D'D (if TYPE = 'O'), or */
/*                     R = D + D' (if TYPE = 'S'); */
/*             = 'B':  Both factors C and D are given, Q = C'C, R = D'D */
/*                     (or R = D + D'). */

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
/*             Using SLICOT Library routine SB02MT to compute the */
/*             corresponding A and Q in this case, before calling SB02OY, */
/*             enables to obtain 2N-by-2N matrix pairs directly. */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether or not the matrix E is identity, as */
/*             follows: */
/*             = 'I':  E is the identity matrix; */
/*             = 'N':  E is a general matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, Q, and E, and the number */
/*             of rows of the matrices B and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             If JOBB = 'B', M is the order of the matrix R, and the */
/*             number of columns of the matrix B.  M >= 0. */
/*             M is not used if JOBB = 'G'. */

/*     P       (input) INTEGER */
/*             If FACT = 'C' or 'D' or 'B', or if TYPE = 'S', P is the */
/*             number of rows of the matrix C and/or D, respectively. */
/*             P >= 0, and if JOBB = 'B' and TYPE = 'S', then P = M. */
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
/*             output weighting matrix Q. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             If FACT = 'C' or 'B', the leading P-by-N part of this */
/*             array must contain the output matrix C of the system. */

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
/*             If FACT = 'D' or 'B', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D of the */
/*             system. */
/*             If JOBB = 'G', this array is not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R. */
/*             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C'; */
/*             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B'; */
/*             LDR >= 1        if JOBB = 'G'. */

/*     L       (input) DOUBLE PRECISION array, dimension (LDL,M) */
/*             If JOBL = 'N' (and JOBB = 'B'), the leading N-by-M part of */
/*             this array must contain the cross weighting matrix L. */
/*             If JOBL = 'Z' or JOBB = 'G', this array is not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N'; */
/*             LDL >= 1        if JOBL = 'Z' or JOBB = 'G'. */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             If JOBE = 'N', the leading N-by-N part of this array must */
/*             contain the matrix E of the descriptor system. */
/*             If JOBE = 'I', E is taken as identity and this array is */
/*             not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of array E. */
/*             LDE >= MAX(1,N) if JOBE = 'N'; */
/*             LDE >= 1        if JOBE = 'I'. */

/*     AF      (output) DOUBLE PRECISION array, dimension (LDAF,*) */
/*             The leading 2N-by-2N part of this array contains the */
/*             matrix A  in the matrix pencil. */
/*                     f */
/*             Array AF must have 2*N+M columns if JOBB = 'B', and 2*N */
/*             columns, otherwise. */

/*     LDAF    INTEGER */
/*             The leading dimension of array AF. */
/*             LDAF >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDAF >= MAX(1,2*N)   if JOBB = 'G'. */

/*     BF      (output) DOUBLE PRECISION array, dimension (LDBF,2*N) */
/*             If DICO = 'D' or JOBB = 'B' or JOBE = 'N', the leading */
/*             2N-by-2N part of this array contains the matrix B  in the */
/*                                                              f */
/*             matrix pencil. */
/*             The last M zero columns are never constructed. */
/*             If DICO = 'C' and JOBB = 'G' and JOBE = 'I', this array */
/*             is not referenced. */

/*     LDBF    INTEGER */
/*             The leading dimension of array BF. */
/*             LDBF >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDBF >= MAX(1,2*N)   if JOBB = 'G' and ( DICO = 'D' or */
/*                                                      JOBE = 'N' ), */
/*             LDBF >= 1            if JOBB = 'G' and ( DICO = 'C' and */
/*                                                      JOBE = 'I' ). */

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
/*             LIWORK >= M if JOBB = 'B', */
/*             LIWORK >= 1 if JOBB = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. If JOBB = 'B', DWORK(2) returns the reciprocal */
/*             of the condition number of the M-by-M lower triangular */
/*             matrix obtained after compression. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1                  if JOBB = 'G', */
/*             LDWORK >= MAX(1,2*N + M,3*M) if JOBB = 'B'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the computed extended matrix pencil is singular, */
/*                   possibly due to rounding errors. */

/*     METHOD */

/*     The extended matrix pairs are constructed, taking various options */
/*     into account. If JOBB = 'B', the problem order is reduced from */
/*     2N+M to 2N (see [1]). */

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

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB02CY by T.G.J. Beelen, Philips, */
/*     Eindhoven, Holland, M. Vanbegin, and P. Van Dooren, Philips */
/*     Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

/*     ****************************************************************** */

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
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    l_dim1 = *ldl;
    l_offset = 1 + l_dim1;
    l -= l_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    bf_dim1 = *ldbf;
    bf_offset = 1 + bf_dim1;
    bf -= bf_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    optc = lsame_(type__, "O", (ftnlen)1, (ftnlen)1);
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    ljobb = lsame_(jobb, "B", (ftnlen)1, (ftnlen)1);
    lfacn = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
    lfacq = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
    lfacr = lsame_(fact, "D", (ftnlen)1, (ftnlen)1);
    lfacb = lsame_(fact, "B", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    ljobe = lsame_(jobe, "I", (ftnlen)1, (ftnlen)1);
    n2 = *n + *n;
    if (ljobb) {
	ljobl = lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1);
	nm = *n + *m;
	nnm = n2 + *m;
    } else {
	nm = *n;
	nnm = n2;
    }
    np1 = *n + 1;
    n2p1 = n2 + 1;

/*     Test the input scalar arguments. */

    if (! optc && ! lsame_(type__, "S", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! ljobb && ! lsame_(jobb, "G", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (! lfacq && ! lfacr && ! lfacb && ! lfacn) {
	*info = -4;
    } else if (! ljobb || lfacn) {
	if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	    *info = -5;
	}
    } else if (ljobb) {
	if (! ljobl && ! lsame_(jobl, "N", (ftnlen)1, (ftnlen)1)) {
	    *info = -6;
	}
    } else if (! ljobe && ! lsame_(jobe, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -7;
    } else if (*n < 0) {
	*info = -8;
    } else if (ljobb) {
	if (*m < 0) {
	    *info = -9;
	}
    } else if (! lfacn || ! optc) {
	if (*p < 0) {
	    *info = -10;
	} else if (ljobb) {
	    if (! optc && *p != *m) {
		*info = -10;
	    }
	}
    } else if (*lda < max(1,*n)) {
	*info = -12;
    } else if (*ldb < max(1,*n)) {
	*info = -14;
    } else if ((lfacn || lfacr) && *ldq < max(1,*n) || (lfacq || lfacb) && *
	    ldq < max(1,*p)) {
	*info = -16;
    } else if (*ldr < 1) {
	*info = -18;
    } else if (ljobb) {
	if ((lfacn || lfacq) && *ldr < *m || (lfacr || lfacb) && *ldr < *p) {
	    *info = -18;
	} else if (! ljobl && *ldl < max(1,*n) || ljobl && *ldl < 1) {
	    *info = -20;
	}
    }
    if (! ljobe && *lde < max(1,*n) || ljobe && *lde < 1) {
	*info = -22;
    } else if (*ldaf < max(1,nnm)) {
	*info = -24;
    } else if ((ljobb || discr || ! ljobe) && *ldbf < nnm || *ldbf < 1) {
	*info = -26;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = nnm, i__2 = *m * 3;
	if (ljobb && *ldwork < max(i__1,i__2) || *ldwork < 1) {
	    *info = -30;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB02OY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    dwork[1] = 1.;
    if (*n == 0) {
	return 0;
    }

/*     Construct the extended matrices in AF and BF, by block-columns. */

    dlacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (ftnlen)4);

    if (! lfacq && ! lfacb) {
	dlacpy_(uplo, n, n, &q[q_offset], ldq, &af[np1 + af_dim1], ldaf, (
		ftnlen)1);
	if (luplo) {

/*           Construct the lower triangle of Q. */

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		dcopy_(&i__2, &q[j + (j + 1) * q_dim1], ldq, &af[np1 + j + j *
			 af_dim1], &c__1);
/* L20: */
	    }

	} else {

/*           Construct the upper triangle of Q. */

	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy_(&i__2, &q[j + q_dim1], ldq, &af[np1 + j * af_dim1], &
			c__1);
/* L40: */
	    }

	}
    } else {
	dsyrk_("Upper", "Transpose", n, p, &c_b26, &q[q_offset], ldq, &c_b27, 
		&af[np1 + af_dim1], ldaf, (ftnlen)5, (ftnlen)9);

	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    dcopy_(&i__2, &af[np1 + j * af_dim1], &c__1, &af[*n + j + af_dim1]
		    , ldaf);
/* L60: */
	}

    }

    if (ljobb) {
	if (ljobl) {
	    dlaset_("Full", m, n, &c_b27, &c_b27, &af[n2p1 + af_dim1], ldaf, (
		    ftnlen)4);
	} else {

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(m, &l[i__ + l_dim1], ldl, &af[n2p1 + i__ * af_dim1], &
			c__1);
/* L80: */
	    }

	}
    }

    if (discr || ljobb) {
	dlaset_("Full", n, n, &c_b27, &c_b27, &af[np1 * af_dim1 + 1], ldaf, (
		ftnlen)4);
    } else {
	if (luplo) {

/*           Construct (1,2) block of AF using the upper triangle of G. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    af[i__ + (*n + j) * af_dim1] = -b[i__ + j * b_dim1];
/* L100: */
		}

		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    af[i__ + (*n + j) * af_dim1] = -b[j + i__ * b_dim1];
/* L120: */
		}

/* L140: */
	    }

	} else {

/*           Construct (1,2) block of AF using the lower triangle of G. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    af[i__ + (*n + j) * af_dim1] = -b[j + i__ * b_dim1];
/* L160: */
		}

		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    af[i__ + (*n + j) * af_dim1] = -b[i__ + j * b_dim1];
/* L180: */
		}

/* L200: */
	    }

	}
    }

    if (discr) {
	if (ljobe) {
	    dlaset_("Full", &nm, n, &c_b27, &c_b46, &af[np1 + np1 * af_dim1], 
		    ldaf, (ftnlen)4);
	} else {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    af[*n + i__ + (*n + j) * af_dim1] = -e[j + i__ * e_dim1];
/* L220: */
		}

/* L240: */
	    }

	    if (ljobb) {
		dlaset_("Full", m, n, &c_b27, &c_b27, &af[n2p1 + np1 * 
			af_dim1], ldaf, (ftnlen)4);
	    }
	}
    } else {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		af[*n + i__ + (*n + j) * af_dim1] = a[j + i__ * a_dim1];
/* L260: */
	    }

/* L280: */
	}

	if (ljobb) {
	    if (optc) {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &b[j + b_dim1], ldb, &af[n2p1 + (*n + j) * 
			    af_dim1], &c__1);
/* L300: */
		}

	    } else {
		dlacpy_("Full", p, n, &q[q_offset], ldq, &af[n2p1 + np1 * 
			af_dim1], ldaf, (ftnlen)4);
	    }
	}
    }

    if (ljobb) {

	if (optc) {
	    dlacpy_("Full", n, m, &b[b_offset], ldb, &af[n2p1 * af_dim1 + 1], 
		    ldaf, (ftnlen)4);
	} else {

	    i__1 = *p;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &q[i__ + q_dim1], ldq, &af[(n2 + i__) * af_dim1 + 1]
			, &c__1);
/* L320: */
	    }

	}

	if (ljobl) {
	    dlaset_("Full", n, m, &c_b27, &c_b27, &af[np1 + n2p1 * af_dim1], 
		    ldaf, (ftnlen)4);
	} else {
	    dlacpy_("Full", n, m, &l[l_offset], ldl, &af[np1 + n2p1 * af_dim1]
		    , ldaf, (ftnlen)4);
	}

	if (! lfacr && ! lfacb) {
	    dlacpy_(uplo, m, m, &r__[r_offset], ldr, &af[n2p1 + n2p1 * 
		    af_dim1], ldaf, (ftnlen)1);
	    if (luplo) {

/*              Construct the lower triangle of R. */

		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m - j;
		    dcopy_(&i__2, &r__[j + (j + 1) * r_dim1], ldr, &af[n2p1 + 
			    j + (n2 + j) * af_dim1], &c__1);
/* L340: */
		}

	    } else {

/*              Construct the upper triangle of R. */

		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = j - 1;
		    dcopy_(&i__2, &r__[j + r_dim1], ldr, &af[n2p1 + (n2 + j) *
			     af_dim1], &c__1);
/* L360: */
		}

	    }
	} else if (optc) {
	    dsyrk_("Upper", "Transpose", m, p, &c_b26, &r__[r_offset], ldr, &
		    c_b27, &af[n2p1 + n2p1 * af_dim1], ldaf, (ftnlen)5, (
		    ftnlen)9);

	    i__1 = *m;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy_(&i__2, &af[n2p1 + (n2 + j) * af_dim1], &c__1, &af[n2 + 
			j + n2p1 * af_dim1], ldaf);
/* L380: */
	    }

	} else {

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    af[n2 + i__ + (n2 + j) * af_dim1] = r__[i__ + j * r_dim1] 
			    + r__[j + i__ * r_dim1];
/* L400: */
		}

/* L420: */
	    }

	}
    }

    if (! ljobb && ! discr && ljobe) {
	return 0;
    }

/*     Construct the first two block columns of BF. */

    if (ljobe) {
	i__1 = *n + nm;
	dlaset_("Full", &i__1, n, &c_b27, &c_b26, &bf[bf_offset], ldbf, (
		ftnlen)4);
    } else {
	dlacpy_("Full", n, n, &e[e_offset], lde, &bf[bf_offset], ldbf, (
		ftnlen)4);
	dlaset_("Full", &nm, n, &c_b27, &c_b27, &bf[np1 + bf_dim1], ldbf, (
		ftnlen)4);
    }

    if (! discr || ljobb) {
	dlaset_("Full", n, n, &c_b27, &c_b27, &bf[np1 * bf_dim1 + 1], ldbf, (
		ftnlen)4);
    } else {
	if (luplo) {

/*           Construct (1,2) block of BF using the upper triangle of G. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    bf[i__ + (*n + j) * bf_dim1] = b[i__ + j * b_dim1];
/* L440: */
		}

		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    bf[i__ + (*n + j) * bf_dim1] = b[j + i__ * b_dim1];
/* L460: */
		}

/* L480: */
	    }

	} else {

/*           Construct (1,2) block of BF using the lower triangle of G. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    bf[i__ + (*n + j) * bf_dim1] = b[j + i__ * b_dim1];
/* L500: */
		}

		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    bf[i__ + (*n + j) * bf_dim1] = b[i__ + j * b_dim1];
/* L520: */
		}

/* L540: */
	    }

	}
    }

    if (discr) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		bf[*n + i__ + (*n + j) * bf_dim1] = -a[j + i__ * a_dim1];
/* L560: */
	    }

/* L580: */
	}

	if (ljobb) {

	    if (optc) {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			bf[n2 + i__ + (*n + j) * bf_dim1] = -b[j + i__ * 
				b_dim1];
/* L600: */
		    }

/* L620: */
		}

	    } else {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

		    i__2 = *p;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			bf[n2 + i__ + (*n + j) * bf_dim1] = -q[i__ + j * 
				q_dim1];
/* L640: */
		    }

/* L660: */
		}

	    }
	}

    } else {
	if (ljobe) {
	    dlaset_("Full", &nm, n, &c_b27, &c_b46, &bf[np1 + np1 * bf_dim1], 
		    ldbf, (ftnlen)4);
	} else {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    bf[*n + i__ + (*n + j) * bf_dim1] = -e[j + i__ * e_dim1];
/* L680: */
		}

/* L700: */
	    }

	    if (ljobb) {
		dlaset_("Full", m, n, &c_b27, &c_b27, &bf[n2p1 + np1 * 
			bf_dim1], ldbf, (ftnlen)4);
	    }
	}
    }

    if (! ljobb) {
	return 0;
    }

/*     Compress the pencil lambda x BF - AF, using QL factorization. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/*     Workspace: need 2*M;  prefer M + M*NB. */

    itau = 1;
    jwork = itau + *m;
    i__1 = *ldwork - jwork + 1;
    dgeqlf_(&nnm, m, &af[n2p1 * af_dim1 + 1], ldaf, &dwork[itau], &dwork[
	    jwork], &i__1, info);
    wrkopt = (integer) dwork[jwork];

/*     Workspace: need 2*N+M;  prefer M + 2*N*NB. */

    i__1 = *ldwork - jwork + 1;
    dormql_("Left", "Transpose", &nnm, &n2, m, &af[n2p1 * af_dim1 + 1], ldaf, 
	    &dwork[itau], &af[af_offset], ldaf, &dwork[jwork], &i__1, info, (
	    ftnlen)4, (ftnlen)9);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

    i__1 = *ldwork - jwork + 1;
    dormql_("Left", "Transpose", &nnm, &n2, m, &af[n2p1 * af_dim1 + 1], ldaf, 
	    &dwork[itau], &bf[bf_offset], ldbf, &dwork[jwork], &i__1, info, (
	    ftnlen)4, (ftnlen)9);

/*     Check the singularity of the L factor in the QL factorization: */
/*     if singular, then the extended matrix pencil is also singular. */
/*     Workspace 3*M. */

    toldef = *tol;
    if (toldef <= 0.) {
	toldef = dlamch_("Epsilon", (ftnlen)7);
    }

    dtrcon_("1-norm", "Lower", "Non unit", m, &af[n2p1 + n2p1 * af_dim1], 
	    ldaf, &rcond, &dwork[1], &iwork[1], info, (ftnlen)6, (ftnlen)5, (
	    ftnlen)8);
/* Computing MAX */
    i__1 = wrkopt, i__2 = *m * 3;
    wrkopt = max(i__1,i__2);

    if (rcond <= toldef) {
	*info = 1;
    }

    dwork[1] = (doublereal) wrkopt;
    dwork[2] = rcond;

    return 0;
/* *** Last line of SB02OY *** */
} /* sb02oy_ */

