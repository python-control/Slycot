/* SB02MU.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static doublereal c_b61 = 1.;

/* Subroutine */ int sb02mu_(char *dico, char *hinv, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	q, integer *ldq, doublereal *s, integer *lds, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen dico_len, 
	ftnlen hinv_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, s_dim1, 
	    s_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, nj, np1;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond, anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lhinv, luplo;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgetri_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *), dgetrs_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer maxwrk;


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

/*     To construct the 2n-by-2n Hamiltonian or symplectic matrix S */
/*     associated to the linear-quadratic optimization problem, used to */
/*     solve the continuous- or discrete-time algebraic Riccati equation, */
/*     respectively. */

/*     For a continuous-time problem, S is defined by */

/*             (  A  -G ) */
/*         S = (        ),                                       (1) */
/*             ( -Q  -A') */

/*     and for a discrete-time problem by */

/*                 -1       -1 */
/*             (  A        A  *G     ) */
/*         S = (   -1           -1   ),                          (2) */
/*             ( QA     A' + Q*A  *G ) */

/*     or */

/*                       -T         -T */
/*             (  A + G*A  *Q   -G*A   ) */
/*         S = (      -T            -T ),                        (3) */
/*             (    -A  *Q         A   ) */

/*     where A, G, and Q are N-by-N matrices, with G and Q symmetric. */
/*     Matrix A must be nonsingular in the discrete-time case. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  Continuous-time system; */
/*             = 'D':  Discrete-time system. */

/*     HINV    CHARACTER*1 */
/*             If DICO = 'D', specifies which of the matrices (2) or (3) */
/*             is constructed, as follows: */
/*             = 'D':  The matrix S in (2) is constructed; */
/*             = 'I':  The (inverse) matrix S in (3) is constructed. */
/*             HINV is not referenced if DICO = 'C'. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, if DICO = 'D', and INFO = 0, the leading N-by-N */
/*                                                     -1 */
/*             part of this array contains the matrix A  . */
/*             Otherwise, the array A is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             The leading N-by-N upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             must contain the upper triangular part or lower triangular */
/*             part, respectively, of the symmetric matrix G. The stricly */
/*             lower triangular part (if UPLO = 'U') or stricly upper */
/*             triangular part (if UPLO = 'L') is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             The leading N-by-N upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             must contain the upper triangular part or lower triangular */
/*             part, respectively, of the symmetric matrix Q. The stricly */
/*             lower triangular part (if UPLO = 'U') or stricly upper */
/*             triangular part (if UPLO = 'L') is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,N). */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N) */
/*             If INFO = 0, the leading 2N-by-2N part of this array */
/*             contains the Hamiltonian or symplectic matrix of the */
/*             problem. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,2*N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK; if DICO = 'D', DWORK(2) returns the reciprocal */
/*             condition number of the given matrix  A. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1          if DICO = 'C'; */
/*             LDWORK >= MAX(2,4*N) if DICO = 'D'. */
/*             For optimum performance LDWORK should be larger, if */
/*             DICO = 'D'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if the leading i-by-i (1 <= i <= N) upper triangular */
/*                   submatrix of A is singular in discrete-time case; */
/*             = N+1:  if matrix A is numerically singular in discrete- */
/*                   time case. */

/*     METHOD */

/*     For a continuous-time problem, the 2n-by-2n Hamiltonian matrix (1) */
/*     is constructed. */
/*     For a discrete-time problem, the 2n-by-2n symplectic matrix (2) or */
/*     (3) - the inverse of the matrix in (2) - is constructed. */

/*     NUMERICAL ASPECTS */

/*     The discrete-time case needs the inverse of the matrix A, hence */
/*     the routine should not be used when A is ill-conditioned. */
/*                               3 */
/*     The algorithm requires 0(n ) floating point operations in the */
/*     discrete-time case. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

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
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    n2 = *n + *n;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    if (discr) {
	lhinv = lsame_(hinv, "D", (ftnlen)1, (ftnlen)1);
    } else {
	lhinv = FALSE_;
    }

/*     Test the input scalar arguments. */

    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (discr) {
	if (! lhinv && ! lsame_(hinv, "I", (ftnlen)1, (ftnlen)1)) {
	    *info = -2;
	}
    }
    if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldg < max(1,*n)) {
	*info = -8;
    } else if (*ldq < max(1,*n)) {
	*info = -10;
    } else if (*lds < max(1,n2)) {
	*info = -12;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 2, i__2 = *n << 2;
	if (*ldwork < 1 || discr && *ldwork < max(i__1,i__2)) {
	    *info = -15;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB02MU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	if (discr) {
	    dwork[2] = 1.;
	}
	return 0;
    }

/*     The code tries to exploit data locality as much as possible. */

    if (! lhinv) {
	dlacpy_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (ftnlen)4)
		;

/*        Construct Hamiltonian matrix in the continuous-time case, or */
/*        prepare symplectic matrix in (3) in the discrete-time case: */

/*        Construct full Q in S(N+1:2*N,1:N) and change the sign, and */
/*        construct full G in S(1:N,N+1:2*N) and change the sign. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    nj = *n + j;
	    if (luplo) {

		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
/* L20: */
		}

		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
/* L40: */
		}

		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
/* L60: */
		}

		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
/* L80: */
		}

	    } else {

		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
/* L100: */
		}

		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
/* L120: */
		}

		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
/* L140: */
		}

		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
/* L180: */
		}

	    }
/* L200: */
	}

	if (! discr) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		nj = *n + j;

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s[*n + i__ + nj * s_dim1] = -a[j + i__ * a_dim1];
/* L220: */
		}

/* L240: */
	    }

	    dwork[1] = 1.;
	}
    }

    if (discr) {

/*        Construct the symplectic matrix (2) or (3) in the discrete-time */
/*        case. */

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of workspace needed at that point in the code, */
/*        as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

/* Computing MAX */
	i__1 = *n << 2, i__2 = *n * ilaenv_(&c__1, "DGETRI", " ", n, &c_n1, &
		c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
	maxwrk = max(i__1,i__2);
	np1 = *n + 1;

	if (lhinv) {

/*           Put  A'  in  S(N+1:2*N,N+1:2*N). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dcopy_(n, &a[i__ + a_dim1], lda, &s[np1 + (*n + i__) * s_dim1]
			, &c__1);
/* L260: */
	    }

	}

/*        Compute the norm of the matrix A. */

	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);

/*        Compute the LU factorization of A. */

	dgetrf_(n, n, &a[a_offset], lda, &iwork[1], info);

/*        Return if INFO is non-zero. */

	if (*info > 0) {
	    dwork[2] = 0.;
	    return 0;
	}

/*        Compute the reciprocal of the condition number of A. */
/*        Workspace: need 4*N. */

	dgecon_("1-norm", n, &a[a_offset], lda, &anorm, &rcond, &dwork[1], &
		iwork[np1], info, (ftnlen)6);

/*        Return if the matrix is singular to working precision. */

	if (rcond < dlamch_("Epsilon", (ftnlen)7)) {
	    *info = *n + 1;
	    dwork[2] = rcond;
	    return 0;
	}

	if (lhinv) {

/*           Compute S in (2). */

/*           Construct full Q in S(N+1:2*N,1:N). */

	    if (luplo) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(&j, &q[j * q_dim1 + 1], &c__1, &s[np1 + j * s_dim1]
			    , &c__1);
		    i__2 = *n - j;
		    dcopy_(&i__2, &q[j + (j + 1) * q_dim1], ldq, &s[np1 + j + 
			    j * s_dim1], &c__1);
/* L270: */
		}
		dcopy_(n, &q[*n * q_dim1 + 1], &c__1, &s[np1 + *n * s_dim1], &
			c__1);
	    } else {
		dcopy_(n, &q[q_dim1 + 1], &c__1, &s[np1 + s_dim1], &c__1);
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = j - 1;
		    dcopy_(&i__2, &q[j + q_dim1], ldq, &s[np1 + j * s_dim1], &
			    c__1);
		    i__2 = *n - j + 1;
		    dcopy_(&i__2, &q[j + j * q_dim1], &c__1, &s[*n + j + j * 
			    s_dim1], &c__1);
/* L280: */
		}
	    }

/*           Compute the solution matrix  X  of the system  X*A = Q  by */
/*                                                                    -1 */
/*           solving  A'*X' = Q and transposing the result to get  Q*A  . */

	    dgetrs_("Transpose", n, n, &a[a_offset], lda, &iwork[1], &s[np1 + 
		    s_dim1], lds, info, (ftnlen)9);

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		dswap_(&i__2, &s[np1 + j + j * s_dim1], &c__1, &s[*n + j + (j 
			+ 1) * s_dim1], lds);
/* L300: */
	    }

/*           Construct full G in S(1:N,N+1:2*N). */

	    if (luplo) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(&j, &g[j * g_dim1 + 1], &c__1, &s[(*n + j) * 
			    s_dim1 + 1], &c__1);
		    i__2 = *n - j;
		    dcopy_(&i__2, &g[j + (j + 1) * g_dim1], ldg, &s[j + 1 + (*
			    n + j) * s_dim1], &c__1);
/* L310: */
		}
		dcopy_(n, &g[*n * g_dim1 + 1], &c__1, &s[n2 * s_dim1 + 1], &
			c__1);
	    } else {
		dcopy_(n, &g[g_dim1 + 1], &c__1, &s[np1 * s_dim1 + 1], &c__1);
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = j - 1;
		    dcopy_(&i__2, &g[j + g_dim1], ldg, &s[(*n + j) * s_dim1 + 
			    1], &c__1);
		    i__2 = *n - j + 1;
		    dcopy_(&i__2, &g[j + j * g_dim1], &c__1, &s[j + (*n + j) *
			     s_dim1], &c__1);
/* L320: */
		}
	    }
/*                            -1 */
/*           Compute  A' + Q*A  *G  in  S(N+1:2N,N+1:2N). */

	    dgemm_("No transpose", "No transpose", n, n, n, &c_b61, &s[np1 + 
		    s_dim1], lds, &s[np1 * s_dim1 + 1], lds, &c_b61, &s[np1 + 
		    np1 * s_dim1], lds, (ftnlen)12, (ftnlen)12);

/*           Compute the solution matrix  Y  of the system  A*Y = G. */

	    dgetrs_("No transpose", n, n, &a[a_offset], lda, &iwork[1], &s[
		    np1 * s_dim1 + 1], lds, info, (ftnlen)12);

/*           Compute the inverse of  A  in situ. */
/*           Workspace: need N;  prefer N*NB. */

	    dgetri_(n, &a[a_offset], lda, &iwork[1], &dwork[1], ldwork, info);
/*                  -1 */
/*           Copy  A    in  S(1:N,1:N). */

	    dlacpy_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (
		    ftnlen)4);

	} else {

/*           Compute S in (3) using the already prepared part. */

/*           Compute the solution matrix  X'  of the system  A*X' = -G */
/*                                                       -T */
/*           and transpose the result to obtain  X = -G*A  . */

	    dgetrs_("No transpose", n, n, &a[a_offset], lda, &iwork[1], &s[
		    np1 * s_dim1 + 1], lds, info, (ftnlen)12);

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		dswap_(&i__2, &s[j + 1 + (*n + j) * s_dim1], &c__1, &s[j + (
			np1 + j) * s_dim1], lds);
/* L340: */
	    }
/*                           -T */
/*           Compute  A + G*A  *Q  in  S(1:N,1:N). */

	    dgemm_("No transpose", "No transpose", n, n, n, &c_b61, &s[np1 * 
		    s_dim1 + 1], lds, &s[np1 + s_dim1], lds, &c_b61, &s[
		    s_offset], lds, (ftnlen)12, (ftnlen)12);

/*           Compute the solution matrix  Y  of the system  A'*Y = -Q. */

	    dgetrs_("Transpose", n, n, &a[a_offset], lda, &iwork[1], &s[np1 + 
		    s_dim1], lds, info, (ftnlen)9);

/*           Compute the inverse of  A  in situ. */
/*           Workspace: need N;  prefer N*NB. */

	    dgetri_(n, &a[a_offset], lda, &iwork[1], &dwork[1], ldwork, info);
/*                  -T */
/*           Copy  A    in  S(N+1:2N,N+1:2N). */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(n, &a[j + a_dim1], lda, &s[np1 + (*n + j) * s_dim1], &
			c__1);
/* L360: */
	    }

	}
	dwork[1] = (doublereal) maxwrk;
	dwork[2] = rcond;
    }

/* *** Last line of SB02MU *** */
    return 0;
} /* sb02mu_ */

