/* NF01BR.f -- translated by f2c (version 20100827).
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
static doublereal c_b19 = 0.;
static integer c__0 = 0;
static doublereal c_b56 = -1.;
static doublereal c_b58 = 1.;

/* Subroutine */ int nf01br_(char *cond, char *uplo, char *trans, integer *n, 
	integer *ipar, integer *lipar, doublereal *r__, integer *ldr, 
	doublereal *sdiag, doublereal *s, integer *lds, doublereal *b, 
	integer *ranks, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen cond_len, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, s_dim1, s_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, i1, bn, nc, st, bsm, bsn;
    static doublereal dum[3];
    static integer rank;
    static logical full;
    static integer nths;
    extern /* Subroutine */ int mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd, ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    static logical lower, tranr;
    static char uplol[1];
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static char transl[1];


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

/*     To solve one of the systems of linear equations */

/*           R*x = b ,  or  R'*x = b , */

/*     in the least squares sense, where R is an n-by-n block upper */
/*     triangular matrix, with the structure */

/*         /   R_1    0    ..   0   |   L_1   \ */
/*         |    0    R_2   ..   0   |   L_2   | */
/*         |    :     :    ..   :   |    :    | , */
/*         |    0     0    ..  R_l  |   L_l   | */
/*         \    0     0    ..   0   |  R_l+1  / */

/*     with the upper triangular submatrices R_k, k = 1:l+1, square, and */
/*     the first l of the same order, BSN. The diagonal elements of each */
/*     block R_k have nonincreasing magnitude. The matrix R is stored in */
/*     the compressed form, as returned by SLICOT Library routine NF01BS, */

/*              /   R_1  |   L_1   \ */
/*              |   R_2  |   L_2   | */
/*       Rc =   |    :   |    :    | , */
/*              |   R_l  |   L_l   | */
/*              \    X   |  R_l+1  / */

/*     where the submatrix X is irrelevant. If the matrix R does not have */
/*     full rank, then a least squares solution is obtained. If l <= 1, */
/*     then R is an upper triangular matrix and its full upper triangle */
/*     is stored. */

/*     Optionally, the transpose of the matrix R can be stored in the */
/*     strict lower triangles of the submatrices R_k, k = 1:l+1, and in */
/*     the arrays SDIAG and S, as described at the parameter UPLO below. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of submatrices R_k should */
/*             be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation and store */
/*                      the numerical rank of R_k in the array entry */
/*                      RANKS(k), for k = 1:l+1; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R_k for zero values; */
/*             = 'U' :  use the ranks already stored in RANKS(1:l+1). */

/*     UPLO    CHARACTER*1 */
/*             Specifies the storage scheme for the matrix R, as follows: */
/*             = 'U' :  the upper triangular part is stored as in Rc; */
/*             = 'L' :  the lower triangular part is stored, namely, */
/*                      - the transpose of the strict upper triangle of */
/*                        R_k is stored in the strict lower triangle of */
/*                        R_k, for k = 1:l+1; */
/*                      - the diagonal elements of R_k, k = 1:l+1, are */
/*                        stored in the array SDIAG; */
/*                      - the transpose of the last block column in R */
/*                        (without R_l+1) is stored in the array S. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of the system of equations, as follows: */
/*             = 'N':  R*x  = b  (No transpose); */
/*             = 'T':  R'*x = b  (Transpose); */
/*             = 'C':  R'*x = b  (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N = BN*BSN + ST >= 0. */
/*             (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix R, as follows: */
/*             IPAR(1) must contain ST, the number of columns of the */
/*                     submatrices L_k and the order of R_l+1.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, l, in the */
/*                     block diagonal part of R.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     R_k, k = 1:l.  BSM >= 0. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks R_k, k = 1:l.  BSN >= 0. */
/*             BSM is not used by this routine, but assumed equal to BSN. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     R       (input) DOUBLE PRECISION array, dimension (LDR, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             If UPLO = 'U', the leading N-by-NC part of this array must */
/*             contain the (compressed) representation (Rc) of the upper */
/*             triangular matrix R. The submatrix X in Rc and the strict */
/*             lower triangular parts of the diagonal blocks R_k, */
/*             k = 1:l+1, are not referenced. If BN <= 1 or BSN = 0, then */
/*             the full upper triangle of R must be stored. */
/*             If UPLO = 'L', BN > 1 and BSN > 0, the leading */
/*             (N-ST)-by-BSN part of this array must contain the */
/*             transposes of the strict upper triangles of R_k, k = 1:l, */
/*             stored in the strict lower triangles of R_k, and the */
/*             strict lower triangle of R_l+1 must contain the transpose */
/*             of the strict upper triangle of R_l+1. The submatrix X */
/*             in Rc is not referenced. The diagonal elements of R_k, */
/*             and, if COND = 'E', the upper triangular parts of R_k, */
/*             k = 1:l+1, are modified internally, but are restored */
/*             on exit. */
/*             If UPLO = 'L' and BN <= 1 or BSN = 0, the leading N-by-N */
/*             strict lower triangular part of this array must contain */
/*             the transpose of the strict upper triangular part of R. */
/*             The diagonal elements and, if COND = 'E', the upper */
/*             triangular elements are modified internally, but are */
/*             restored on exit. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R.  LDR >= MAX(1,N). */

/*     SDIAG   (input) DOUBLE PRECISION array, dimension (N) */
/*             If UPLO = 'L', this array must contain the diagonal */
/*             entries of R_k, k = 1:l+1. This array is modified */
/*             internally, but is restored on exit. */
/*             This parameter is not referenced if UPLO = 'U'. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N-ST) */
/*             If UPLO = 'L', BN > 1, and BSN > 0, the leading */
/*             ST-by-(N-ST) part of this array must contain the transpose */
/*             of the rectangular part of the last block column in R, */
/*             that is [ L_1' L_2' ... L_l' ] . If COND = 'E', S is */
/*             modified internally, but is restored on exit. */
/*             This parameter is not referenced if UPLO = 'U', or */
/*             BN <= 1, or BSN = 0. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= 1,         if UPLO = 'U', or BN <= 1, or BSN = 0; */
/*             LDS >= MAX(1,ST), if UPLO = 'L', BN > 1, and BSN > 0. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the right hand side */
/*             vector b. */
/*             On exit, this array contains the (least squares) solution */
/*             of the system R*x = b or R'*x = b. */

/*     RANKS   (input or output) INTEGER array, dimension (r), where */
/*             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1; */
/*             r = BN,      if ST = 0 and BSN > 0; */
/*             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 ); */
/*             r = 0,       if ST = 0 and BSN = 0. */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical ranks of the submatrices R_k, k = 1:l(+1). */
/*             On exit, if COND = 'E' or 'N' and N > 0, this array */
/*             contains the numerical ranks of the submatrices R_k, */
/*             k = 1:l(+1), estimated according to the value of COND. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the submatrices R_k. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank. If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Denote Full = ( BN <= 1 or  BSN = 0 ); */
/*                    Comp = ( BN >  1 and BSN > 0 ). */
/*             LDWORK >= 2*N,           if Full and COND = 'E'; */
/*             LDWORK >= 2*MAX(BSN,ST), if Comp and COND = 'E'; */
/*             LDWORK >= 0,   in the remaining cases. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Block back or forward substitution is used (depending on TRANS */
/*     and UPLO), exploiting the special structure and storage scheme of */
/*     the matrix R. If a submatrix R_k, k = 1:l+1, is singular, a local */
/*     basic least squares solution is computed. Therefore, the returned */
/*     result is not the basic least squares solution for the whole */
/*     problem, but a concatenation of (least squares) solutions of the */
/*     individual subproblems involving R_k, k = 1:l+1 (with adapted */
/*     right hand sides). */

/*     NUMERICAL ASPECTS */
/*                                    2    2 */
/*     The algorithm requires 0(BN*BSN + ST + N*ST) operations and is */
/*     backward stable, if R is nonsingular. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --ipar;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --sdiag;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    --b;
    --ranks;
    --dwork;

    /* Function Body */
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
    tranr = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

    *info = 0;
    if (! (econd || ncond || lsame_(cond, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (tranr || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lipar < 4) {
	*info = -6;
    } else {
	st = ipar[1];
	bn = ipar[2];
	bsm = ipar[3];
	bsn = ipar[4];
	nths = bn * bsn;
	full = bn <= 1 || bsn == 0;
/* Computing MIN */
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
	if (min(i__1,bsn) < 0) {
	    *info = -5;
	} else if (*n != nths + st) {
	    *info = -4;
	} else if (*ldr < max(1,*n)) {
	    *info = -8;
	} else if (*lds < 1 || lower && ! full && *lds < st) {
	    *info = -11;
	} else {
	    if (econd) {
		if (full) {
		    l = *n << 1;
		} else {
		    l = max(bsn,st) << 1;
		}
	    } else {
		l = 0;
	    }
	    if (*ldwork < l) {
		*info = -16;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BR", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    if (econd) {
	toldef = *tol;
	if (toldef <= 0.) {

/*           Use the default tolerance in rank determination. */

	    toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	}
    }

    nc = bsn + st;
    if (full) {

/*        Special case: l <= 1 or BSN = 0; R is just an upper triangular */
/*        matrix. */

	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and, if COND = 'E', swap the upper and lower triangular */
/*           parts of R, in order to find the numerical rank. */

	    i__1 = *ldr + 1;
	    dswap_(n, &r__[r_offset], &i__1, &sdiag[1], &c__1);
	    if (econd) {
		*(unsigned char *)uplol = 'U';
		*(unsigned char *)transl = *(unsigned char *)trans;

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - j + 1;
		    dswap_(&i__2, &r__[j + j * r_dim1], ldr, &r__[j + j * 
			    r_dim1], &c__1);
/* L10: */
		}

	    } else {
		*(unsigned char *)uplol = *(unsigned char *)uplo;
		if (tranr) {
		    *(unsigned char *)transl = 'N';
		} else {
		    *(unsigned char *)transl = 'T';
		}
	    }
	} else {
	    *(unsigned char *)uplol = *(unsigned char *)uplo;
	    *(unsigned char *)transl = *(unsigned char *)trans;
	}

	if (econd) {

/*           Estimate the reciprocal condition number and set the rank. */
/*           Workspace: 2*N. */

	    mb03od_("No QR", n, n, &r__[r_offset], ldr, &ipar[1], &toldef, &
		    c_b19, &dwork[1], &rank, dum, &dwork[1], ldwork, info, (
		    ftnlen)5);
	    ranks[1] = rank;

	} else if (ncond) {

/*           Determine rank(R) by checking zero diagonal entries. */

	    rank = *n;

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (r__[j + j * r_dim1] == 0. && rank == *n) {
		    rank = j - 1;
		}
/* L20: */
	    }

	    ranks[1] = rank;

	} else {

/*           Use the stored rank. */

	    rank = ranks[1];
	}

/*        Solve R*x = b, or R'*x = b using back or forward substitution. */

	dum[0] = 0.;
	if (rank < *n) {
	    i__1 = *n - rank;
	    dcopy_(&i__1, dum, &c__0, &b[rank + 1], &c__1);
	}
	dtrsv_(uplol, transl, "NonUnit", &rank, &r__[r_offset], ldr, &b[1], &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)7);

	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and, if COND = 'E', swap back the upper and lower triangular */
/*           parts of R. */

	    i__1 = *ldr + 1;
	    dswap_(n, &r__[r_offset], &i__1, &sdiag[1], &c__1);
	    if (econd) {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - j + 1;
		    dswap_(&i__2, &r__[j + j * r_dim1], ldr, &r__[j + j * 
			    r_dim1], &c__1);
/* L30: */
		}

	    }

	}
	return 0;
    }

/*     General case: l > 1 and BSN > 0. */

    i__ = 1;
    l = bn;
    if (econd) {

/*        Estimate the reciprocal condition numbers and set the ranks. */

	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and swap the upper and lower triangular parts of R, in order */
/*           to find the numerical rank. Swap S and the transpose of the */
/*           rectangular part of the last block column of R. */

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *ldr + 1;
		dswap_(&bsn, &r__[i__ + r_dim1], &i__2, &sdiag[i__], &c__1);

		i__2 = bsn;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = bsn - j + 1;
		    dswap_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			    r_dim1], &c__1);
		    ++i__;
/* L40: */
		}

/* L50: */
	    }

	    if (st > 0) {
		i__1 = *ldr + 1;
		dswap_(&st, &r__[i__ + (bsn + 1) * r_dim1], &i__1, &sdiag[i__]
			, &c__1);

		i__1 = nc;
		for (j = bsn + 1; j <= i__1; ++j) {
		    dswap_(&nths, &r__[j * r_dim1 + 1], &c__1, &s[j - bsn + 
			    s_dim1], lds);
		    i__2 = nc - j + 1;
		    dswap_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			    r_dim1], &c__1);
		    ++i__;
/* L60: */
		}

	    }

	}

	i1 = 1;

/*        Determine rank(R_k) using incremental condition estimation. */
/*        Workspace 2*MAX(BSN,ST). */

	i__1 = bn;
	for (k = 1; k <= i__1; ++k) {
	    mb03od_("No QR", &bsn, &bsn, &r__[i1 + r_dim1], ldr, &ipar[1], &
		    toldef, &c_b19, &dwork[1], &ranks[k], dum, &dwork[1], 
		    ldwork, info, (ftnlen)5);
	    i1 += bsn;
/* L70: */
	}

	if (st > 0) {
	    ++l;
	    mb03od_("No QR", &st, &st, &r__[i1 + (bsn + 1) * r_dim1], ldr, &
		    ipar[1], &toldef, &c_b19, &dwork[1], &ranks[l], dum, &
		    dwork[1], ldwork, info, (ftnlen)5);
	}

    } else if (ncond) {

/*        Determine rank(R_k) by checking zero diagonal entries. */

	if (lower) {

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		rank = bsn;

		i__2 = bsn;
		for (j = 1; j <= i__2; ++j) {
		    if (sdiag[i__] == 0. && rank == bsn) {
			rank = j - 1;
		    }
		    ++i__;
/* L80: */
		}

		ranks[k] = rank;
/* L90: */
	    }

	    if (st > 0) {
		++l;
		rank = st;

		i__1 = st;
		for (j = 1; j <= i__1; ++j) {
		    if (sdiag[i__] == 0. && rank == st) {
			rank = j - 1;
		    }
		    ++i__;
/* L100: */
		}

		ranks[l] = rank;
	    }

	} else {

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		rank = bsn;

		i__2 = bsn;
		for (j = 1; j <= i__2; ++j) {
		    if (r__[i__ + j * r_dim1] == 0. && rank == bsn) {
			rank = j - 1;
		    }
		    ++i__;
/* L110: */
		}

		ranks[k] = rank;
/* L120: */
	    }

	    if (st > 0) {
		++l;
		rank = st;

		i__1 = nc;
		for (j = bsn + 1; j <= i__1; ++j) {
		    if (r__[i__ + j * r_dim1] == 0. && rank == st) {
			rank = j - bsn - 1;
		    }
		    ++i__;
/* L130: */
		}

		ranks[l] = rank;
	    }
	}

    } else {

/*        Set the number of elements of RANKS. Then use the stored ranks. */

	if (st > 0) {
	    ++l;
	}
    }

/*     Solve the triangular system for x. If the system is singular, */
/*     then obtain a basic least squares solution. */

    dum[0] = 0.;
    if (lower && ! econd) {

	if (! tranr) {

/*           Solve R*x = b using back substitution, with R' stored in */
/*           the arrays R, SDIAG and S. Swap diag(R) and SDIAG. */

	    i1 = nths + 1;
	    if (st > 0) {
		rank = ranks[l];
		if (rank < st) {
		    i__1 = st - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		i__1 = *ldr + 1;
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
		dtrsv_("Lower", "Transpose", "NonUnit", &rank, &r__[i1 + (bsn 
			+ 1) * r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (
			ftnlen)9, (ftnlen)7);
		i__1 = *ldr + 1;
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
		dgemv_("Transpose", &st, &nths, &c_b56, &s[s_offset], lds, &b[
			nths + 1], &c__1, &c_b58, &b[1], &c__1, (ftnlen)9);
	    }

	    for (k = bn; k >= 1; --k) {
		i1 -= bsn;
		rank = ranks[k];
		if (rank < bsn) {
		    i__1 = bsn - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		i__1 = *ldr + 1;
		dswap_(&bsn, &r__[i1 + r_dim1], &i__1, &sdiag[i1], &c__1);
		dtrsv_("Lower", "Transpose", "NonUnit", &rank, &r__[i1 + 
			r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)9, (
			ftnlen)7);
		i__1 = *ldr + 1;
		dswap_(&bsn, &r__[i1 + r_dim1], &i__1, &sdiag[i1], &c__1);
/* L140: */
	    }

	} else {

/*           Solve R'*x = b using forward substitution, with R' stored in */
/*           the arrays R, SDIAG and S. Swap diag(R) and SDIAG. */

	    i1 = 1;
	    if (tranr) {
		*(unsigned char *)transl = 'N';
	    } else {
		*(unsigned char *)transl = 'T';
	    }

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		rank = ranks[k];
		if (rank < bsn) {
		    i__2 = bsn - rank;
		    dcopy_(&i__2, dum, &c__0, &b[i1 + rank], &c__1);
		}
		i__2 = *ldr + 1;
		dswap_(&bsn, &r__[i1 + r_dim1], &i__2, &sdiag[i1], &c__1);
		dtrsv_("Lower", transl, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
		i__2 = *ldr + 1;
		dswap_(&bsn, &r__[i1 + r_dim1], &i__2, &sdiag[i1], &c__1);
		i1 += bsn;
/* L150: */
	    }

	    if (st > 0) {
		rank = ranks[l];
		if (rank < st) {
		    i__1 = st - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		dgemv_("NoTranspose", &st, &nths, &c_b56, &s[s_offset], lds, &
			b[1], &c__1, &c_b58, &b[i1], &c__1, (ftnlen)11);
		i__1 = *ldr + 1;
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
		dtrsv_("Lower", transl, "NonUnit", &rank, &r__[i1 + (bsn + 1) 
			* r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
		i__1 = *ldr + 1;
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
	    }

	}

    } else {

	if (! tranr) {

/*           Solve R*x = b using back substitution. */

	    i1 = nths + 1;
	    if (st > 0) {
		rank = ranks[l];
		if (rank < st) {
		    i__1 = st - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + (bsn + 1) *
			 r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
		dgemv_(trans, &nths, &st, &c_b56, &r__[(bsn + 1) * r_dim1 + 1]
			, ldr, &b[nths + 1], &c__1, &c_b58, &b[1], &c__1, (
			ftnlen)1);
	    }

	    for (k = bn; k >= 1; --k) {
		i1 -= bsn;
		rank = ranks[k];
		if (rank < bsn) {
		    i__1 = bsn - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
/* L160: */
	    }

	} else {

/*           Solve R'*x = b using forward substitution. */

	    i1 = 1;

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		rank = ranks[k];
		if (rank < bsn) {
		    i__2 = bsn - rank;
		    dcopy_(&i__2, dum, &c__0, &b[i1 + rank], &c__1);
		}
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
		i1 += bsn;
/* L170: */
	    }

	    if (st > 0) {
		rank = ranks[l];
		if (rank < st) {
		    i__1 = st - rank;
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
		}
		dgemv_(trans, &nths, &st, &c_b56, &r__[(bsn + 1) * r_dim1 + 1]
			, ldr, &b[1], &c__1, &c_b58, &b[i1], &c__1, (ftnlen)1)
			;
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + (bsn + 1) *
			 r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
	    }

	}
    }

    if (econd && lower) {
	i__ = 1;

/*        If COND = 'E' and UPLO = 'L', swap the diagonal elements of R */
/*        and the elements of SDIAG and swap back the upper and lower */
/*        triangular parts of R, including the part corresponding to S. */

	i__1 = bn;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *ldr + 1;
	    dswap_(&bsn, &r__[i__ + r_dim1], &i__2, &sdiag[i__], &c__1);

	    i__2 = bsn;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = bsn - j + 1;
		dswap_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
		++i__;
/* L180: */
	    }

/* L190: */
	}

	if (st > 0) {
	    i__1 = *ldr + 1;
	    dswap_(&st, &r__[i__ + (bsn + 1) * r_dim1], &i__1, &sdiag[i__], &
		    c__1);

	    i__1 = nc;
	    for (j = bsn + 1; j <= i__1; ++j) {
		dswap_(&nths, &r__[j * r_dim1 + 1], &c__1, &s[j - bsn + 
			s_dim1], lds);
		i__2 = nc - j + 1;
		dswap_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
		++i__;
/* L200: */
	    }

	}

    }

    return 0;

/* *** Last line of NF01BR *** */
} /* nf01br_ */

