/* NF01BS.f -- translated by f2c (version 20100827).
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
static logical c_true = TRUE_;

/* Subroutine */ int nf01bs_(integer *n, integer *ipar, integer *lipar, 
	doublereal *fnorm, doublereal *j, integer *ldj, doublereal *e, 
	doublereal *jnorms, doublereal *gnorm, integer *ipvt, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, k, l, m, bn, jl, st, bsm, bsn, jlm, mmn;
    static doublereal sum;
    static integer ibsm, ibsn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau, nths;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int md03bx_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer ibsni;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dlapmt_(logical *, integer *, integer *, 
	    doublereal *, integer *, integer *), dormqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
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

/*     To compute the QR factorization of the Jacobian matrix J, as */
/*     received in compressed form from SLICOT Library routine NF01BD, */

/*            /  dy(1)/dwb(1)  |  dy(1)/ dtheta  \ */
/*       Jc = |       :        |       :         | , */
/*            \  dy(L)/dwb(L)  |  dy(L)/ dtheta  / */

/*     and to apply the transformation Q on the error vector e (in-situ). */
/*     The factorization is J*P = Q*R, where Q is a matrix with */
/*     orthogonal columns, P a permutation matrix, and R an upper */
/*     trapezoidal matrix with diagonal elements of nonincreasing */
/*     magnitude for each block column (see below). The 1-norm of the */
/*     scaled gradient is also returned. */

/*     Actually, the Jacobian J has the block form */

/*       dy(1)/dwb(1)       0         .....       0        dy(1)/dtheta */
/*            0        dy(2)/dwb(2)   .....       0        dy(2)/dtheta */
/*          .....         .....       .....     .....         ..... */
/*            0           .....         0    dy(L)/dwb(L)  dy(L)/dtheta */

/*     but the zero blocks are omitted. The diagonal blocks have the */
/*     same size and correspond to the nonlinear part. The last block */
/*     column corresponds to the linear part. It is assumed that the */
/*     Jacobian matrix has at least as many rows as columns. The linear */
/*     or nonlinear parts can be empty. If L <= 1, the Jacobian is */
/*     represented as a full matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J. */
/*             N = BN*BSN + ST >= 0.  (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain ST, the number of parameters */
/*                     corresponding to the linear part.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, BN = L, */
/*                     for the parameters corresponding to the nonlinear */
/*                     part.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     J_k = dy(k)/dwb(k), k = 1:BN, if BN > 0, or the */
/*                     number of rows of the matrix J, if BN <= 1. */
/*                     BN*BSM >= N, if BN > 0; */
/*                     BSM >= N,    if BN = 0. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     FNORM   (input) DOUBLE PRECISION */
/*             The Euclidean norm of the vector e.  FNORM >= 0. */

/*     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             On entry, the leading NR-by-NC part of this array must */
/*             contain the (compressed) representation (Jc) of the */
/*             Jacobian matrix J, where NR = BSM if BN <= 1, and */
/*             NR = BN*BSM, if BN > 1. */
/*             On exit, the leading N-by-NC part of this array contains */
/*             a (compressed) representation of the upper triangular */
/*             factor R of the Jacobian matrix. The matrix R has the same */
/*             structure as the Jacobian matrix J, but with an additional */
/*             diagonal block. Note that for efficiency of the later */
/*             calculations, the matrix R is delivered with the leading */
/*             dimension MAX(1,N), possibly much smaller than the value */
/*             of LDJ on entry. */

/*     LDJ     (input/output) INTEGER */
/*             The leading dimension of array J. */
/*             On entry, LDJ >= MAX(1,NR). */
/*             On exit,  LDJ >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (NR) */
/*             On entry, this array contains the vector e, */
/*             e = vec( Y - y ), where Y is set of output samples, and */
/*             vec denotes the concatenation of the columns of a matrix. */
/*             On exit, this array contains the updated vector Z*Q'*e, */
/*             where Z is the block row permutation matrix used in the */
/*             QR factorization of J (see METHOD). */

/*     JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the Euclidean norms of the columns */
/*             of the Jacobian matrix, considered in the initial order. */

/*     GNORM   (output) DOUBLE PRECISION */
/*             If FNORM > 0, the 1-norm of the scaled vector J'*e/FNORM, */
/*             with each element i further divided by JNORMS(i) (if */
/*             JNORMS(i) is nonzero). */
/*             If FNORM = 0, the returned value of GNORM is 0. */

/*     IPVT    (output) INTEGER array, dimension (N) */
/*             This array defines the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1,      if N = 0 or BN <= 1 and BSM = N = 1; */
/*                               otherwise, */
/*             LDWORK >= 4*N+1,  if BN <= 1 or  BSN = 0; */
/*             LDWORK >= JWORK,  if BN >  1 and BSN > 0, where JWORK is */
/*                               given by the following procedure: */
/*              JWORK  = BSN + MAX(3*BSN+1,ST); */
/*              JWORK  = MAX(JWORK,4*ST+1),         if BSM > BSN; */
/*              JWORK  = MAX(JWORK,(BSM-BSN)*(BN-1)), */
/*                                                  if BSN < BSM < 2*BSN. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     A QR factorization with column pivoting of the matrix J is */
/*     computed, J*P = Q*R. */

/*     If l = L > 1, the R factor of the QR factorization has the same */
/*     structure as the Jacobian, but with an additional diagonal block. */
/*     Denote */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_l  |  L_l  / */

/*     The algorithm consists in two phases. In the first phase, the */
/*     algorithm uses QR factorizations with column pivoting for each */
/*     block J_k, k = 1:l, and applies the orthogonal matrix Q'_k to the */
/*     corresponding part of the last block column and of e. After all */
/*     block rows have been processed, the block rows are interchanged */
/*     so that the zeroed submatrices in the first l block columns are */
/*     moved to the bottom part. The same block row permutation Z is */
/*     also applied to the vector e. At the end of the first phase, */
/*     the structure of the processed matrix J is */

/*         /   R_1    0    ..   0   |  L^1_1  \ */
/*         |    0    R_2   ..   0   |  L^1_2  | */
/*         |    :     :    ..   :   |    :    | . */
/*         |    :     :    ..   :   |    :    | */
/*         |    0     0    ..  R_l  |  L^1_l  | */
/*         |    0     0    ..   0   |  L^2_1  | */
/*         |    :     :    ..   :   |    :    | */
/*         \    0     0    ..   0   |  L^2_l  / */

/*     In the second phase, the submatrix L^2_1:l is triangularized */
/*     using an additional QR factorization with pivoting. (The columns */
/*     of L^1_1:l are also permuted accordingly.) Therefore, the column */
/*     pivoting is restricted to each such local block column. */

/*     If l <= 1, the matrix J is triangularized in one phase, by one */
/*     QR factorization with pivoting. In this case, the column */
/*     pivoting is global. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     Feb. 22, 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Jacobian matrix, matrix algebra, */
/*     matrix operations, Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --dwork;
    --ipvt;
    --jnorms;
    --e;
    --j;
    --ipar;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*lipar < 4) {
	*info = -3;
    } else if (*fnorm < 0.) {
	*info = -4;
    } else if (*ldj < max(1,*n)) {
	*info = -6;
    } else {
	st = ipar[1];
	bn = ipar[2];
	bsm = ipar[3];
	bsn = ipar[4];
	nths = bn * bsn;
	mmn = bsm - bsn;
	if (bn > 0) {
	    m = bn * bsm;
	} else {
	    m = *n;
	}
/* Computing MIN */
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
	if (min(i__1,bsn) < 0) {
	    *info = -2;
	} else if (*n != nths + st) {
	    *info = -1;
	} else if (m < *n) {
	    *info = -2;
	} else if (*ldj < max(1,m)) {
	    *info = -6;
	} else {
	    if (*n == 0) {
		jwork = 1;
	    } else if (bn <= 1 || bsn == 0) {
		if (bn <= 1 && bsm == 1 && *n == 1) {
		    jwork = 1;
		} else {
		    jwork = (*n << 2) + 1;
		}
	    } else {
/* Computing MAX */
		i__1 = bsn * 3 + 1;
		jwork = bsn + max(i__1,st);
		if (bsm > bsn) {
/* Computing MAX */
		    i__1 = jwork, i__2 = (st << 2) + 1;
		    jwork = max(i__1,i__2);
		    if (bsm < bsn << 1) {
/* Computing MAX */
			i__1 = jwork, i__2 = mmn * (bn - 1);
			jwork = max(i__1,i__2);
		    }
		}
	    }
	    if (*ldwork < jwork) {
		*info = -12;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("NF01BS", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *gnorm = 0.;
    if (*n == 0) {
	*ldj = 1;
	dwork[1] = 1.;
	return 0;
    }

    if (bn <= 1 || bsn == 0) {

/*        Special case, l <= 1 or BSN = 0: the Jacobian is represented */
/*        as a full matrix. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

/*        Workspace: need:    4*N + 1; */
/*                   prefer:  3*N + ( N+1 )*NB. */

	md03bx_(&m, n, fnorm, &j[1], ldj, &e[1], &jnorms[1], gnorm, &ipvt[1], 
		&dwork[1], ldwork, info);
	return 0;
    }

/*     General case: l > 1 and BSN > 0. */
/*     Initialize the column pivoting indices. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipvt[i__] = 0;
/* L10: */
    }

/*     Compute the QR factorization with pivoting of J. */
/*     Pivoting is done separately on each block column of J. */

    wrkopt = 1;
    ibsn = 1;
    jl = *ldj * bsn + 1;
    jwork = bsn + 1;

    i__1 = m;
    i__2 = bsm;
    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += i__2) {

/*        Compute the QR factorization with pivoting of J_k, and apply Q' */
/*        to the corresponding part of the last block-column and of e. */
/*        Workspace: need:    4*BSN + 1; */
/*                   prefer:  3*BSN + ( BSN+1 )*NB. */

	i__3 = *ldwork - jwork + 1;
	dgeqp3_(&bsm, &bsn, &j[ibsm], ldj, &ipvt[ibsn], &dwork[1], &dwork[
		jwork], &i__3, info);
/* Computing MAX */
	i__3 = wrkopt, i__4 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__3,i__4);
	if (ibsm > 1) {

/*           Adjust the column pivoting indices. */

	    i__3 = ibsn + bsn - 1;
	    for (i__ = ibsn; i__ <= i__3; ++i__) {
		ipvt[i__] = ipvt[i__] + ibsn - 1;
/* L20: */
	    }

	}

	if (st > 0) {

/*           Workspace: need:    BSN + ST; */
/*                      prefer:  BSN + ST*NB. */

	    i__3 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &bsm, &st, &bsn, &j[ibsm], ldj, &
		    dwork[1], &j[jl], ldj, &dwork[jwork], &i__3, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
	    i__3 = wrkopt, i__4 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__3,i__4);
	}

/*        Workspace: need:    BSN + 1; */
/*                   prefer:  BSN + NB. */

	i__3 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &bsm, &c__1, &bsn, &j[ibsm], ldj, &dwork[
		1], &e[ibsm], &bsm, &dwork[jwork], &i__3, info, (ftnlen)4, (
		ftnlen)9);
	jl += bsm;
	ibsn += bsn;
/* L30: */
    }

    if (mmn > 0) {

/*        Case BSM > BSN. */
/*        Compute the original column norms for the first block column */
/*        of Jc. */
/*        Permute the rows of the first block column to move the zeroed */
/*        submatrices to the bottom. In the same loops, reshape the */
/*        first block column of R to have the leading dimension N. */

	l = ipvt[1];
	jnorms[l] = abs(j[1]);
	ibsm = bsm + 1;
	ibsn = bsn + 1;

	i__2 = bn - 1;
	for (k = 1; k <= i__2; ++k) {
	    j[ibsn] = j[ibsm];
	    l = ipvt[ibsn];
	    jnorms[l] = (d__1 = j[ibsn], abs(d__1));
	    ibsm += bsm;
	    ibsn += bsn;
/* L40: */
	}

	ibsn += st;

	i__2 = bsn;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    ibsm = (i__ - 1) * *ldj + 1;
	    jl = i__;

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {

		i__3 = i__ - 1;
		for (l = 0; l <= i__3; ++l) {
		    j[ibsn + l] = j[ibsm + l];
/* L45: */
		}

		l = ipvt[jl];
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
		ibsm += bsm;
		ibsn += bsn;
		jl += bsn;
/* L50: */
	    }

	    ibsn += st;
/* L60: */
	}

/*        Permute the rows of the second block column of Jc and of */
/*        the vector e. */

	jl = *ldj * bsn;
	if (bsm >= bsn << 1) {

/*           A swap operation can be used. */

	    i__2 = st;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ibsn = bsn + 1;

		i__1 = m;
		i__3 = bsm;
		for (ibsm = bsm + 1; i__3 < 0 ? ibsm >= i__1 : ibsm <= i__1; 
			ibsm += i__3) {
		    dswap_(&mmn, &j[jl + ibsm], &c__1, &j[jl + ibsn], &c__1);
		    ibsn += bsn;
/* L70: */
		}

		jl += *ldj;
/* L80: */
	    }

/*           Permute the rows of e. */

	    ibsn = bsn + 1;

	    i__2 = m;
	    i__3 = bsm;
	    for (ibsm = bsm + 1; i__3 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm 
		    += i__3) {
		dswap_(&mmn, &e[ibsm], &c__1, &e[ibsn], &c__1);
		ibsn += bsn;
/* L90: */
	    }

	} else {

/*           A swap operation cannot be used. */
/*           Workspace: need:    ( BSM-BSN )*( BN-1 ). */

	    i__3 = st;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ibsn = bsn + 1;
		jlm = jl + ibsn;
		jwork = 1;

		i__2 = m;
		i__1 = bsm;
		for (ibsm = bsm + 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; 
			ibsm += i__1) {
		    dcopy_(&mmn, &j[jlm], &c__1, &dwork[jwork], &c__1);

		    i__4 = jl + bsn - 1;
		    for (k = jl; k <= i__4; ++k) {
			j[ibsn + k] = j[ibsm + k];
/* L105: */
		    }

		    jlm += bsm;
		    ibsn += bsn;
		    jwork += mmn;
/* L100: */
		}

		i__1 = mmn * (bn - 1);
		dcopy_(&i__1, &dwork[1], &c__1, &j[jl + ibsn], &c__1);
		jl += *ldj;
/* L110: */
	    }

/*           Permute the rows of e. */

	    ibsn = bsn + 1;
	    jlm = ibsn;
	    jwork = 1;

	    i__3 = m;
	    i__1 = bsm;
	    for (ibsm = bsm + 1; i__1 < 0 ? ibsm >= i__3 : ibsm <= i__3; ibsm 
		    += i__1) {
		dcopy_(&mmn, &e[jlm], &c__1, &dwork[jwork], &c__1);

		i__2 = bsn - 1;
		for (k = 0; k <= i__2; ++k) {
		    e[ibsn + k] = e[ibsm + k];
/* L115: */
		}

		jlm += bsm;
		ibsn += bsn;
		jwork += mmn;
/* L120: */
	    }

	    i__1 = mmn * (bn - 1);
	    dcopy_(&i__1, &dwork[1], &c__1, &e[ibsn], &c__1);
	}

	if (st > 0) {

/*           Compute the QR factorization with pivoting of the submatrix */
/*           L^2_1:l, and apply Q' to the corresponding part of e. */

/*           Workspace: need:    4*ST + 1; */
/*                      prefer:  3*ST + ( ST+1 )*NB. */

	    jl = (*ldj + bn) * bsn + 1;
	    itau = 1;
	    jwork = itau + st;
	    i__1 = mmn * bn;
	    i__3 = *ldwork - jwork + 1;
	    dgeqp3_(&i__1, &st, &j[jl], ldj, &ipvt[nths + 1], &dwork[itau], &
		    dwork[jwork], &i__3, info);
/* Computing MAX */
	    i__1 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__1,i__3);

/*           Permute columns of the upper part of the second block */
/*           column of Jc. */

	    dlapmt_(&c_true, &nths, &st, &j[jl - nths], ldj, &ipvt[nths + 1]);

/*           Adjust the column pivoting indices. */

	    i__1 = *n;
	    for (i__ = nths + 1; i__ <= i__1; ++i__) {
		ipvt[i__] += nths;
/* L130: */
	    }

/*           Workspace: need:    ST + 1; */
/*                      prefer:  ST + NB. */

	    i__1 = mmn * bn;
	    i__3 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &i__1, &c__1, &st, &j[jl], ldj, &
		    dwork[itau], &e[ibsn], ldj, &dwork[jwork], &i__3, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__1,i__3);

/*           Reshape the second block column of R to have the leading */
/*           dimension N. */

	    ibsn = *n * bsn + 1;
	    dlacpy_("Full", n, &st, &j[*ldj * bsn + 1], ldj, &j[ibsn], n, (
		    ftnlen)4);

/*           Compute the original column norms for the second block */
/*           column. */

	    i__1 = *n;
	    for (i__ = nths + 1; i__ <= i__1; ++i__) {
		l = ipvt[i__];
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
		ibsn += *n;
/* L140: */
	    }

	}

    } else {

/*        Case BSM = BSN. */
/*        Compute the original column norms for the first block column */
/*        of Jc. */

	ibsn = 1;

	i__1 = bsn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jl = i__;

	    i__3 = bn;
	    for (k = 1; k <= i__3; ++k) {
		l = ipvt[jl];
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
		ibsn += bsn;
		jl += bsn;
/* L150: */
	    }

	    ibsn += st;
/* L160: */
	}

	i__1 = *n;
	for (i__ = nths + 1; i__ <= i__1; ++i__) {
	    ipvt[i__] = i__;
/* L170: */
	}

    }

/*     Compute the norm of the scaled gradient. */

    if (*fnorm != 0.) {

	i__1 = nths;
	i__3 = bsn;
	for (ibsn = 1; i__3 < 0 ? ibsn >= i__1 : ibsn <= i__1; ibsn += i__3) {
	    ibsni = ibsn;

	    i__2 = bsn;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		l = ipvt[ibsn + i__ - 1];
		if (jnorms[l] != 0.) {
		    sum = ddot_(&i__, &j[ibsni], &c__1, &e[ibsn], &c__1) / *
			    fnorm;
/* Computing MAX */
		    d__2 = *gnorm, d__3 = (d__1 = sum / jnorms[l], abs(d__1));
		    *gnorm = max(d__2,d__3);
		}
		ibsni += *n;
/* L180: */
	    }

/* L190: */
	}

	ibsni = *n * bsn + 1;

	i__3 = *n;
	for (i__ = nths + 1; i__ <= i__3; ++i__) {
	    l = ipvt[i__];
	    if (jnorms[l] != 0.) {
		sum = ddot_(&i__, &j[ibsni], &c__1, &e[1], &c__1) / *fnorm;
/* Computing MAX */
		d__2 = *gnorm, d__3 = (d__1 = sum / jnorms[l], abs(d__1));
		*gnorm = max(d__2,d__3);
	    }
	    ibsni += *n;
/* L200: */
	}

    }

    *ldj = *n;
    dwork[1] = (doublereal) wrkopt;
    return 0;

/* *** Last line of NF01BS *** */
} /* nf01bs_ */

