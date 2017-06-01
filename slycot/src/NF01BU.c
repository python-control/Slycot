/* NF01BU.f -- translated by f2c (version 20100827).
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
static doublereal c_b9 = 0.;
static doublereal c_b11 = 1.;
static integer c__0 = 0;

/* Subroutine */ int nf01bu_(char *stor, char *uplo, integer *n, integer *
	ipar, integer *lipar, doublereal *dpar, integer *ldpar, doublereal *j,
	 integer *ldj, doublereal *jtj, integer *ldjtj, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen stor_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal c__;
    static integer k, m, i1, bn, ii, jl, st, bsm, bsn;
    static doublereal tmp[1];
    static integer ibsm, ibsn, nbsn;
    static logical full;
    static integer itmp[1], nths;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     nf01bv_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dlaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To compute the matrix J'*J + c*I, for the Jacobian J as received */
/*     from SLICOT Library routine NF01BD: */

/*          /  dy(1)/dwb(1)  |  dy(1)/dtheta  \ */
/*     Jc = |       :        |       :        | . */
/*          \  dy(L)/dwb(L)  |  dy(L)/dtheta  / */

/*     This is a compressed representation of the actual structure */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_L  |  L_L  / */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STOR    CHARACTER*1 */
/*             Specifies the storage scheme for the symmetric */
/*             matrix J'*J + c*I, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the matrix J'*J + c*I is stored, */
/*             as follows: */
/*             = 'U' :  the upper triagular part is stored; */
/*             = 'L' :  the lower triagular part is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix J'*J + c*I. */
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
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             The leading NR-by-NC part of this array must contain */
/*             the (compressed) representation (Jc) of the Jacobian */
/*             matrix J, where NR = BSM if BN <= 1, and NR = BN*BSM, */
/*             if BN > 1. */

/*     LDJ     (input) INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NR). */

/*     JTJ     (output) DOUBLE PRECISION array, */
/*                      dimension (LDJTJ,N),    if STOR = 'F', */
/*                      dimension (N*(N+1)/2),  if STOR = 'P'. */
/*             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if */
/*             STOR = 'P') part of this array contains the upper or */
/*             lower triangle of the matrix J'*J + c*I, depending on */
/*             UPLO = 'U', or UPLO = 'L', respectively, stored either as */
/*             a two-dimensional, or one-dimensional array, depending */
/*             on STOR. */

/*     LDJTJ   INTEGER */
/*             The leading dimension of the array JTJ. */
/*             LDJTJ >= MAX(1,N), if STOR = 'F'. */
/*             LDJTJ >= 1,        if STOR = 'P'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             Currently, this array is not used. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix product is computed columnn-wise, exploiting the */
/*     symmetry. BLAS 3 routines DGEMM and DSYRK are used if STOR = 'F', */
/*     and BLAS 2 routine DGEMV is used if STOR = 'P'. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001. */

/*     REVISIONS */

/*     V. Sima, Dec. 2001, Mar. 2002. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations, */
/*     Wiener system. */

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

    /* Parameter adjustments */
    --ipar;
    --dpar;
    j_dim1 = *ldj;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --jtj;
    --dwork;

    /* Function Body */
    *info = 0;

    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

    if (! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lipar < 4) {
	*info = -5;
    } else if (*ldpar < 1) {
	*info = -7;
    } else if (*ldjtj < 1 || full && *ldjtj < *n) {
	*info = -11;
    } else if (*ldwork < 0) {
	*info = -13;
    } else {
	st = ipar[1];
	bn = ipar[2];
	bsm = ipar[3];
	bsn = ipar[4];
	nths = bn * bsn;
	if (bn > 1) {
	    m = bn * bsm;
	} else {
	    m = bsm;
	}
/* Computing MIN */
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
	if (min(i__1,bsn) < 0) {
	    *info = -4;
	} else if (*n != nths + st) {
	    *info = -3;
	} else if (*ldj < max(1,m)) {
	    *info = -9;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    c__ = dpar[1];

    if (bn <= 1 || bsn == 0 || bsm == 0) {

/*        Special case, l <= 1 or BSN = 0 or BSM = 0: the Jacobian is */
/*        represented as a full matrix. */

	itmp[0] = m;
	nf01bv_(stor, uplo, n, itmp, &c__1, &dpar[1], &c__1, &j[j_offset], 
		ldj, &jtj[1], ldjtj, &dwork[1], ldwork, info, (ftnlen)1, (
		ftnlen)1);
	return 0;
    }

/*     General case: l > 1, BSN > 0, BSM > 0. */

    jl = bsn + 1;

    if (full) {

	nbsn = *n * bsn;

	if (upper) {

/*           Compute the leading upper triangular part (full storage). */

	    dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[1], ldjtj, (ftnlen)1);
	    dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[j_offset], ldj, &
		    c_b11, &jtj[1], ldjtj, (ftnlen)1, (ftnlen)9);
	    ibsn = bsn;
	    i1 = nbsn + 1;

	    i__1 = m;
	    i__2 = bsm;
	    for (ibsm = bsm + 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm 
		    += i__2) {
		ii = i1 + ibsn;
		dlaset_("Full", &ibsn, &bsn, &c_b9, &c_b9, &jtj[i1], ldjtj, (
			ftnlen)4);
		i1 += nbsn;
		dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[ii], ldjtj, (
			ftnlen)1);
		dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[ibsm + 
			j_dim1], ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (
			ftnlen)9);
		ibsn += bsn;
/* L10: */
	    }

	    if (st > 0) {

/*              Compute the last block column. */

		i__2 = m;
		i__1 = bsm;
		for (ibsm = 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm +=
			 i__1) {
		    dgemm_("Transpose", "NoTranspose", &bsn, &st, &bsm, &
			    c_b11, &j[ibsm + j_dim1], ldj, &j[ibsm + jl * 
			    j_dim1], ldj, &c_b9, &jtj[i1], ldjtj, (ftnlen)9, (
			    ftnlen)11);
		    i1 += bsn;
/* L20: */
		}

		dlaset_(uplo, &st, &st, &c_b9, &c__, &jtj[i1], ldjtj, (ftnlen)
			1);
		dsyrk_(uplo, "Transpose", &st, &m, &c_b11, &j[jl * j_dim1 + 1]
			, ldj, &c_b11, &jtj[i1], ldjtj, (ftnlen)1, (ftnlen)9);
	    }

	} else {

/*           Compute the leading lower triangular part (full storage). */

	    ibsn = nths;
	    ii = 1;

	    i__1 = m;
	    i__2 = bsm;
	    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += 
		    i__2) {
		i1 = ii + bsn;
		dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[ii], ldjtj, (
			ftnlen)1);
		dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[ibsm + 
			j_dim1], ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (
			ftnlen)9);
		ibsn -= bsn;
		dlaset_("Full", &ibsn, &bsn, &c_b9, &c_b9, &jtj[i1], ldjtj, (
			ftnlen)4);
		ii = i1 + nbsn;
		if (st > 0) {
		    dgemm_("Transpose", "NoTranspose", &st, &bsn, &bsm, &
			    c_b11, &j[ibsm + jl * j_dim1], ldj, &j[ibsm + 
			    j_dim1], ldj, &c_b9, &jtj[i1 + ibsn], ldjtj, (
			    ftnlen)9, (ftnlen)11);
		}
/* L30: */
	    }

	    if (st > 0) {

/*              Compute the last diagonal block. */

		dlaset_(uplo, &st, &st, &c_b9, &c__, &jtj[ii], ldjtj, (ftnlen)
			1);
		dsyrk_(uplo, "Transpose", &st, &m, &c_b11, &j[jl * j_dim1 + 1]
			, ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (ftnlen)9);
	    }

	}

    } else {

	tmp[0] = 0.;

	if (upper) {

/*           Compute the leading upper triangular part (packed storage). */

	    ibsn = 0;
	    i1 = 1;

	    i__2 = m;
	    i__1 = bsm;
	    for (ibsm = 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm += 
		    i__1) {

		i__3 = bsn;
		for (k = 1; k <= i__3; ++k) {
		    ii = i1 + ibsn;
		    dcopy_(&ibsn, tmp, &c__0, &jtj[i1], &c__1);
		    dgemv_("Transpose", &bsm, &k, &c_b11, &j[ibsm + j_dim1], 
			    ldj, &j[ibsm + k * j_dim1], &c__1, &c_b9, &jtj[ii]
			    , &c__1, (ftnlen)9);
		    i1 = ii + k;
		    jtj[i1 - 1] += c__;
/* L40: */
		}

		ibsn += bsn;
/* L50: */
	    }

/*           Compute the last block column. */

	    i__1 = st;
	    for (k = 1; k <= i__1; ++k) {

		i__2 = m;
		i__3 = bsm;
		for (ibsm = 1; i__3 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm +=
			 i__3) {
		    dgemv_("Transpose", &bsm, &bsn, &c_b11, &j[ibsm + j_dim1],
			     ldj, &j[ibsm + (bsn + k) * j_dim1], &c__1, &c_b9,
			     &jtj[i1], &c__1, (ftnlen)9);
		    i1 += bsn;
/* L60: */
		}

		dgemv_("Transpose", &m, &k, &c_b11, &j[jl * j_dim1 + 1], ldj, 
			&j[(bsn + k) * j_dim1 + 1], &c__1, &c_b9, &jtj[i1], &
			c__1, (ftnlen)9);
		i1 += k;
		jtj[i1 - 1] += c__;
/* L70: */
	    }

	} else {

/*           Compute the leading lower triangular part (packed storage). */

	    ibsn = nths;
	    ii = 1;

	    i__1 = m;
	    i__3 = bsm;
	    for (ibsm = 1; i__3 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += 
		    i__3) {
		ibsn -= bsn;

		i__2 = bsn;
		for (k = 1; k <= i__2; ++k) {
		    i1 = ii + bsn - k + 1;
		    dcopy_(&ibsn, tmp, &c__0, &jtj[i1], &c__1);
		    i__4 = bsn - k + 1;
		    dgemv_("Transpose", &bsm, &i__4, &c_b11, &j[ibsm + k * 
			    j_dim1], ldj, &j[ibsm + k * j_dim1], &c__1, &c_b9,
			     &jtj[ii], &c__1, (ftnlen)9);
		    jtj[ii] += c__;
		    i1 += ibsn;
		    ii = i1 + st;
		    if (st > 0) {
			dgemv_("Transpose", &bsm, &st, &c_b11, &j[ibsm + jl * 
				j_dim1], ldj, &j[ibsm + k * j_dim1], &c__1, &
				c_b9, &jtj[i1], &c__1, (ftnlen)9);
		    }
/* L80: */
		}

/* L90: */
	    }

/*           Compute the last diagonal block. */

	    i__3 = st;
	    for (k = 1; k <= i__3; ++k) {
		i__1 = st - k + 1;
		dgemv_("Transpose", &m, &i__1, &c_b11, &j[(bsn + k) * j_dim1 
			+ 1], ldj, &j[(bsn + k) * j_dim1 + 1], &c__1, &c_b9, &
			jtj[ii], &c__1, (ftnlen)9);
		jtj[ii] += c__;
		ii = ii + st - k + 1;
/* L100: */
	    }

	}

    }

    return 0;

/* *** Last line of NF01BU *** */
} /* nf01bu_ */

