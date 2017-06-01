/* MB02CD.f -- translated by f2c (version 20100827).
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

static doublereal c_b17 = 1.;
static doublereal c_b19 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb02cd_(char *job, char *typet, integer *k, integer *n, 
	doublereal *t, integer *ldt, doublereal *g, integer *ldg, doublereal *
	r__, integer *ldr, doublereal *l, integer *ldl, doublereal *cs, 
	integer *lcs, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen job_len, ftnlen typet_len)
{
    /* System generated locals */
    integer g_dim1, g_offset, l_dim1, l_offset, r_dim1, r_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ierr;
    extern /* Subroutine */ int mb02cx_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), mb02cy_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical compg, compl, compr;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer starti, maxwrk, startr, startt;


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

/*     To compute the Cholesky factor and the generator and/or the */
/*     Cholesky factor of the inverse of a symmetric positive definite */
/*     (s.p.d.) block Toeplitz matrix T, defined by either its first */
/*     block row, or its first block column, depending on the routine */
/*     parameter TYPET. Transformation information is stored. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the output of the routine, as follows: */
/*             = 'G':  only computes the generator G of the inverse; */
/*             = 'R':  computes the generator G of the inverse and the */
/*                     Cholesky factor R of T, i.e., if TYPET = 'R', */
/*                     then R'*R = T, and if TYPET = 'C', then R*R' = T; */
/*             = 'L':  computes the generator G and the Cholesky factor L */
/*                     of the inverse, i.e., if TYPET = 'R', then */
/*                     L'*L = inv(T), and if TYPET = 'C', then */
/*                     L*L' = inv(T); */
/*             = 'A':  computes the generator G, the Cholesky factor L */
/*                     of the inverse and the Cholesky factor R of T; */
/*             = 'O':  only computes the Cholesky factor R of T. */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; if demanded, the Cholesky factors */
/*                     R and L are upper and lower triangular, */
/*                     respectively, and G contains the transposed */
/*                     generator of the inverse; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; if demanded, the Cholesky */
/*                     factors R and L are lower and upper triangular, */
/*                     respectively, and G contains the generator of the */
/*                     inverse. This choice results in a column oriented */
/*                     algorithm which is usually faster. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,N*K) / (LDT,K) */
/*             On entry, the leading K-by-N*K / N*K-by-K part of this */
/*             array must contain the first block row / column of an */
/*             s.p.d. block Toeplitz matrix. */
/*             On exit, if INFO = 0, then the leading K-by-N*K / N*K-by-K */
/*             part of this array contains, in the first K-by-K block, */
/*             the upper / lower Cholesky factor of T(1:K,1:K), and in */
/*             the remaining part, the Householder transformations */
/*             applied during the process. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),    if TYPET = 'R'; */
/*             LDT >= MAX(1,N*K),  if TYPET = 'C'. */

/*     G       (output)  DOUBLE PRECISION array, dimension */
/*             (LDG,N*K) / (LDG,2*K) */
/*             If INFO = 0 and JOB = 'G', 'R', 'L', or 'A', the leading */
/*             2*K-by-N*K / N*K-by-2*K part of this array contains, in */
/*             the first K-by-K block of the second block row / column, */
/*             the lower right block of L (necessary for updating */
/*             factorizations in SLICOT Library routine MB02DD), and */
/*             in the remaining part, the generator of the inverse of T. */
/*             Actually, to obtain a generator one has to set */
/*                 G(K+1:2*K, 1:K) = 0,    if TYPET = 'R'; */
/*                 G(1:K, K+1:2*K) = 0,    if TYPET = 'C'. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G. */
/*             LDG >= MAX(1,2*K),  if TYPET = 'R' and */
/*                                    JOB = 'G', 'R', 'L', or 'A'; */
/*             LDG >= MAX(1,N*K),  if TYPET = 'C' and */
/*                                    JOB = 'G', 'R', 'L', or 'A'; */
/*             LDG >= 1,           if JOB = 'O'. */

/*     R       (output)  DOUBLE PRECISION array, dimension (LDR,N*K) */
/*             If INFO = 0 and JOB = 'R', 'A', or 'O', then the leading */
/*             N*K-by-N*K part of this array contains the upper / lower */
/*             Cholesky factor of T. */
/*             The elements in the strictly lower / upper triangular part */
/*             are not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1,N*K),  if JOB = 'R', 'A', or 'O'; */
/*             LDR >= 1,           if JOB = 'G', or 'L'. */

/*     L       (output)  DOUBLE PRECISION array, dimension (LDL,N*K) */
/*             If INFO = 0 and JOB = 'L', or 'A', then the leading */
/*             N*K-by-N*K part of this array contains the lower / upper */
/*             Cholesky factor of the inverse of T. */
/*             The elements in the strictly upper / lower triangular part */
/*             are not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of the array L. */
/*             LDL >= MAX(1,N*K),  if JOB = 'L', or 'A'; */
/*             LDL >= 1,           if JOB = 'G', 'R', or 'O'. */

/*     CS      (output)  DOUBLE PRECISION array, dimension (LCS) */
/*             If INFO = 0, then the leading 3*(N-1)*K part of this */
/*             array contains information about the hyperbolic rotations */
/*             and Householder transformations applied during the */
/*             process. This information is needed for updating the */
/*             factorizations in SLICOT Library routine MB02DD. */

/*     LCS     INTEGER */
/*             The length of the array CS.  LCS >= 3*(N-1)*K. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,(N-1)*K). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The Toeplitz matrix */
/*                   associated with T is not (numerically) positive */
/*                   definite. */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */
/*                               3 2 */
/*     The algorithm requires 0(K N ) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2000, */
/*     February 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    l_dim1 = *ldl;
    l_offset = 1 + l_dim1;
    l -= l_offset;
    --cs;
    --dwork;

    /* Function Body */
    *info = 0;
    compl = lsame_(job, "L", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
    compg = lsame_(job, "G", (ftnlen)1, (ftnlen)1) || lsame_(job, "R", (
	    ftnlen)1, (ftnlen)1) || compl;
    compr = lsame_(job, "R", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1) || lsame_(job, "O", (ftnlen)1, (ftnlen)1);
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (compg || compr)) {
	*info = -1;
    } else if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*k < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < *n * *k) {
	*info = -6;
    } else if (*ldg < 1 || compg && (isrow && *ldg < *k << 1 || ! isrow && *
	    ldg < *n * *k)) {
	*info = -8;
    } else if (*ldr < 1 || compr && *ldr < *n * *k) {
	*info = -10;
    } else if (*ldl < 1 || compl && *ldl < *n * *k) {
	*info = -12;
    } else if (*lcs < (*n - 1) * 3 * *k) {
	*info = -14;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = (*n - 1) * *k;
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
	    i__1 = 1, i__2 = (*n - 1) * *k;
	    dwork[1] = (doublereal) max(i__1,i__2);
	    *info = -16;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02CD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*k,*n) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    maxwrk = 1;
    if (isrow) {

/*        T is the first block row of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

	dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

	    *info = 1;
	    return 0;
	}

	if (*n > 1) {
	    i__1 = (*n - 1) * *k;
	    dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &c_b17, 
		    &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], ldt, (
		    ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
	}

/*        Initialize the output matrices. */

	if (compg) {
	    i__1 = *k << 1;
	    i__2 = *n * *k;
	    dlaset_("All", &i__1, &i__2, &c_b19, &c_b19, &g[g_offset], ldg, (
		    ftnlen)3);
	    i__1 = *ldg + 1;
	    dlaset_("All", &c__1, k, &c_b17, &c_b17, &g[*k + 1 + g_dim1], &
		    i__1, (ftnlen)3);
	    dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, k, &c_b17, &t[
		    t_offset], ldt, &g[*k + 1 + g_dim1], ldg, (ftnlen)4, (
		    ftnlen)5, (ftnlen)9, (ftnlen)7);
	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &g[*k + 1 + (*k 
			+ 1) * g_dim1], ldg, (ftnlen)5);
	    }
	    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &g[g_offset], 
		    ldg, (ftnlen)5);
	}

	if (compl) {
	    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &l[l_offset], 
		    ldl, (ftnlen)5);
	}

	if (compr) {
	    i__1 = *n * *k;
	    dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);
	}

/*        Processing the generator. */

	if (compg) {

/*           Here we use G as working array for holding the generator. */
/*           T contains the second row of the generator. */
/*           G contains in its first block row the second row of the */
/*           inverse generator. */
/*           The second block row of G is partitioned as follows: */

/*           [ First block of the inverse generator, ... */
/*             First row of the generator, ... */
/*             The rest of the blocks of the inverse generator ] */

/*           The reason for the odd partitioning is that the first block */
/*           of the inverse generator will be thrown out at the end and */
/*           we want to avoid reordering. */

/*           (N-1)*K locations of DWORK are used by SLICOT Library */
/*           routine MB02CY. */

	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		startr = (i__ - 1) * *k + 1;
		starti = (*n - i__ + 1) * *k + 1;
		startt = (i__ - 2) * 3 * *k + 1;

/*              Transformations acting on the generator: */

		i__2 = *k * 3;
		mb02cx_("Row", k, k, k, &g[*k + 1 + (*k + 1) * g_dim1], ldg, &
			t[startr * t_dim1 + 1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)3);

		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);
		if (*n > i__) {
		    i__2 = (*n - i__) * *k;
		    i__3 = *k * 3;
		    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &g[*k + 1 + 
			    ((*k << 1) + 1) * g_dim1], ldg, &t[(startr + *k) *
			     t_dim1 + 1], ldt, &t[startr * t_dim1 + 1], ldt, &
			    cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
			    ftnlen)3, (ftnlen)11);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
		    maxwrk = max(i__2,i__3);
		}

		if (compr) {
		    i__2 = (*n - i__ + 1) * *k;
		    dlacpy_("Upper", k, &i__2, &g[*k + 1 + (*k + 1) * g_dim1],
			     ldg, &r__[startr + startr * r_dim1], ldr, (
			    ftnlen)5);
		}

/*              Transformations acting on the inverse generator: */

		dlaset_("All", k, k, &c_b19, &c_b19, &g[*k + 1 + starti * 
			g_dim1], ldg, (ftnlen)3);
		i__2 = *k * 3;
		mb02cy_("Row", "Triangular", k, k, k, k, &g[*k + 1 + g_dim1], 
			ldg, &g[startr * g_dim1 + 1], ldg, &t[startr * t_dim1 
			+ 1], ldt, &cs[startt], &i__2, &dwork[1], ldwork, &
			ierr, (ftnlen)3, (ftnlen)10);
/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);

		i__2 = (i__ - 1) * *k;
		i__3 = *k * 3;
		mb02cy_("Row", "NoStructure", k, k, &i__2, k, &g[*k + 1 + 
			starti * g_dim1], ldg, &g[g_offset], ldg, &t[startr * 
			t_dim1 + 1], ldt, &cs[startt], &i__3, &dwork[1], 
			ldwork, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);

		if (compl) {
		    i__2 = (i__ - 1) * *k;
		    dlacpy_("All", k, &i__2, &g[*k + 1 + starti * g_dim1], 
			    ldg, &l[startr + l_dim1], ldl, (ftnlen)3);
		    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &l[
			    startr + ((i__ - 1) * *k + 1) * l_dim1], ldl, (
			    ftnlen)5);
		}
/* L10: */
	    }

	} else {

/*           Here R is used as working array for holding the generator. */
/*           Again, T contains the second row of the generator. */
/*           The current row of R contains the first row of the */
/*           generator. */

	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[*k + 1 + (*
			k + 1) * r_dim1], ldr, (ftnlen)5);
	    }

	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		startr = (i__ - 1) * *k + 1;
		startt = (i__ - 2) * 3 * *k + 1;
		i__2 = *k * 3;
		mb02cx_("Row", k, k, k, &r__[startr + startr * r_dim1], ldr, &
			t[startr * t_dim1 + 1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)3);
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);
		if (*n > i__) {
		    i__2 = (*n - i__) * *k;
		    i__3 = *k * 3;
		    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &r__[startr 
			    + (startr + *k) * r_dim1], ldr, &t[(startr + *k) *
			     t_dim1 + 1], ldt, &t[startr * t_dim1 + 1], ldt, &
			    cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
			    ftnlen)3, (ftnlen)11);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
		    maxwrk = max(i__2,i__3);

		    i__2 = (*n - i__) * *k;
		    dlacpy_("Upper", k, &i__2, &r__[startr + startr * r_dim1],
			     ldr, &r__[startr + *k + (startr + *k) * r_dim1], 
			    ldr, (ftnlen)5);
		}
/* L20: */
	    }

	}

    } else {

/*        T is the first block column of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

	dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

	    *info = 1;
	    return 0;
	}

	if (*n > 1) {
	    i__1 = (*n - 1) * *k;
	    dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &c_b17,
		     &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)7);
	}

/*        Initialize the output matrices. */

	if (compg) {
	    i__1 = *n * *k;
	    i__2 = *k << 1;
	    dlaset_("All", &i__1, &i__2, &c_b19, &c_b19, &g[g_offset], ldg, (
		    ftnlen)3);
	    i__1 = *ldg + 1;
	    dlaset_("All", &c__1, k, &c_b17, &c_b17, &g[(*k + 1) * g_dim1 + 1]
		    , &i__1, (ftnlen)3);
	    dtrsm_("Right", "Lower", "Transpose", "NonUnit", k, k, &c_b17, &t[
		    t_offset], ldt, &g[(*k + 1) * g_dim1 + 1], ldg, (ftnlen)5,
		     (ftnlen)5, (ftnlen)9, (ftnlen)7);
	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &g[*k + 1 + (*k 
			+ 1) * g_dim1], ldg, (ftnlen)5);
	    }
	    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &g[
		    g_offset], ldg, (ftnlen)5);
	}

	if (compl) {
	    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &l[
		    l_offset], ldl, (ftnlen)5);
	}

	if (compr) {
	    i__1 = *n * *k;
	    dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);
	}

/*        Processing the generator. */

	if (compg) {

/*           Here we use G as working array for holding the generator. */
/*           T contains the second column of the generator. */
/*           G contains in its first block column the second column of */
/*           the inverse generator. */
/*           The second block column of G is partitioned as follows: */

/*           [ First block of the inverse generator; ... */
/*             First column of the generator; ... */
/*             The rest of the blocks of the inverse generator ] */

/*           The reason for the odd partitioning is that the first block */
/*           of the inverse generator will be thrown out at the end and */
/*           we want to avoid reordering. */

/*           (N-1)*K locations of DWORK are used by SLICOT Library */
/*           routine MB02CY. */

	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		startr = (i__ - 1) * *k + 1;
		starti = (*n - i__ + 1) * *k + 1;
		startt = (i__ - 2) * 3 * *k + 1;

/*              Transformations acting on the generator: */

		i__2 = *k * 3;
		mb02cx_("Column", k, k, k, &g[*k + 1 + (*k + 1) * g_dim1], 
			ldg, &t[startr + t_dim1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)6);

		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);
		if (*n > i__) {
		    i__2 = (*n - i__) * *k;
		    i__3 = *k * 3;
		    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &g[(*k <<
			     1) + 1 + (*k + 1) * g_dim1], ldg, &t[startr + *k 
			    + t_dim1], ldt, &t[startr + t_dim1], ldt, &cs[
			    startt], &i__3, &dwork[1], ldwork, &ierr, (ftnlen)
			    6, (ftnlen)11);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
		    maxwrk = max(i__2,i__3);
		}

		if (compr) {
		    i__2 = (*n - i__ + 1) * *k;
		    dlacpy_("Lower", &i__2, k, &g[*k + 1 + (*k + 1) * g_dim1],
			     ldg, &r__[startr + startr * r_dim1], ldr, (
			    ftnlen)5);
		}

/*              Transformations acting on the inverse generator: */

		dlaset_("All", k, k, &c_b19, &c_b19, &g[starti + (*k + 1) * 
			g_dim1], ldg, (ftnlen)3);
		i__2 = *k * 3;
		mb02cy_("Column", "Triangular", k, k, k, k, &g[(*k + 1) * 
			g_dim1 + 1], ldg, &g[startr + g_dim1], ldg, &t[startr 
			+ t_dim1], ldt, &cs[startt], &i__2, &dwork[1], ldwork,
			 &ierr, (ftnlen)6, (ftnlen)10);
/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);

		i__2 = (i__ - 1) * *k;
		i__3 = *k * 3;
		mb02cy_("Column", "NoStructure", k, k, &i__2, k, &g[starti + (
			*k + 1) * g_dim1], ldg, &g[g_offset], ldg, &t[startr 
			+ t_dim1], ldt, &cs[startt], &i__3, &dwork[1], ldwork,
			 &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);

		if (compl) {
		    i__2 = (i__ - 1) * *k;
		    dlacpy_("All", &i__2, k, &g[starti + (*k + 1) * g_dim1], 
			    ldg, &l[startr * l_dim1 + 1], ldl, (ftnlen)3);
		    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &l[
			    (i__ - 1) * *k + 1 + startr * l_dim1], ldl, (
			    ftnlen)5);
		}
/* L30: */
	    }

	} else {

/*           Here R is used as working array for holding the generator. */
/*           Again, T contains the second column of the generator. */
/*           The current column of R contains the first column of the */
/*           generator. */

	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[*k + 1 + (*
			k + 1) * r_dim1], ldr, (ftnlen)5);
	    }

	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		startr = (i__ - 1) * *k + 1;
		startt = (i__ - 2) * 3 * *k + 1;
		i__2 = *k * 3;
		mb02cx_("Column", k, k, k, &r__[startr + startr * r_dim1], 
			ldr, &t[startr + t_dim1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)6);
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

		    *info = 1;
		    return 0;
		}

/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[1];
		maxwrk = max(i__2,i__3);
		if (*n > i__) {
		    i__2 = (*n - i__) * *k;
		    i__3 = *k * 3;
		    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &r__[
			    startr + *k + startr * r_dim1], ldr, &t[startr + *
			    k + t_dim1], ldt, &t[startr + t_dim1], ldt, &cs[
			    startt], &i__3, &dwork[1], ldwork, &ierr, (ftnlen)
			    6, (ftnlen)11);
/* Computing MAX */
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
		    maxwrk = max(i__2,i__3);

		    i__2 = (*n - i__) * *k;
		    dlacpy_("Lower", &i__2, k, &r__[startr + startr * r_dim1],
			     ldr, &r__[startr + *k + (startr + *k) * r_dim1], 
			    ldr, (ftnlen)5);
		}
/* L40: */
	    }

	}
    }

    dwork[1] = (doublereal) maxwrk;

    return 0;

/* *** Last line of MB02CD *** */
} /* mb02cd_ */

