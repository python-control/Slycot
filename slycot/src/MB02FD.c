/* MB02FD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 1.;

/* Subroutine */ int mb02fd_(char *typet, integer *k, integer *n, integer *p, 
	integer *s, doublereal *t, integer *ldt, doublereal *r__, integer *
	ldr, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	typet_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, st, ierr;
    extern /* Subroutine */ int mb02cx_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), mb02cy_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpotrf_(char *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer maxwrk, countr, startr;


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

/*     To compute the incomplete Cholesky (ICC) factor of a symmetric */
/*     positive definite (s.p.d.) block Toeplitz matrix T, defined by */
/*     either its first block row, or its first block column, depending */
/*     on the routine parameter TYPET. */

/*     By subsequent calls of this routine, further rows / columns of */
/*     the Cholesky factor can be added. */
/*     Furthermore, the generator of the Schur complement of the leading */
/*     (P+S)*K-by-(P+S)*K block in T is available, which can be used, */
/*     e.g., for measuring the quality of the ICC factorization. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; the ICC factor R is upper */
/*                     trapezoidal; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; the ICC factor R is lower */
/*                     trapezoidal; this choice leads to better */
/*                     localized memory references and hence a faster */
/*                     algorithm. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     P       (input)  INTEGER */
/*             The number of previously computed block rows / columns */
/*             of R.  0 <= P <= N. */

/*     S       (input)  INTEGER */
/*             The number of block rows / columns of R to compute. */
/*             0 <= S <= N-P. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,(N-P)*K) / (LDT,K) */
/*             On entry, if P = 0, then the leading K-by-N*K / N*K-by-K */
/*             part of this array must contain the first block row / */
/*             column of an s.p.d. block Toeplitz matrix. */
/*             If P > 0, the leading K-by-(N-P)*K / (N-P)*K-by-K must */
/*             contain the negative generator of the Schur complement of */
/*             the leading P*K-by-P*K part in T, computed from previous */
/*             calls of this routine. */
/*             On exit, if INFO = 0, then the leading K-by-(N-P)*K / */
/*             (N-P)*K-by-K part of this array contains, in the first */
/*             K-by-K block, the upper / lower Cholesky factor of */
/*             T(1:K,1:K), in the following S-1 K-by-K blocks, the */
/*             Householder transformations applied during the process, */
/*             and in the remaining part, the negative generator of the */
/*             Schur complement of the leading (P+S)*K-by(P+S)*K part */
/*             in T. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),        if TYPET = 'R'; */
/*             LDT >= MAX(1,(N-P)*K),  if TYPET = 'C'. */

/*     R       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDR, N*K)       / (LDR, S*K )     if P = 0; */
/*             (LDR, (N-P+1)*K) / (LDR, (S+1)*K ) if P > 0. */
/*             On entry, if P > 0, then the leading K-by-(N-P+1)*K / */
/*             (N-P+1)*K-by-K part of this array must contain the */
/*             nonzero blocks of the last block row / column in the */
/*             ICC factor from a previous call of this routine. Note that */
/*             this part is identical with the positive generator of */
/*             the Schur complement of the leading P*K-by-P*K part in T. */
/*             If P = 0, then R is only an output parameter. */
/*             On exit, if INFO = 0 and P = 0, then the leading */
/*             S*K-by-N*K / N*K-by-S*K part of this array contains the */
/*             upper / lower trapezoidal ICC factor. */
/*             On exit, if INFO = 0 and P > 0, then the leading */
/*             (S+1)*K-by-(N-P+1)*K / (N-P+1)*K-by-(S+1)*K part of this */
/*             array contains the upper / lower trapezoidal part of the */
/*             P-th to (P+S)-th block rows / columns of the ICC factor. */
/*             The elements in the strictly lower / upper trapezoidal */
/*             part are not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1, S*K ),        if TYPET = 'R' and P = 0; */
/*             LDR >= MAX(1, (S+1)*K ),    if TYPET = 'R' and P > 0; */
/*             LDR >= MAX(1, N*K ),        if TYPET = 'C' and P = 0; */
/*             LDR >= MAX(1, (N-P+1)*K ),  if TYPET = 'C' and P > 0. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -11,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,(N+1)*K,4*K),   if P = 0; */
/*             LDWORK >= MAX(1,(N-P+2)*K,4*K), if P > 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed; the Toeplitz matrix */
/*                   associated with T is not (numerically) positive */
/*                   definite in its leading (P+S)*K-by-(P+S)*K part. */

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
/*                               3 */
/*     The algorithm requires 0(K S (N-P)) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001, */
/*     Mar. 2004. */

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
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*k < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*p < 0 || *p > *n) {
	*info = -4;
    } else if (*s < 0 || *s > *n - *p) {
	*info = -5;
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < (*n - *p) *
	     *k) {
	*info = -7;
    } else if (*ldr < 1 || isrow && *p == 0 && *ldr < *s * *k || isrow && *p 
	    > 0 && *ldr < (*s + 1) * *k || ! isrow && *p == 0 && *ldr < *n * *
	    k || ! isrow && *p > 0 && *ldr < (*n - *p + 1) * *k) {
	*info = -9;
    } else {
	if (*p == 0) {
	    countr = (*n + 1) * *k;
	} else {
	    countr = (*n - *p + 2) * *k;
	}
/* Computing MAX */
	i__1 = countr, i__2 = *k << 2;
	countr = max(i__1,i__2);
	if (*ldwork < max(1,countr)) {
	    dwork[1] = (doublereal) max(1,countr);
	    *info = -11;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02FD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*k,*n);
    if (min(i__1,*s) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    maxwrk = 1;

    if (isrow) {

	if (*p == 0) {

/*           T is the first block row of a block Toeplitz matrix. */
/*           Bring T to proper form by triangularizing its first block. */

	    dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }

	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &
			c_b10, &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], 
			ldt, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
	    }
	    i__1 = *n * *k;
	    dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);

	    if (*s == 1) {
		dwork[1] = 1.;
		return 0;
	    }

	    st = 2;
	    countr = (*n - 1) * *k;
	} else {
	    st = 1;
	    countr = (*n - *p) * *k;
	}

	startr = 1;

	i__1 = *s;
	for (i__ = st; i__ <= i__1; ++i__) {
	    dlacpy_("Upper", k, &countr, &r__[startr + startr * r_dim1], ldr, 
		    &r__[startr + *k + (startr + *k) * r_dim1], ldr, (ftnlen)
		    5);
	    startr += *k;
	    countr -= *k;
	    i__2 = *k * 3;
	    i__3 = *ldwork - *k * 3;
	    mb02cx_("Row", k, k, k, &r__[startr + startr * r_dim1], ldr, &t[
		    startr * t_dim1 + 1], ldt, &dwork[1], &i__2, &dwork[*k * 
		    3 + 1], &i__3, &ierr, (ftnlen)3);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }

/* Computing MAX */
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
	    maxwrk = max(i__2,i__3);
	    i__2 = *k * 3;
	    i__3 = *ldwork - *k * 3;
	    mb02cy_("Row", "NoStructure", k, k, &countr, k, &r__[startr + (
		    startr + *k) * r_dim1], ldr, &t[(startr + *k) * t_dim1 + 
		    1], ldt, &t[startr * t_dim1 + 1], ldt, &dwork[1], &i__2, &
		    dwork[*k * 3 + 1], &i__3, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
	    maxwrk = max(i__2,i__3);
/* L10: */
	}

    } else {

	if (*p == 0) {

/*           T is the first block column of a block Toeplitz matrix. */
/*           Bring T to proper form by triangularizing its first block. */

	    dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }

	    if (*n > 1) {
		i__1 = (*n - 1) * *k;
		dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &
			c_b10, &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (
			ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)7);
	    }
	    i__1 = *n * *k;
	    dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);

	    if (*s == 1) {
		dwork[1] = 1.;
		return 0;
	    }

	    st = 2;
	    countr = (*n - 1) * *k;
	} else {
	    st = 1;
	    countr = (*n - *p) * *k;
	}

	startr = 1;

	i__1 = *s;
	for (i__ = st; i__ <= i__1; ++i__) {
	    dlacpy_("Lower", &countr, k, &r__[startr + startr * r_dim1], ldr, 
		    &r__[startr + *k + (startr + *k) * r_dim1], ldr, (ftnlen)
		    5);
	    startr += *k;
	    countr -= *k;
	    i__2 = *k * 3;
	    i__3 = *ldwork - *k * 3;
	    mb02cx_("Column", k, k, k, &r__[startr + startr * r_dim1], ldr, &
		    t[startr + t_dim1], ldt, &dwork[1], &i__2, &dwork[*k * 3 
		    + 1], &i__3, &ierr, (ftnlen)6);
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

		*info = 1;
		return 0;
	    }

/* Computing MAX */
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
	    maxwrk = max(i__2,i__3);
	    i__2 = *k * 3;
	    i__3 = *ldwork - *k * 3;
	    mb02cy_("Column", "NoStructure", k, k, &countr, k, &r__[startr + *
		    k + startr * r_dim1], ldr, &t[startr + *k + t_dim1], ldt, 
		    &t[startr + t_dim1], ldt, &dwork[1], &i__2, &dwork[*k * 3 
		    + 1], &i__3, &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
	    maxwrk = max(i__2,i__3);
/* L20: */
	}

    }

    dwork[1] = (doublereal) maxwrk;

    return 0;

/* *** Last line of MB02FD *** */
} /* mb02fd_ */

