/* MB03RD.f -- translated by f2c (version 20100827).
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

static doublereal c_b20 = 0.;
static doublereal c_b28 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb03rd_(char *jobx, char *sort, integer *n, doublereal *
	pmax, doublereal *a, integer *lda, doublereal *x, integer *ldx, 
	integer *nblcks, integer *blsize, doublereal *wr, doublereal *wi, 
	doublereal *tol, doublereal *dwork, integer *info, ftnlen jobx_len, 
	ftnlen sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, d__;
    static integer i__, j, k, l, l11, l22;
    static doublereal sc;
    static integer da11, da22;
    static doublereal cav, rav;
    static integer l22m1;
    static doublereal edif, emax;
    static char jobv[1];
    static integer ierr;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03qx_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), mb03rx_(char *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, ftnlen), 
	    mb03ry_(integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, integer *);
    static logical ljobx, lsorn, lsors, lsort;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safemn;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal thresh;


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

/*     To reduce a matrix A in real Schur form to a block-diagonal form */
/*     using well-conditioned non-orthogonal similarity transformations. */
/*     The condition numbers of the transformations used for reduction */
/*     are roughly bounded by PMAX*PMAX, where PMAX is a given value. */
/*     The transformations are optionally postmultiplied in a given */
/*     matrix X. The real Schur form is optionally ordered, so that */
/*     clustered eigenvalues are grouped in the same block. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBX    CHARACTER*1 */
/*             Specifies whether or not the transformations are */
/*             accumulated, as follows: */
/*             = 'N':  The transformations are not accumulated; */
/*             = 'U':  The transformations are accumulated in X (the */
/*                     given matrix X is updated). */

/*     SORT    CHARACTER*1 */
/*             Specifies whether or not the diagonal blocks of the real */
/*             Schur form are reordered, as follows: */
/*             = 'N':  The diagonal blocks are not reordered; */
/*             = 'S':  The diagonal blocks are reordered before each */
/*                     step of reduction, so that clustered eigenvalues */
/*                     appear in the same block; */
/*             = 'C':  The diagonal blocks are not reordered, but the */
/*                     "closest-neighbour" strategy is used instead of */
/*                     the standard "closest to the mean" strategy */
/*                     (see METHOD); */
/*             = 'B':  The diagonal blocks are reordered before each */
/*                     step of reduction, and the "closest-neighbour" */
/*                     strategy is used (see METHOD). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and X.  N >= 0. */

/*     PMAX    (input) DOUBLE PRECISION */
/*             An upper bound for the infinity norm of elementary */
/*             submatrices of the individual transformations used for */
/*             reduction (see METHOD).  PMAX >= 1.0D0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A to be block-diagonalized, in real */
/*             Schur form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the computed block-diagonal matrix, in real Schur */
/*             canonical form. The non-diagonal blocks are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if JOBX = 'U', the leading N-by-N part of this */
/*             array must contain a given matrix X. */
/*             On exit, if JOBX = 'U', the leading N-by-N part of this */
/*             array contains the product of the given matrix X and the */
/*             transformation matrix that reduced A to block-diagonal */
/*             form. The transformation matrix is itself a product of */
/*             non-orthogonal similarity transformations having elements */
/*             with magnitude less than or equal to PMAX. */
/*             If JOBX = 'N', this array is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X. */
/*             LDX >= 1,        if JOBX = 'N'; */
/*             LDX >= MAX(1,N), if JOBX = 'U'. */

/*     NBLCKS  (output) INTEGER */
/*             The number of diagonal blocks of the matrix A. */

/*     BLSIZE  (output) INTEGER array, dimension (N) */
/*             The first NBLCKS elements of this array contain the orders */
/*             of the resulting diagonal blocks of the matrix A. */

/*     WR,     (output) DOUBLE PRECISION arrays, dimension (N) */
/*     WI      These arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the matrix A. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in the ordering of the diagonal */
/*             blocks of the real Schur form matrix. */
/*             If the user sets TOL > 0, then the given value of TOL is */
/*             used as an absolute tolerance: a block i and a temporarily */
/*             fixed block 1 (the first block of the current trailing */
/*             submatrix to be reduced) are considered to belong to the */
/*             same cluster if their eigenvalues satisfy */

/*               | lambda_1 - lambda_i | <= TOL. */

/*             If the user sets TOL < 0, then the given value of TOL is */
/*             used as a relative tolerance: a block i and a temporarily */
/*             fixed block 1 are considered to belong to the same cluster */
/*             if their eigenvalues satisfy, for j = 1, ..., N, */

/*               | lambda_1 - lambda_i | <= | TOL | * max | lambda_j |. */

/*             If the user sets TOL = 0, then an implicitly computed, */
/*             default tolerance, defined by TOL = SQRT( SQRT( EPS ) ) */
/*             is used instead, as a relative tolerance, where EPS is */
/*             the machine precision (see LAPACK Library routine DLAMCH). */
/*             If SORT = 'N' or 'C', this parameter is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Consider first that SORT = 'N'. Let */

/*            ( A    A   ) */
/*            (  11   12 ) */
/*        A = (          ), */
/*            ( 0    A   ) */
/*            (       22 ) */

/*     be the given matrix in real Schur form, where initially A   is the */
/*                                                              11 */
/*     first diagonal block of dimension 1-by-1 or 2-by-2. An attempt is */
/*     made to compute a transformation matrix X of the form */

/*            ( I   P ) */
/*        X = (       )                                               (1) */
/*            ( 0   I ) */

/*     (partitioned as A), so that */

/*                 ( A     0  ) */
/*         -1      (  11      ) */
/*        X  A X = (          ), */
/*                 ( 0    A   ) */
/*                 (       22 ) */

/*     and the elements of P do not exceed the value PMAX in magnitude. */
/*     An adaptation of the standard method for solving Sylvester */
/*     equations [1], which controls the magnitude of the individual */
/*     elements of the computed solution [2], is used to obtain matrix P. */
/*     When this attempt failed, an 1-by-1 (or 2-by-2) diagonal block of */
/*     A  , whose eigenvalue(s) is (are) the closest to the mean of those */
/*      22 */
/*     of A   is selected, and moved by orthogonal similarity */
/*         11 */
/*     transformations in the leading position of A  ; the moved diagonal */
/*                                                 22 */
/*     block is then added to the block A  , increasing its order by 1 */
/*                                       11 */
/*     (or 2). Another attempt is made to compute a suitable */
/*     transformation matrix X with the new definitions of the blocks A */
/*                                                                     11 */
/*     and A  . After a successful transformation matrix X has been */
/*          22 */
/*     obtained, it postmultiplies the current transformation matrix */
/*     (if JOBX = 'U'), and the whole procedure is repeated for the */
/*     matrix A  . */
/*             22 */

/*     When SORT = 'S', the diagonal blocks of the real Schur form are */
/*     reordered before each step of the reduction, so that each cluster */
/*     of eigenvalues, defined as specified in the definition of TOL, */
/*     appears in adjacent blocks. The blocks for each cluster are merged */
/*     together, and the procedure described above is applied to the */
/*     larger blocks. Using the option SORT = 'S' will usually provide */
/*     better efficiency than the standard option (SORT = 'N'), proposed */
/*     in [2], because there could be no or few unsuccessful attempts */
/*     to compute individual transformation matrices X of the form (1). */
/*     However, the resulting dimensions of the blocks are usually */
/*     larger; this could make subsequent calculations less efficient. */

/*     When SORT = 'C' or 'B', the procedure is similar to that for */
/*     SORT = 'N' or 'S', respectively, but the block of A   whose */
/*                                                        22 */
/*     eigenvalue(s) is (are) the closest to those of A   (not to their */
/*                                                     11 */
/*     mean) is selected and moved to the leading position of A  . This */
/*                                                             22 */
/*     is called the "closest-neighbour" strategy. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Bavely, C. and Stewart, G.W. */
/*         An Algorithm for Computing Reducing Subspaces by Block */
/*         Diagonalization. */
/*         SIAM J. Numer. Anal., 16, pp. 359-367, 1979. */

/*     [3] Demmel, J. */
/*         The Condition Number of Equivalence Transformations that */
/*         Block Diagonalize Matrix Pencils. */
/*         SIAM J. Numer. Anal., 20, pp. 599-610, 1983. */

/*     NUMERICAL ASPECTS */
/*                                       3                     4 */
/*     The algorithm usually requires 0(N ) operations, but 0(N ) are */
/*     possible in the worst case, when all diagonal blocks in the real */
/*     Schur form of A are 1-by-1, and the matrix cannot be diagonalized */
/*     by well-conditioned transformations. */

/*     FURTHER COMMENTS */

/*     The individual non-orthogonal transformation matrices used in the */
/*     reduction of A to a block-diagonal form have condition numbers */
/*     of the order PMAX*PMAX. This does not guarantee that their product */
/*     is well-conditioned enough. The routine can be easily modified to */
/*     provide estimates for the condition numbers of the clusters of */
/*     eigenvalues. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Partly based on the RASP routine BDIAG by A. Varga, German */
/*     Aerospace Center, DLR Oberpfaffenhofen. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003. */

/*     KEYWORDS */

/*     Diagonalization, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --blsize;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobx = lsame_(jobx, "U", (ftnlen)1, (ftnlen)1);
    lsorn = lsame_(sort, "N", (ftnlen)1, (ftnlen)1);
    lsors = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);
    lsort = lsame_(sort, "B", (ftnlen)1, (ftnlen)1) || lsors;
    if (! ljobx && ! lsame_(jobx, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! lsorn && ! lsort && ! lsame_(sort, "C", (ftnlen)1, (ftnlen)1)
	    ) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*pmax < 1.) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldx < 1 || ljobx && *ldx < *n) {
	*info = -8;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB03RD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *nblcks = 0;
    if (*n == 0) {
	return 0;
    }

/*     Set the "safe" minimum positive number with representable */
/*     reciprocal, and set JOBV parameter for MB03RX routine. */

    safemn = dlamch_("Safe minimum", (ftnlen)12);
    sc = 1. / safemn;
    dlabad_(&safemn, &sc);
    safemn /= dlamch_("Precision", (ftnlen)9);
    *(unsigned char *)jobv = *(unsigned char *)jobx;
    if (ljobx) {
	*(unsigned char *)jobv = 'V';
    }

/*     Compute the eigenvalues of A and set the tolerance for reordering */
/*     the eigenvalues in clusters, if needed. */

    mb03qx_(n, &a[a_offset], lda, &wr[1], &wi[1], info);

    if (lsort) {
	thresh = abs(*tol);
	if (thresh == 0.) {

/*           Use the default tolerance in ordering the blocks. */

	    thresh = sqrt(sqrt(dlamch_("Epsilon", (ftnlen)7)));
	}

	if (*tol <= 0.) {

/*           Use a relative tolerance. Find max | lambda_j |, j = 1 : N. */

	    emax = 0.;
	    l = 1;
/*           WHILE ( L.LE.N ) DO */
L10:
	    if (l <= *n) {
		if (wi[l] == 0.) {
/* Computing MAX */
		    d__2 = emax, d__3 = (d__1 = wr[l], abs(d__1));
		    emax = max(d__2,d__3);
		    ++l;
		} else {
/* Computing MAX */
		    d__1 = emax, d__2 = dlapy2_(&wr[l], &wi[l]);
		    emax = max(d__1,d__2);
		    l += 2;
		}
		goto L10;
	    }
/*           END WHILE 10 */
	    thresh *= emax;
	}
    }

/*     Define the following submatrices of A: */
/*     A11, the DA11-by-DA11 block in position (L11,L11); */
/*     A22, the DA22-by-DA22 block in position (L22,L22); */
/*     A12, the DA11-by-DA22 block in position (L11,L22); */
/*     A21, the DA22-by-DA11 block in position (L22,L11) (null initially */
/*                                                        and finally). */
/*     The following loop uses L11 as loop variable and try to separate a */
/*     block in position (L11,L11), with possibly clustered eigenvalues, */
/*     separated by the other eigenvalues (in the block A22). */

    l11 = 1;
/*     WHILE ( L11.LE.N ) DO */
L20:
    if (l11 <= *n) {
	++(*nblcks);
	if (wi[l11] == 0.) {
	    da11 = 1;
	} else {
	    da11 = 2;
	}

	if (lsort) {

/*           The following loop, using K as loop variable, finds the */
/*           blocks whose eigenvalues are close to those of A11 and */
/*           moves these blocks (if any) to the leading position of A22. */

	    l22 = l11 + da11;
	    k = l22;
/*           WHILE ( K.LE.N ) DO */
L30:
	    if (k <= *n) {
		d__1 = wr[l11] - wr[k];
		d__2 = wi[l11] - wi[k];
		edif = dlapy2_(&d__1, &d__2);
		if (edif <= thresh) {

/*                 An 1x1 or a 2x2 block of A22 has been found so that */

/*                    abs( lambda_1 - lambda_k ) <= THRESH */

/*                 where lambda_1 and lambda_k denote an eigenvalue */
/*                 of A11 and of that block in A22, respectively. */
/*                 Try to move that block to the leading position of A22. */

		    mb03rx_(jobv, n, &l22, &k, &a[a_offset], lda, &x[x_offset]
			    , ldx, &wr[1], &wi[1], &dwork[1], (ftnlen)1);

/*                 Extend A11 with the leading block of A22. */

		    if (wi[l22] == 0.) {
			++da11;
		    } else {
			da11 += 2;
		    }
		    l22 = l11 + da11;
		}
		if (wi[k] == 0.) {
		    ++k;
		} else {
		    k += 2;
		}
		goto L30;
	    }
/*           END WHILE 30 */
	}

/*        The following loop uses L22 as loop variable and forms a */
/*        separable DA11-by-DA11 block A11 in position (L11,L11). */

	l22 = l11 + da11;
	l22m1 = l22 - 1;
/*        WHILE ( L22.LE.N ) DO */
L40:
	if (l22 <= *n) {
	    da22 = *n - l22m1;

/*           Try to separate the block A11 of order DA11 by using a */
/*           well-conditioned similarity transformation. */

/*           First save A12' in the block A21. */

	    ma02ad_("Full", &da11, &da22, &a[l11 + l22 * a_dim1], lda, &a[l22 
		    + l11 * a_dim1], lda, (ftnlen)4);

/*           Solve  -A11*P + P*A22 = A12. */

	    mb03ry_(&da11, &da22, pmax, &a[l11 + l11 * a_dim1], lda, &a[l22 + 
		    l22 * a_dim1], lda, &a[l11 + l22 * a_dim1], lda, &ierr);

	    if (ierr == 1) {

/*              The annihilation of A12 failed. Restore A12 and A21. */

		ma02ad_("Full", &da22, &da11, &a[l22 + l11 * a_dim1], lda, &a[
			l11 + l22 * a_dim1], lda, (ftnlen)4);
		dlaset_("Full", &da22, &da11, &c_b20, &c_b20, &a[l22 + l11 * 
			a_dim1], lda, (ftnlen)4);

		if (lsorn || lsors) {

/*                 Extend A11 with an 1x1 or 2x2 block of A22 having the */
/*                 nearest eigenvalues to the mean of eigenvalues of A11 */
/*                 and resume the loop. */
/*                 First compute the mean of eigenvalues of A11. */

		    rav = 0.;
		    cav = 0.;

		    i__1 = l22m1;
		    for (i__ = l11; i__ <= i__1; ++i__) {
			rav += wr[i__];
			cav += (d__1 = wi[i__], abs(d__1));
/* L50: */
		    }

		    rav /= da11;
		    cav /= da11;

/*                 Loop to find the eigenvalue of A22 nearest to the */
/*                 above computed mean. */

		    d__1 = rav - wr[l22];
		    d__2 = cav - wi[l22];
		    d__ = dlapy2_(&d__1, &d__2);
		    k = l22;
		    if (wi[l22] == 0.) {
			l = l22 + 1;
		    } else {
			l = l22 + 2;
		    }
/*                 WHILE ( L.LE.N ) DO */
L60:
		    if (l <= *n) {
			d__1 = rav - wr[l];
			d__2 = cav - wi[l];
			c__ = dlapy2_(&d__1, &d__2);
			if (c__ < d__) {
			    d__ = c__;
			    k = l;
			}
			if (wi[l] == 0.) {
			    ++l;
			} else {
			    l += 2;
			}
			goto L60;
		    }
/*                 END WHILE 60 */

		} else {

/*                 Extend A11 with an 1x1 or 2x2 block of A22 having the */
/*                 nearest eigenvalues to the cluster of eigenvalues of */
/*                 A11 and resume the loop. */

/*                 Loop to find the eigenvalue of A22 of minimum distance */
/*                 to the cluster. */

		    d__ = sc;
		    l = l22;
		    k = l22;
/*                 WHILE ( L.LE.N ) DO */
L70:
		    if (l <= *n) {
			i__ = l11;
/*                    WHILE ( I.LE.L22M1 ) DO */
L80:
			if (i__ <= l22m1) {
			    d__1 = wr[i__] - wr[l];
			    d__2 = wi[i__] - wi[l];
			    c__ = dlapy2_(&d__1, &d__2);
			    if (c__ < d__) {
				d__ = c__;
				k = l;
			    }
			    if (wi[i__] == 0.) {
				++i__;
			    } else {
				i__ += 2;
			    }
			    goto L80;
			}
/*                    END WHILE 80 */
			if (wi[l] == 0.) {
			    ++l;
			} else {
			    l += 2;
			}
			goto L70;
		    }
/*                 END WHILE 70 */
		}

/*              Try to move block found to the leading position of A22. */

		mb03rx_(jobv, n, &l22, &k, &a[a_offset], lda, &x[x_offset], 
			ldx, &wr[1], &wi[1], &dwork[1], (ftnlen)1);

/*              Extend A11 with the leading block of A22. */

		if (wi[l22] == 0.) {
		    ++da11;
		} else {
		    da11 += 2;
		}
		l22 = l11 + da11;
		l22m1 = l22 - 1;
		goto L40;
	    }
	}
/*        END WHILE 40 */

	if (ljobx) {

/*           Accumulate the transformation in X. */
/*           Only columns L22, ..., N are modified. */

	    if (l22 <= *n) {
		dgemm_("No transpose", "No transpose", n, &da22, &da11, &
			c_b28, &x[l11 * x_dim1 + 1], ldx, &a[l11 + l22 * 
			a_dim1], lda, &c_b28, &x[l22 * x_dim1 + 1], ldx, (
			ftnlen)12, (ftnlen)12);
	    }

/*           Scale to unity the (non-zero) columns of X which will be */
/*           no more modified and transform A11 accordingly. */

	    i__1 = l22m1;
	    for (j = l11; j <= i__1; ++j) {
		sc = dnrm2_(n, &x[j * x_dim1 + 1], &c__1);
		if (sc > safemn) {
		    dscal_(&da11, &sc, &a[j + l11 * a_dim1], lda);
		    sc = 1. / sc;
		    dscal_(n, &sc, &x[j * x_dim1 + 1], &c__1);
		    dscal_(&da11, &sc, &a[l11 + j * a_dim1], &c__1);
		}
/* L90: */
	    }

	}
	if (l22 <= *n) {

/*           Set A12 and A21 to zero. */

	    dlaset_("Full", &da11, &da22, &c_b20, &c_b20, &a[l11 + l22 * 
		    a_dim1], lda, (ftnlen)4);
	    dlaset_("Full", &da22, &da11, &c_b20, &c_b20, &a[l22 + l11 * 
		    a_dim1], lda, (ftnlen)4);
	}

/*        Store the orders of the diagonal blocks in BLSIZE. */

	blsize[*nblcks] = da11;
	l11 = l22;
	goto L20;
    }
/*     END WHILE 20 */

    return 0;
/* *** Last line of MB03RD *** */
} /* mb03rd_ */

