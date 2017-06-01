/* MB03SD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 1.;
static doublereal c_b10 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb03sd_(char *jobscl, integer *n, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *wr, doublereal *wi, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen jobscl_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m;
    static doublereal x, y;
    static integer n2, bl, jw, ihi, ilo;
    static doublereal swap;
    static logical blas3;
    extern /* Subroutine */ int ma01ad_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), ma02ed_(char *, integer *, 
	    doublereal *, integer *, ftnlen);
    static logical scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical block;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dsymm_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal dummy[1];
    static integer jwork;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dgebal_(char *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, ftnlen), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer ignore;
    extern /* Subroutine */ int dhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical sorted;


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

/*     To compute the eigenvalues of an N-by-N square-reduced Hamiltonian */
/*     matrix */

/*                ( A'   G'  ) */
/*         H'  =  (        T ).                                       (1) */
/*                ( Q'  -A'  ) */

/*     Here, A' is an N-by-N matrix, and G' and Q' are symmetric N-by-N */
/*     matrices.  It is assumed without a check that H' is square- */
/*     reduced, i.e., that */

/*           2    ( A''   G'' ) */
/*         H'  =  (         T )    with A'' upper Hessenberg.         (2) */
/*                ( 0    A''  ) */

/*                            T                2 */
/*     (Equivalently, Q'A'- A' Q' = 0, A'' = A' + G'Q', and for i > j+1, */
/*      A''(i,j) = 0.)  Ordinarily, H' is the output from SLICOT Library */
/*     routine MB04ZD. The eigenvalues of H' are computed as the square */
/*     roots of the eigenvalues of A''. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBSCL  CHARACTER*1 */
/*             Specifies whether or not balancing operations should */
/*             be performed by the LAPACK subroutine DGEBAL on the */
/*             Hessenberg matrix A'' in (2), as follows: */
/*             = 'N':  do not use balancing; */
/*             = 'S':  do scaling in order to equilibrate the rows */
/*                     and columns of A''. */
/*             See LAPACK subroutine DGEBAL and Section METHOD below. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             upper left block A' of the square-reduced Hamiltonian */
/*             matrix H' in (1), as produced by SLICOT Library routine */
/*             MB04ZD. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1) */
/*             The leading N-by-N lower triangular part of this array */
/*             must contain the lower triangle of the lower left */
/*             symmetric block Q' of the square-reduced Hamiltonian */
/*             matrix H' in (1), and the N-by-N upper triangular part of */
/*             the submatrix in the columns 2 to N+1 of this array must */
/*             contain the upper triangle of the upper right symmetric */
/*             block G' of the square-reduced Hamiltonian matrix H' */
/*             in (1), as produced by SLICOT Library routine MB04ZD. */
/*             So, if i >= j, then Q'(i,j) is stored in QG(i,j) and */
/*             G'(i,j) is stored in QG(j,i+1). */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             The arrays WR and WI contain the real and imaginary parts, */
/*             respectively, of the N eigenvalues of H' with non-negative */
/*             real part.  The remaining N eigenvalues are the negatives */
/*             of these eigenvalues. */
/*             Eigenvalues are stored in WR and WI in decreasing order of */
/*             magnitude of the real parts, i.e., WR(I) >= WR(I+1). */
/*             (In particular, an eigenvalue closest to the imaginary */
/*              axis is WR(N)+WI(N)i.) */
/*             In addition, eigenvalues with zero real part are sorted in */
/*             decreasing order of magnitude of imaginary parts.  Note */
/*             that non-real eigenvalues with non-zero real part appear */
/*             in complex conjugate pairs, but eigenvalues with zero real */
/*             part do not, in general, appear in complex conjugate */
/*             pairs. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= MAX(1,N*(N+1)). */
/*             For good performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, then the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO =  i, i <= N, then LAPACK subroutine DHSEQR */
/*                   failed to converge while computing the i-th */
/*                   eigenvalue. */

/*     METHOD */

/*     The routine forms the upper Hessenberg matrix A'' in (2) and calls */
/*     LAPACK subroutines to calculate its eigenvalues.  The eigenvalues */
/*     of H' are the square roots of the eigenvalues of A''. */

/*     REFERENCES */

/*     [1] Van Loan, C. F. */
/*         A Symplectic Method for Approximating All the Eigenvalues of */
/*         a Hamiltonian Matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     [2] Byers, R. */
/*         Hamiltonian and Symplectic Algorithms for the Algebraic */
/*         Riccati Equation. */
/*         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983. */

/*     [3] Benner, P., Byers, R., and Barth, E. */
/*         Fortran 77 Subroutines for Computing the Eigenvalues of */
/*         Hamiltonian Matrices. I: The Square-Reduced Method. */
/*         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires (32/3)*N**3 + O(N**2) floating point */
/*     operations. */
/*     Eigenvalues computed by this subroutine are exact eigenvalues */
/*     of a perturbed Hamiltonian matrix  H' + E  where */

/*                 || E || <= c sqrt(eps) || H' ||, */

/*     c is a modest constant depending on the dimension N and eps is the */
/*     machine precision. Moreover, if the norm of H' and an eigenvalue */
/*     are of roughly the same magnitude, the computed eigenvalue is */
/*     essentially as accurate as the computed eigenvalue obtained by */
/*     traditional methods. See [1] or [2]. */

/*     CONTRIBUTOR */

/*     P. Benner, Universitaet Bremen, Germany, and */
/*     R. Byers, University of Kansas, Lawrence, USA. */
/*     Aug. 1998, routine DHAEVS. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998, SLICOT Library version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2002, */
/*     May 2009. */

/*     KEYWORDS */

/*     Eigenvalues, (square-reduced) Hamiltonian matrix, symplectic */
/*     similarity transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    *info = 0;
    n2 = *n * *n;
    scale = lsame_(jobscl, "S", (ftnlen)1, (ftnlen)1);
    if (! (scale || lsame_(jobscl, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldqg < max(1,*n)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = n2 + *n;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -10;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03SD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

    chunk = (*ldwork - n2) / *n;
    block = min(chunk,*n) > 1;
    blas3 = chunk >= *n;

    if (blas3) {
	jwork = n2 + 1;
    } else {
	jwork = 1;
    }
/*                             2 */
/*     Form the matrix A'' = A'  + G'Q'. */

    dlacpy_("Lower", n, n, &qg[qg_offset], ldqg, &dwork[jwork], n, (ftnlen)5);
    ma02ed_("Lower", n, &dwork[jwork], n, (ftnlen)5);

    if (blas3) {

/*        Use BLAS 3 calculation. */

	dsymm_("Left", "Upper", n, n, &c_b9, &qg[(qg_dim1 << 1) + 1], ldqg, &
		dwork[jwork], n, &c_b10, &dwork[1], n, (ftnlen)4, (ftnlen)5);

    } else if (block) {
	jw = n2 + 1;

/*        Use BLAS 3 for as many columns of Q' as possible. */

	i__1 = *n;
	i__2 = chunk;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = *n - j + 1;
	    bl = min(i__3,chunk);
	    dsymm_("Left", "Upper", n, &bl, &c_b9, &qg[(qg_dim1 << 1) + 1], 
		    ldqg, &dwork[*n * (j - 1) + 1], n, &c_b10, &dwork[jw], n, 
		    (ftnlen)4, (ftnlen)5);
	    dlacpy_("Full", n, &bl, &dwork[jw], n, &dwork[*n * (j - 1) + 1], 
		    n, (ftnlen)4);
/* L10: */
	}

    } else {

/*        Use BLAS 2 calculation. */

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    dsymv_("Upper", n, &c_b9, &qg[(qg_dim1 << 1) + 1], ldqg, &dwork[*
		    n * (j - 1) + 1], &c__1, &c_b10, &wr[1], &c__1, (ftnlen)5)
		    ;
	    dcopy_(n, &wr[1], &c__1, &dwork[*n * (j - 1) + 1], &c__1);
/* L20: */
	}

    }

    dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b9, &a[a_offset], lda, &
	    a[a_offset], lda, &c_b9, &dwork[1], n, (ftnlen)11, (ftnlen)11);
    if (scale && *n > 2) {
	i__2 = *n - 2;
	i__1 = *n - 2;
	dlaset_("Lower", &i__2, &i__1, &c_b10, &c_b10, &dwork[3], n, (ftnlen)
		5);
    }
/*                               2 */
/*     Find the eigenvalues of A' + G'Q'. */

    dgebal_(jobscl, n, &dwork[1], n, &ilo, &ihi, &dwork[n2 + 1], &ignore, (
	    ftnlen)1);
    dhseqr_("Eigenvalues", "NoSchurVectors", n, &ilo, &ihi, &dwork[1], n, &wr[
	    1], &wi[1], dummy, &c__1, &dwork[n2 + 1], n, info, (ftnlen)11, (
	    ftnlen)14);
    if (*info == 0) {

/*        Eigenvalues of H' are the square roots of those computed above. */

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = wr[i__];
	    y = wi[i__];
	    ma01ad_(&x, &y, &wr[i__], &wi[i__]);
/* L30: */
	}

/*        Sort eigenvalues into decreasing order by real part and, for */
/*        eigenvalues with zero real part only, decreasing order of */
/*        imaginary part.  (This simple bubble sort preserves the */
/*        relative order of eigenvalues with equal but nonzero real part. */
/*        This ensures that complex conjugate pairs remain */
/*        together.) */

	sorted = FALSE_;

	for (m = *n; m >= 1; --m) {
	    if (sorted) {
		goto L60;
	    }
	    sorted = TRUE_;

	    i__2 = m - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (wr[i__] < wr[i__ + 1] || wr[i__] == 0. && wr[i__ + 1] == 
			0. && wi[i__] < wi[i__ + 1]) {
		    swap = wr[i__];
		    wr[i__] = wr[i__ + 1];
		    wr[i__ + 1] = swap;
		    swap = wi[i__];
		    wi[i__] = wi[i__ + 1];
		    wi[i__ + 1] = swap;

		    sorted = FALSE_;

		}
/* L40: */
	    }

/* L50: */
	}

L60:

	;
    }

    dwork[1] = (doublereal) (n2 << 1);
    return 0;
/*     *** Last line of MB03SD *** */
} /* mb03sd_ */

