/* MB04UD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb04ud_(char *jobq, char *jobz, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ranke, 
	integer *istair, doublereal *tol, doublereal *dwork, integer *info, 
	ftnlen jobq_len, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l, lk, km1, nr1, mnk;
    static doublereal emx, tau;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal toler;
    static logical lzero;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobqi, ljobzi, updatq;
    static doublereal emxnrm;
    static logical updatz;


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

/*     To compute orthogonal transformations Q and Z such that the */
/*     transformed pencil Q'(sE-A)Z has the E matrix in column echelon */
/*     form, where E and A are M-by-N matrices. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBQ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the unitary row permutations, as follows: */
/*             = 'N':  Do not form Q; */
/*             = 'I':  Q is initialized to the unit matrix and the */
/*                     unitary row permutation matrix Q is returned; */
/*             = 'U':  The given matrix Q is updated by the unitary */
/*                     row permutations used in the reduction. */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the unitary column transformations, as follows: */
/*             = 'N':  Do not form Z; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     unitary transformation matrix Z is returned; */
/*             = 'U':  The given matrix Z is updated by the unitary */
/*                     transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in the matrices A, E and the order of */
/*             the matrix Q.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in the matrices A, E and the order */
/*             of the matrix Z.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the A matrix of the pencil sE-A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the unitary transformed matrix Q' * A * Z. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the E matrix of the pencil sE-A, to be reduced to */
/*             column echelon form. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the unitary transformed matrix Q' * E * Z, which is in */
/*             column echelon form. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,M). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             On entry, if JOBQ = 'U', then the leading M-by-M part of */
/*             this array must contain a given matrix Q (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading M-by-M part of this array contains the product of */
/*             the input matrix Q and the row permutation matrix used to */
/*             transform the rows of matrix E. */
/*             On exit, if JOBQ = 'I', then the leading M-by-M part of */
/*             this array contains the matrix of accumulated unitary */
/*             row transformations performed. */
/*             If JOBQ = 'N', the array Q is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDQ = 1 and */
/*             declare this array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. If JOBQ = 'U' or */
/*             JOBQ = 'I', LDQ >= MAX(1,M); if JOBQ = 'N', LDQ >= 1. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*) */
/*             On entry, if JOBZ = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix Z (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix Z and the column transformation matrix */
/*             used to transform the columns of matrix E. */
/*             On exit, if JOBZ = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated unitary */
/*             column transformations performed. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'U' or */
/*             JOBZ = 'I', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     RANKE   (output) INTEGER */
/*             The computed rank of the unitary transformed matrix E. */

/*     ISTAIR  (output) INTEGER array, dimension (M) */
/*             This array contains information on the column echelon form */
/*             of the unitary transformed matrix E. Specifically, */
/*             ISTAIR(i) = +j if the first non-zero element E(i,j) */
/*             is a corner point and -j otherwise, for i = 1,2,...,M. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance below which matrix elements are considered */
/*             to be zero. If the user sets TOL to be less than (or */
/*             equal to) zero then the tolerance is taken as */
/*             EPS * MAX(ABS(E(I,J))), where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH), */
/*             I = 1,2,...,M and J = 1,2,...,N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension MAX(M,N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an M-by-N matrix pencil sE-A with E not necessarily regular, */
/*     the routine computes a unitary transformed pencil Q'(sE-A)Z such */
/*     that the matrix Q' * E * Z is in column echelon form (trapezoidal */
/*     form).  Further details can be found in [1]. */

/*     [An M-by-N matrix E with rank(E) = r is said to be in column */
/*     echelon form if the following conditions are satisfied: */
/*     (a) the first (N - r) columns contain only zero elements; and */
/*     (b) if E(i(k),k) is the last nonzero element in column k for */
/*         k = N-r+1,...,N, i.e. E(i(k),k) <> 0 and E(j,k) = 0 for */
/*         j > i(k), then 1 <= i(N-r+1) < i(N-r+2) < ... < i(N) <= M.] */

/*     REFERENCES */

/*     [1] Beelen, Th. and Van Dooren, P. */
/*         An improved algorithm for the computation of Kronecker's */
/*         canonical form of a singular pencil. */
/*         Linear Algebra and Applications, 105, pp. 9-65, 1988. */

/*     NUMERICAL ASPECTS */

/*     It is shown in [1] that the algorithm is numerically backward */
/*     stable. The operations count is proportional to (MAX(M,N))**3. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Based on Release 3.0 routine MB04SD modified by A. Varga, */
/*     German Aerospace Research Establishment, Oberpfaffenhofen, */
/*     Germany, Dec. 1997, to transform also the matrix A. */

/*     REVISIONS */

/*     A. Varga, DLR Oberpfaffenhofen, June 2005. */

/*     KEYWORDS */

/*     Echelon form, orthogonal transformation, staircase form. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --istair;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobqi = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
    updatq = ljobqi || lsame_(jobq, "U", (ftnlen)1, (ftnlen)1);
    ljobzi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
    updatz = ljobzi || lsame_(jobz, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! updatq && ! lsame_(jobq, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! updatz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else if (*lde < max(1,*m)) {
	*info = -8;
    } else if (! updatq && *ldq < 1 || updatq && *ldq < max(1,*m)) {
	*info = -10;
    } else if (! updatz && *ldz < 1 || updatz && *ldz < max(1,*n)) {
	*info = -12;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB04UD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize Q and Z to the identity matrices, if needed. */

    if (ljobqi) {
	dlaset_("Full", m, m, &c_b10, &c_b11, &q[q_offset], ldq, (ftnlen)4);
    }
    if (ljobzi) {
	dlaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz, (ftnlen)4);
    }

/*     Quick return if possible. */

    *ranke = min(*m,*n);

    if (*ranke == 0) {
	return 0;
    }

    toler = *tol;
    if (toler <= 0.) {
	toler = dlamch_("Epsilon", (ftnlen)7) * dlange_("M", m, n, &e[
		e_offset], lde, &dwork[1], (ftnlen)1);
    }

    k = *n;
    lzero = FALSE_;

/*     WHILE ( ( K > 0 ) AND ( NOT a zero submatrix encountered ) ) DO */
L20:
    if (k > 0 && ! lzero) {

/*         Intermediate form of E */

/*                     <--k--><--n-k-> */
/*                l=1 |x....x|       | */
/*                    |      |       | */
/*                    |  Ek  |   X   | */
/*                    |      |       | */
/*            l=m-n+k |x....x|       | */
/*                    ---------------- */
/*                    |      |x ... x|  } */
/*                    |  O   |  x x x|  } */
/*                    |      |    x x|  } n-k */
/*                    |      |      x|  } */

/*        where submatrix Ek = E[1:m-n+k;1:k]. */

/*        Determine row LK in submatrix Ek with largest max-norm */
/*        (starting with row m-n+k). */

	mnk = *m - *n + k;
	emxnrm = 0.;
	lk = mnk;

	for (l = mnk; l >= 1; --l) {
	    emx = (d__1 = e[l + idamax_(&k, &e[l + e_dim1], lde) * e_dim1], 
		    abs(d__1));
	    if (emx > emxnrm) {
		emxnrm = emx;
		lk = l;
	    }
/* L40: */
	}

	if (emxnrm <= toler) {

/*           Set submatrix Ek to zero. */

	    dlaset_("Full", &mnk, &k, &c_b10, &c_b10, &e[e_offset], lde, (
		    ftnlen)4);
	    lzero = TRUE_;
	    *ranke = *n - k;
	} else {

/*           Submatrix Ek is not considered to be identically zero. */
/*           Check whether rows have to be interchanged. */

	    if (lk != mnk) {

/*              Interchange rows lk and m-n+k in whole A- and E-matrix */
/*              and update the row transformation matrix Q, if needed. */
/*              (For Q, the number of elements involved is m.) */

		dswap_(n, &e[lk + e_dim1], lde, &e[mnk + e_dim1], lde);
		dswap_(n, &a[lk + a_dim1], lda, &a[mnk + a_dim1], lda);
		if (updatq) {
		    dswap_(m, &q[lk * q_dim1 + 1], &c__1, &q[mnk * q_dim1 + 1]
			    , &c__1);
		}
	    }

	    km1 = k - 1;

/*           Determine a Householder transformation to annihilate */
/*           E(m-n+k,1:k-1) using E(m-n+k,k) as pivot. */
/*           Apply the transformation to the columns of A and Ek */
/*           (number of elements involved is m for A and m-n+k for Ek). */
/*           Update the column transformation matrix Z, if needed */
/*           (number of elements involved is n). */

	    dlarfg_(&k, &e[mnk + k * e_dim1], &e[mnk + e_dim1], lde, &tau);
	    emx = e[mnk + k * e_dim1];
	    e[mnk + k * e_dim1] = 1.;
	    i__1 = mnk - 1;
	    dlarf_("Right", &i__1, &k, &e[mnk + e_dim1], lde, &tau, &e[
		    e_offset], lde, &dwork[1], (ftnlen)5);
	    dlarf_("Right", m, &k, &e[mnk + e_dim1], lde, &tau, &a[a_offset], 
		    lda, &dwork[1], (ftnlen)5);
	    if (updatz) {
		dlarf_("Right", n, &k, &e[mnk + e_dim1], lde, &tau, &z__[
			z_offset], ldz, &dwork[1], (ftnlen)5);
	    }
	    e[mnk + k * e_dim1] = emx;
	    dlaset_("Full", &c__1, &km1, &c_b10, &c_b10, &e[mnk + e_dim1], 
		    lde, (ftnlen)4);

	    k = km1;
	}
	goto L20;
    }
/*     END WHILE 20 */

/*     Initialise administration staircase form, i.e. */
/*     ISTAIR(i) =  j  if E(i,j) is a nonzero corner point */
/*               = -j  if E(i,j) is on the boundary but is no corner */
/*                     point. */
/*     Thus, */
/*     ISTAIR(m-k) =   n-k           for k=0,...,rank(E)-1 */
/*                 = -(n-rank(E)+1)  for k=rank(E),...,m-1. */

    i__1 = *ranke - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	istair[*m - i__] = *n - i__;
/* L60: */
    }

    nr1 = -(*n - *ranke + 1);

    i__1 = *m - *ranke;
    for (i__ = 1; i__ <= i__1; ++i__) {
	istair[i__] = nr1;
/* L80: */
    }

    return 0;
/* *** Last line of MB04UD *** */
} /* mb04ud_ */

