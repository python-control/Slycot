/* TB01UD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static logical c_false = FALSE_;

/* Subroutine */ int tb01ud_(char *jobz, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, integer *ncont, integer *indcon, integer *nblk, 
	doublereal *z__, integer *ldz, doublereal *tau, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, ni, nj, nbl, iqr, rank, itau;
    static doublereal fnrm;
    static integer mcrt, ncrt;
    static doublereal sval[3];
    extern /* Subroutine */ int mb01pd_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical ljobf, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal anorm, bnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical ljobz;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlapmt_(logical *, integer *, 
	    integer *, doublereal *, integer *, integer *), dorgqr_(integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
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

/*     To find a controllable realization for the linear time-invariant */
/*     multi-input system */

/*             dX/dt = A * X + B * U, */
/*                Y  = C * X, */

/*     where A, B, and C are N-by-N, N-by-M, and P-by-N matrices, */
/*     respectively, and A and B are reduced by this routine to */
/*     orthogonal canonical form using (and optionally accumulating) */
/*     orthogonal similarity transformations, which are also applied */
/*     to C.  Specifically, the system (A, B, C) is reduced to the */
/*     triplet (Ac, Bc, Cc), where Ac = Z' * A * Z, Bc = Z' * B, */
/*     Cc = C * Z,  with */

/*             [ Acont     *    ]         [ Bcont ] */
/*        Ac = [                ],   Bc = [       ], */
/*             [   0    Auncont ]         [   0   ] */

/*        and */

/*                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ] */
/*                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ] */
/*                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ] */
/*        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ], */
/*                [  .   .     . .    .     .  ]         [ .  ] */
/*                [  .   .       .    .     .  ]         [ .  ] */
/*                [  0   0   . . .  Ap,p-1 App ]         [ 0  ] */

/*     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and */
/*     p is the controllability index of the pair.  The size of the */
/*     block  Auncont is equal to the dimension of the uncontrollable */
/*     subspace of the pair (A, B). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal similarity transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form Z and do not store the orthogonal */
/*                     transformations; */
/*             = 'F':  Do not form Z, but store the orthogonal */
/*                     transformations in the factored form; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NCONT-by-NCONT part contains the */
/*             upper block Hessenberg state dynamics matrix Acont in Ac, */
/*             given by Z' * A * Z, of a controllable realization for */
/*             the original system. The elements below the first block- */
/*             subdiagonal are set to zero. The leading N-by-N part */
/*             contains the matrix Ac. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the leading NCONT-by-M part of this array */
/*             contains the transformed input matrix Bcont in Bc, given */
/*             by Z' * B, with all elements but the first block set to */
/*             zero. The leading N-by-M part contains the matrix Bc. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output matrix Cc, given by C * Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NCONT   (output) INTEGER */
/*             The order of the controllable state-space representation. */

/*     INDCON  (output) INTEGER */
/*             The controllability index of the controllable part of the */
/*             system representation. */

/*     NBLK    (output) INTEGER array, dimension (N) */
/*             The leading INDCON elements of this array contain the */
/*             the orders of the diagonal blocks of Acont. */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If JOBZ = 'I', then the leading N-by-N part of this */
/*             array contains the matrix of accumulated orthogonal */
/*             similarity transformations which reduces the given system */
/*             to orthogonal canonical form. */
/*             If JOBZ = 'F', the elements below the diagonal, with the */
/*             array TAU, represent the orthogonal transformation matrix */
/*             as a product of elementary reflectors. The transformation */
/*             matrix can then be obtained by calling the LAPACK Library */
/*             routine DORGQR. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'I' or */
/*             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The elements of TAU contain the scalar factors of the */
/*             elementary reflectors used in the reduction of B and A. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS */
/*             is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N, 3*M, P). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Matrix B is first QR-decomposed and the appropriate orthogonal */
/*     similarity transformation applied to the matrix A. Leaving the */
/*     first rank(B) states unchanged, the remaining lower left block */
/*     of A is then QR-decomposed and the new orthogonal matrix, Q1, */
/*     is also applied to the right of A to complete the similarity */
/*     transformation. By continuing in this manner, a completely */
/*     controllable state-space pair (Acont, Bcont) is found for the */
/*     given (A, B), where Acont is upper block Hessenberg with each */
/*     subdiagonal block of full row rank, and Bcont is zero apart from */
/*     its (independent) first rank(B) rows. */
/*     All orthogonal transformations determined in this process are also */
/*     applied to the matrix C, from the right. */
/*     NOTE that the system controllability indices are easily */
/*     calculated from the dimensions of the blocks of Acont. */

/*     REFERENCES */

/*     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D. */
/*         Orthogonal Invariants and Canonical Forms for Linear */
/*         Controllable Systems. */
/*         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981. */

/*     [2] Paige, C.C. */
/*         Properties of numerical algorithms related to computing */
/*         controllablity. */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981. */

/*     [3] Petkov, P.Hr., Konstantinov, M.M., Gu, D.W. and */
/*         Postlethwaite, I. */
/*         Optimal Pole Assignment Design of Linear Multi-Input Systems. */
/*         Leicester University, Report 99-11, May 1996. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     If the system matrices A and B are badly scaled, it would be */
/*     useful to scale them with SLICOT routine TB01ID, before calling */
/*     the routine. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, Nov. 2003. */
/*     A. Varga, DLR Oberpfaffenhofen, March 2002, Nov. 2003. */

/*     KEYWORDS */

/*     Controllability, minimal realization, orthogonal canonical form, */
/*     orthogonal transformation. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --nblk;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --tau;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobf = lsame_(jobz, "F", (ftnlen)1, (ftnlen)1);
    ljobi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
    ljobz = ljobf || ljobi;

/*     Test the input scalar arguments. */

    if (! ljobz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldc < max(1,*p)) {
	*info = -10;
    } else if (! ljobz && *ldz < 1 || ljobz && *ldz < max(1,*n)) {
	*info = -15;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*n), i__2 = *m * 3, i__1 = max(i__1,i__2);
	if (*ldwork < max(i__1,*p)) {
	    *info = -20;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB01UD", &i__1, (ftnlen)6);
	return 0;
    }

    *ncont = 0;
    *indcon = 0;

/*     Calculate the absolute norms of A and B (used for scaling). */

    anorm = dlange_("M", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
    bnorm = dlange_("M", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)1);

/*     Quick return if possible. */

    if (min(*n,*m) == 0 || bnorm == 0.) {
	if (*n > 0) {
	    if (ljobi) {
		dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, (
			ftnlen)4);
	    } else if (ljobf) {
		dlaset_("Full", n, n, &c_b9, &c_b9, &z__[z_offset], ldz, (
			ftnlen)4);
		dlaset_("Full", n, &c__1, &c_b9, &c_b9, &tau[1], n, (ftnlen)4)
			;
	    }
	}
	dwork[1] = 1.;
	return 0;
    }

/*     Scale (if needed) the matrices A and B. */

    mb01pd_("S", "G", n, n, &c__0, &c__0, &anorm, &c__0, &nblk[1], &a[
	    a_offset], lda, info, (ftnlen)1, (ftnlen)1);
    mb01pd_("S", "G", n, m, &c__0, &c__0, &bnorm, &c__0, &nblk[1], &b[
	    b_offset], ldb, info, (ftnlen)1, (ftnlen)1);

/*     Compute the Frobenius norm of [ B  A ] (used for rank estimation). */

    d__1 = dlange_("F", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)1);
    d__2 = dlange_("F", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
    fnrm = dlapy2_(&d__1, &d__2);

    toldef = *tol;
    if (toldef <= 0.) {

/*        Use the default tolerance in controllability determination. */

	toldef = (doublereal) (*n * *n) * dlamch_("EPSILON", (ftnlen)7);
    }

    if (fnrm < toldef) {
	fnrm = 1.;
    }

    wrkopt = 1;
    ni = 0;
    itau = 1;
    ncrt = *n;
    mcrt = *m;
    iqr = 1;

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

L10:

/*        Rank-revealing QR decomposition with column pivoting. */
/*        The calculation is performed in NCRT rows of B starting from */
/*        the row IQR (initialized to 1 and then set to rank(B)+1). */
/*        Workspace: 3*MCRT. */

    mb03oy_(&ncrt, &mcrt, &b[iqr + b_dim1], ldb, &toldef, &fnrm, &rank, sval, 
	    &iwork[1], &tau[itau], &dwork[1], info);

    if (rank != 0) {
	nj = ni;
	ni = *ncont;
	*ncont += rank;
	++(*indcon);
	nblk[*indcon] = rank;

/*           Premultiply and postmultiply the appropriate block row */
/*           and block column of A by Q' and Q, respectively. */
/*           Workspace: need   NCRT; */
/*                      prefer NCRT*NB. */

	dormqr_("Left", "Transpose", &ncrt, &ncrt, &rank, &b[iqr + b_dim1], 
		ldb, &tau[itau], &a[ni + 1 + (ni + 1) * a_dim1], lda, &dwork[
		1], ldwork, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);

/*           Workspace: need   N; */
/*                      prefer N*NB. */

	dormqr_("Right", "No transpose", n, &ncrt, &rank, &b[iqr + b_dim1], 
		ldb, &tau[itau], &a[(ni + 1) * a_dim1 + 1], lda, &dwork[1], 
		ldwork, info, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);

/*           Postmultiply the appropriate block column of C by Q. */
/*           Workspace: need   P; */
/*                      prefer P*NB. */

	dormqr_("Right", "No transpose", p, &ncrt, &rank, &b[iqr + b_dim1], 
		ldb, &tau[itau], &c__[(ni + 1) * c_dim1 + 1], ldc, &dwork[1], 
		ldwork, info, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);

/*           If required, save transformations. */

	if (ljobz && ncrt > 1) {
	    i__1 = ncrt - 1;
/* Computing MIN */
	    i__3 = rank, i__4 = ncrt - 1;
	    i__2 = min(i__3,i__4);
	    dlacpy_("L", &i__1, &i__2, &b[iqr + 1 + b_dim1], ldb, &z__[ni + 2 
		    + itau * z_dim1], ldz, (ftnlen)1);
	}

/*           Zero the subdiagonal elements of the current matrix. */

	if (rank > 1) {
	    i__1 = rank - 1;
	    i__2 = rank - 1;
	    dlaset_("L", &i__1, &i__2, &c_b9, &c_b9, &b[iqr + 1 + b_dim1], 
		    ldb, (ftnlen)1);
	}

/*           Backward permutation of the columns of B or A. */

	if (*indcon == 1) {
	    dlapmt_(&c_false, &rank, m, &b[iqr + b_dim1], ldb, &iwork[1]);
	    iqr = rank + 1;
	} else {
	    i__1 = mcrt;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(&rank, &b[iqr + j * b_dim1], &c__1, &a[ni + 1 + (nj + 
			iwork[j]) * a_dim1], &c__1);
/* L20: */
	    }
	}

	itau += rank;
	if (rank != ncrt) {
	    mcrt = rank;
	    ncrt -= rank;
	    dlacpy_("G", &ncrt, &mcrt, &a[*ncont + 1 + (ni + 1) * a_dim1], 
		    lda, &b[iqr + b_dim1], ldb, (ftnlen)1);
	    dlaset_("G", &ncrt, &mcrt, &c_b9, &c_b9, &a[*ncont + 1 + (ni + 1) 
		    * a_dim1], lda, (ftnlen)1);
	    goto L10;
	}
    }

/*     If required, accumulate transformations. */
/*     Workspace: need N;  prefer N*NB. */

    if (ljobi) {
	i__1 = itau - 1;
	dorgqr_(n, n, &i__1, &z__[z_offset], ldz, &tau[1], &dwork[1], ldwork, 
		info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);
    }

/*     Annihilate the trailing blocks of B. */

    if (iqr <= *n) {
	i__1 = *n - iqr + 1;
	dlaset_("G", &i__1, m, &c_b9, &c_b9, &b[iqr + b_dim1], ldb, (ftnlen)1)
		;
    }

/*     Annihilate the trailing elements of TAU, if JOBZ = 'F'. */

    if (ljobf) {
	i__1 = *n;
	for (j = itau; j <= i__1; ++j) {
	    tau[j] = 0.;
/* L30: */
	}
    }

/*     Undo scaling of A and B. */

    if (*indcon < *n) {
	nbl = *indcon + 1;
	nblk[nbl] = *n - *ncont;
    } else {
	nbl = 0;
    }
    mb01pd_("U", "H", n, n, &c__0, &c__0, &anorm, &nbl, &nblk[1], &a[a_offset]
	    , lda, info, (ftnlen)1, (ftnlen)1);
    mb01pd_("U", "G", &nblk[1], m, &c__0, &c__0, &bnorm, &c__0, &nblk[1], &b[
	    b_offset], ldb, info, (ftnlen)1, (ftnlen)1);

/*     Set optimal workspace dimension. */

    dwork[1] = (doublereal) wrkopt;
    return 0;
/* *** Last line of TB01UD *** */
} /* tb01ud_ */

