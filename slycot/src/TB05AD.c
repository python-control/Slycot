/* TB05AD.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;

/* Subroutine */ int tb05ad_(char *baleig, char *inita, integer *n, integer *
	m, integer *p, doublecomplex *freq, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *rcond, doublecomplex *g, integer *ldg, doublereal *evre, 
	doublereal *evim, doublecomplex *hinvb, integer *ldhinv, integer *
	iwork, doublereal *dwork, integer *ldwork, doublecomplex *zwork, 
	integer *lzwork, integer *info, ftnlen baleig_len, ftnlen inita_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, g_dim1, 
	    g_offset, hinvb_dim1, hinvb_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer ij, jj, jp, igh, low, itau;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb02rz_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen), mb02sz_(integer *, doublecomplex *, integer *
	    , integer *, integer *), dswap_(integer *, doublereal *, integer *
	    , doublereal *, integer *), mb02tz_(char *, integer *, doublereal 
	    *, doublecomplex *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, ftnlen);
    static doublereal hnorm;
    static integer jwork;
    static logical lbalba;
    extern /* Subroutine */ int dgebal_(char *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen);
    static char balanc[1];
    static logical lbalea, lbaleb, lbalec;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    static logical linita;
    extern /* Subroutine */ int dhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dormhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
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

/*     To find the complex frequency response matrix (transfer matrix) */
/*     G(freq) of the state-space representation (A,B,C) given by */
/*                                   -1 */
/*        G(freq) = C * ((freq*I - A)  ) * B */

/*     where A, B and C are real N-by-N, N-by-M and P-by-N matrices */
/*     respectively and freq is a complex scalar. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALEIG  CHARACTER*1 */
/*             Determines whether the user wishes to balance matrix A */
/*             and/or compute its eigenvalues and/or estimate the */
/*             condition number of the problem as follows: */
/*             = 'N':  The matrix A should not be balanced and neither */
/*                     the eigenvalues of A nor the condition number */
/*                     estimate of the problem are to be calculated; */
/*             = 'C':  The matrix A should not be balanced and only an */
/*                     estimate of the condition number of the problem */
/*                     is to be calculated; */
/*             = 'B' or 'E' and INITA = 'G':  The matrix A is to be */
/*                     balanced and its eigenvalues calculated; */
/*             = 'A' and INITA = 'G':  The matrix A is to be balanced, */
/*                     and its eigenvalues and an estimate of the */
/*                     condition number of the problem are to be */
/*                     calculated. */

/*     INITA   CHARACTER*1 */
/*             Specifies whether or not the matrix A is already in upper */
/*             Hessenberg form as follows: */
/*             = 'G':  The matrix A is a general matrix; */
/*             = 'H':  The matrix A is in upper Hessenberg form and */
/*                     neither balancing nor the eigenvalues of A are */
/*                     required. */
/*             INITA must be set to 'G' for the first call to the */
/*             routine, unless the matrix A is already in upper */
/*             Hessenberg form and neither balancing nor the eigenvalues */
/*             of A are required. Thereafter, it must be set to 'H' for */
/*             all subsequent calls. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of states, i.e. the order of the state */
/*             transition matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of inputs, i.e. the number of columns in the */
/*             matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of outputs, i.e. the number of rows in the */
/*             matrix C.  P >= 0. */

/*     FREQ    (input) COMPLEX*16 */
/*             The frequency freq at which the frequency response matrix */
/*             (transfer matrix) is to be evaluated. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A. */
/*             If INITA = 'G', then, on exit, the leading N-by-N part of */
/*             this array contains an upper Hessenberg matrix similar to */
/*             (via an orthogonal matrix consisting of a sequence of */
/*             Householder transformations) the original state transition */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             If INITA = 'G', then, on exit, the leading N-by-M part of */
/*             this array contains the product of the transpose of the */
/*             orthogonal transformation matrix used to reduce A to upper */
/*             Hessenberg form and the original input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             If INITA = 'G', then, on exit, the leading P-by-N part of */
/*             this array contains the product of the original output/ */
/*             state matrix C and the orthogonal transformation matrix */
/*             used to reduce A to upper Hessenberg form. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If BALEIG = 'C' or BALEIG = 'A', then RCOND contains an */
/*             estimate of the reciprocal of the condition number of */
/*             matrix H with respect to inversion (see METHOD). */

/*     G       (output) COMPLEX*16 array, dimension (LDG,M) */
/*             The leading P-by-M part of this array contains the */
/*             frequency response matrix G(freq). */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,P). */

/*     EVRE,   (output) DOUBLE PRECISION arrays, dimension (N) */
/*     EVIM    If INITA = 'G' and BALEIG = 'B' or 'E' or BALEIG = 'A', */
/*             then these arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the matrix A. */
/*             Otherwise, these arrays are not referenced. */

/*     HINVB   (output) COMPLEX*16 array, dimension (LDHINV,M) */
/*             The leading N-by-M part of this array contains the */
/*                      -1 */
/*             product H  B. */

/*     LDHINV  INTEGER */
/*             The leading dimension of array HINVB.  LDHINV >= MAX(1,N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N - 1 + MAX(N,M,P)), */
/*                       if INITA = 'G' and BALEIG = 'N', or 'B', or 'E'; */
/*             LDWORK >= MAX(1, N + MAX(N,M-1,P-1)), */
/*                       if INITA = 'G' and BALEIG = 'C', or 'A'; */
/*             LDWORK >= MAX(1, 2*N), */
/*                       if INITA = 'H' and BALEIG = 'C', or 'A'; */
/*             LDWORK >= 1, otherwise. */
/*             For optimum performance when INITA = 'G' LDWORK should be */
/*             larger. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX(1,N*N+2*N), if BALEIG = 'C', or 'A'; */
/*             LZWORK >= MAX(1,N*N),     otherwise. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if more than 30*N iterations are required to */
/*                   isolate all the eigenvalues of the matrix A; the */
/*                   computations are continued; */
/*             = 2:  if either FREQ is too near to an eigenvalue of the */
/*                   matrix A, or RCOND is less than EPS, where EPS is */
/*                   the machine  precision (see LAPACK Library routine */
/*                   DLAMCH). */

/*     METHOD */

/*     The matrix A is first balanced (if BALEIG = 'B' or 'E', or */
/*     BALEIG = 'A') and then reduced to upper Hessenberg form; the same */
/*     transformations are applied to the matrix B and the matrix C. */
/*     The complex Hessenberg matrix  H = (freq*I - A) is then used */
/*                       -1 */
/*     to solve for C * H  * B. */

/*     Depending on the input values of parameters BALEIG and INITA, */
/*     the eigenvalues of matrix A and the condition number of */
/*     matrix H with respect to inversion are also calculated. */

/*     REFERENCES */

/*     [1] Laub, A.J. */
/*         Efficient Calculation of Frequency Response Matrices from */
/*         State-Space Models. */
/*         ACM TOMS, 12, pp. 26-33, 1986. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01FD by A.J.Laub, University of */
/*     Southern California, Los Angeles, CA 90089, United States of */
/*     America, June 1982. */

/*     REVISIONS */

/*     V. Sima, February 22, 1998 (changed the name of TB01RD). */
/*     V. Sima, February 12, 1999, August 7, 2003. */
/*     A. Markovski, Technical University of Sofia, September 30, 2003. */
/*     V. Sima, October 1, 2003. */

/*     KEYWORDS */

/*     Frequency response, Hessenberg form, matrix algebra, input output */
/*     description, multivariable system, orthogonal transformation, */
/*     similarity transformation, state-space representation, transfer */
/*     matrix. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --evre;
    --evim;
    hinvb_dim1 = *ldhinv;
    hinvb_offset = 1 + hinvb_dim1;
    hinvb -= hinvb_offset;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    *info = 0;
    lbalec = lsame_(baleig, "C", (ftnlen)1, (ftnlen)1);
    lbaleb = lsame_(baleig, "B", (ftnlen)1, (ftnlen)1) || lsame_(baleig, 
	    "E", (ftnlen)1, (ftnlen)1);
    lbalea = lsame_(baleig, "A", (ftnlen)1, (ftnlen)1);
    lbalba = lbaleb || lbalea;
    linita = lsame_(inita, "G", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lbalec && ! lbalba && ! lsame_(baleig, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! linita && ! lsame_(inita, "H", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldc < max(1,*p)) {
	*info = -12;
    } else if (*ldg < max(1,*p)) {
	*info = -15;
    } else if (*ldhinv < max(1,*n)) {
	*info = -19;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(*n,*m);
/* Computing MAX */
	i__2 = *n, i__3 = *m - 1, i__2 = max(i__2,i__3), i__3 = *p - 1;
	if (linita && ! lbalec && ! lbalea && *ldwork < *n - 1 + max(i__1,*p) 
		|| linita && (lbalec || lbalea) && *ldwork < *n + max(i__2,
		i__3) || ! linita && (lbalec || lbalea) && *ldwork < *n << 1 
		|| *ldwork < 1) {
	    *info = -22;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 1, i__2 = *n * *n;
	    if ((lbalec || lbalea) && *lzwork < *n * (*n + 2) || *lzwork < 
		    max(i__1,i__2)) {
		*info = -24;
	    }
	}
    }

    if (*info != 0) {

/*        Error return */

	i__1 = -(*info);
	xerbla_("TB05AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (min(*m,*p) > 0) {
	    zlaset_("Full", p, m, &c_b1, &c_b1, &g[g_offset], ldg, (ftnlen)4);
	}
	*rcond = 1.;
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    wrkopt = 1;

    if (linita) {
	*(unsigned char *)balanc = 'N';
	if (lbalba) {
	    *(unsigned char *)balanc = 'B';
	}

/*        Workspace: need N. */

	dgebal_(balanc, n, &a[a_offset], lda, &low, &igh, &dwork[1], info, (
		ftnlen)1);
	if (lbalba) {

/*           Adjust B and C matrices based on information in the */
/*           vector DWORK which describes the balancing of A and is */
/*           defined in the subroutine DGEBAL. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		jj = j;
		if (jj < low || jj > igh) {
		    if (jj < low) {
			jj = low - jj;
		    }
		    jp = (integer) dwork[jj];
		    if (jp != jj) {

/*                    Permute rows of B. */

			if (*m > 0) {
			    dswap_(m, &b[jj + b_dim1], ldb, &b[jp + b_dim1], 
				    ldb);
			}

/*                    Permute columns of C. */

			if (*p > 0) {
			    dswap_(p, &c__[jj * c_dim1 + 1], &c__1, &c__[jp * 
				    c_dim1 + 1], &c__1);
			}
		    }
		}
/* L10: */
	    }

	    if (igh != low) {

		i__1 = igh;
		for (j = low; j <= i__1; ++j) {
		    t = dwork[j];

/*                 Scale rows of permuted B. */

		    if (*m > 0) {
			d__1 = 1. / t;
			dscal_(m, &d__1, &b[j + b_dim1], ldb);
		    }

/*                 Scale columns of permuted C. */

		    if (*p > 0) {
			dscal_(p, &t, &c__[j * c_dim1 + 1], &c__1);
		    }
/* L20: */
		}

	    }
	}

/*        Reduce A to Hessenberg form by orthogonal similarities and */
/*        accumulate the orthogonal transformations into B and C. */
/*        Workspace: need 2*N - 1;  prefer N - 1 + N*NB. */

	itau = 1;
	jwork = itau + *n - 1;
	i__1 = *ldwork - jwork + 1;
	dgehrd_(n, &low, &igh, &a[a_offset], lda, &dwork[itau], &dwork[jwork],
		 &i__1, info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Workspace: need N - 1 + M;  prefer N - 1 + M*NB. */

	i__1 = *ldwork - jwork + 1;
	dormhr_("Left", "Transpose", n, m, &low, &igh, &a[a_offset], lda, &
		dwork[itau], &b[b_offset], ldb, &dwork[jwork], &i__1, info, (
		ftnlen)4, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Workspace: need N - 1 + P;  prefer N - 1 + P*NB. */

	i__1 = *ldwork - jwork + 1;
	dormhr_("Right", "No transpose", p, n, &low, &igh, &a[a_offset], lda, 
		&dwork[itau], &c__[c_offset], ldc, &dwork[jwork], &i__1, info,
		 (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	if (lbalba) {

/*           Temporarily store Hessenberg form of A in array ZWORK. */

	    ij = 0;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ++ij;
		    i__3 = ij;
		    i__4 = i__ + j * a_dim1;
		    z__1.r = a[i__4], z__1.i = 0.;
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L30: */
		}

/* L40: */
	    }

/*           Compute the eigenvalues of A if that option is requested. */
/*           Workspace: need N. */

	    dhseqr_("Eigenvalues", "No Schur", n, &low, &igh, &a[a_offset], 
		    lda, &evre[1], &evim[1], &dwork[1], &c__1, &dwork[1], 
		    ldwork, info, (ftnlen)11, (ftnlen)8);

/*           Restore upper Hessenberg form of A. */

	    ij = 0;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ++ij;
		    i__3 = ij;
		    a[i__ + j * a_dim1] = zwork[i__3].r;
/* L50: */
		}

/* L60: */
	    }

	    if (*info > 0) {

/*              DHSEQR could not evaluate the eigenvalues of A. */

		*info = 1;
	    }
	}
    }

/*     Update  H := (FREQ * I) - A   with appropriate value of FREQ. */

    ij = 0;
    jj = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++ij;
	    i__3 = ij;
	    i__4 = i__ + j * a_dim1;
	    z__2.r = a[i__4], z__2.i = 0.;
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L70: */
	}

	i__2 = jj;
	i__3 = jj;
	z__1.r = freq->r + zwork[i__3].r, z__1.i = freq->i + zwork[i__3].i;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
	jj = jj + *n + 1;
/* L80: */
    }

    if (lbalec || lbalea) {

/*        Efficiently compute the 1-norm of the matrix for condition */
/*        estimation. */

	hnorm = 0.;
	jj = 1;

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 1;
	    t = z_abs(&zwork[jj]) + dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
	    if (j < *n) {
		t += (d__1 = a[j + 1 + j * a_dim1], abs(d__1));
	    }
	    hnorm = max(hnorm,t);
	    jj = jj + *n + 1;
/* L90: */
	}

    }

/*     Factor the complex Hessenberg matrix. */

    mb02sz_(n, &zwork[1], n, &iwork[1], info);
    if (*info != 0) {
	*info = 2;
    }

    if (lbalec || lbalea) {

/*        Estimate the condition of the matrix. */

/*        Workspace: need 2*N. */

	mb02tz_("1-norm", n, &hnorm, &zwork[1], n, &iwork[1], rcond, &dwork[1]
		, &zwork[*n * *n + 1], info, (ftnlen)6);
/* Computing MAX */
	i__1 = wrkopt, i__2 = *n << 1;
	wrkopt = max(i__1,i__2);
	if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
	    *info = 2;
	}
    }

    if (*info != 0) {

/*        Error return: Linear system is numerically or exactly singular. */

	return 0;
    }

/*     Compute  (H-INVERSE)*B. */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * hinvb_dim1;
	    i__4 = i__ + j * b_dim1;
	    z__1.r = b[i__4], z__1.i = 0.;
	    hinvb[i__3].r = z__1.r, hinvb[i__3].i = z__1.i;
/* L100: */
	}

/* L110: */
    }

    mb02rz_("No transpose", n, m, &zwork[1], n, &iwork[1], &hinvb[
	    hinvb_offset], ldhinv, info, (ftnlen)12);

/*     Compute  C*(H-INVERSE)*B. */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * g_dim1;
	    g[i__3].r = 0., g[i__3].i = 0.;
/* L120: */
	}

	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {

	    i__3 = *p;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = i__ + j * g_dim1;
		i__5 = i__ + j * g_dim1;
		i__6 = i__ + k * c_dim1;
		z__3.r = c__[i__6], z__3.i = 0.;
		i__7 = k + j * hinvb_dim1;
		z__2.r = z__3.r * hinvb[i__7].r - z__3.i * hinvb[i__7].i, 
			z__2.i = z__3.r * hinvb[i__7].i + z__3.i * hinvb[i__7]
			.r;
		z__1.r = g[i__5].r + z__2.r, z__1.i = g[i__5].i + z__2.i;
		g[i__4].r = z__1.r, g[i__4].i = z__1.i;
/* L130: */
	    }

/* L140: */
	}

/* L150: */
    }

/*     G now contains the desired frequency response matrix. */
/*     Set the optimal workspace. */

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of TB05AD *** */
} /* tb05ad_ */

