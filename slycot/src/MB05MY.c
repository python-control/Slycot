/* MB05MY.f -- translated by f2c (version 20100827).
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
static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__8 = 8;
static integer c__4 = 4;

/* Subroutine */ int mb05my_(char *balanc, integer *n, doublereal *a, integer 
	*lda, doublereal *wr, doublereal *wi, doublereal *r__, integer *ldr, 
	doublereal *q, integer *ldq, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen balanc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, r_dim1, r_offset, i__1, i__2, 
	    i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, ihi, ilo;
    static doublereal dum[1], eps;
    static integer ibal, maxb;
    static doublereal anrm;
    static integer ierr, itau, nout;
    static logical scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebal_(
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    static logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical select[1];
    static doublereal bignum;
    extern /* Subroutine */ int dorghr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dhseqr_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtrevc_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen);
    static integer hsdwor, minwrk, maxwrk;
    static doublereal smlnum;


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

/*     To compute, for an N-by-N real nonsymmetric matrix A, the */
/*     orthogonal matrix Q reducing it to real Schur form T, the */
/*     eigenvalues, and the right eigenvectors of T. */

/*     The right eigenvector r(j) of T satisfies */
/*                      T * r(j) = lambda(j) * r(j) */
/*     where lambda(j) is its eigenvalue. */

/*     The matrix of right eigenvectors R is upper triangular, by */
/*     construction. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Indicates how the input matrix should be diagonally scaled */
/*             to improve the conditioning of its eigenvalues as follows: */
/*             = 'N':  Do not diagonally scale; */
/*             = 'S':  Diagonally scale the matrix, i.e. replace A by */
/*                     D*A*D**(-1), where D is a diagonal matrix chosen */
/*                     to make the rows and columns of A more equal in */
/*                     norm. Do not permute. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the leading N-by-N upper quasi-triangular part of */
/*             this array contains the real Schur canonical form of A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,N). */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             WR and WI contain the real and imaginary parts, */
/*             respectively, of the computed eigenvalues. Complex */
/*             conjugate pairs of eigenvalues appear consecutively */
/*             with the eigenvalue having the positive imaginary part */
/*             first. */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the matrix of right eigenvectors R, in the same */
/*             order as their eigenvalues. The real and imaginary parts */
/*             of a complex eigenvector corresponding to an eigenvalue */
/*             with positive imaginary part are stored in consecutive */
/*             columns. (The corresponding conjugate eigenvector is not */
/*             stored.) The eigenvectors are not backward transformed */
/*             for balancing (when BALANC = 'S'). */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= max(1,N). */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             The leading N-by-N part of this array contains the */
/*             orthogonal matrix Q which has reduced A to real Schur */
/*             form. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK. */
/*             If BALANC = 'S', DWORK(2),...,DWORK(N+1) return the */
/*             scaling factors used for balancing. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= max(1,4*N). */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, the QR algorithm failed to compute all */
/*                   the eigenvalues, and no eigenvectors have been */
/*                   computed; elements i+1:N of WR and WI contain */
/*                   eigenvalues which have converged. */

/*     METHOD */

/*     This routine uses the QR algorithm to obtain the real Schur form */
/*     T of matrix A. Then, the right eigenvectors of T are computed, */
/*     but they are not backtransformed into the eigenvectors of A. */
/*     MB05MY is a modification of the LAPACK driver routine DGEEV. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05AY. */

/*     REVISIONS */

/*     V. Sima, April 25, 2003, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvector decomposition, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
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

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wr;
    --wi;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    scale = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1);
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || scale)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldr < max(1,*n)) {
	*info = -8;
    } else if (*ldq < max(1,*n)) {
	*info = -10;
    }

/*     Compute workspace. */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV. */
/*       HSDWOR refers to the workspace preferred by DHSEQR, as */
/*       calculated below. HSDWOR is computed assuming ILO=1 and IHI=N, */
/*       the worst case.) */

    minwrk = 1;
    if (*info == 0 && *ldwork >= 1) {
	maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "DGEHRD", " ", n, &c__1, n, &
		c__0, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
	i__1 = 1, i__2 = *n << 2;
	minwrk = max(i__1,i__2);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "DORGHR", 
		" ", n, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)1);
	maxwrk = max(i__1,i__2);
/* Computing MAX */
	i__1 = ilaenv_(&c__8, "DHSEQR", "SV", n, &c__1, n, &c_n1, (ftnlen)6, (
		ftnlen)2);
	maxb = max(i__1,2);
/* Computing MIN */
/* Computing MAX */
	i__3 = 2, i__4 = ilaenv_(&c__4, "DHSEQR", "SV", n, &c__1, n, &c_n1, (
		ftnlen)6, (ftnlen)2);
	i__1 = min(maxb,*n), i__2 = max(i__3,i__4);
	k = min(i__1,i__2);
/* Computing MAX */
	i__1 = k * (k + 2), i__2 = *n << 1;
	hsdwor = max(i__1,i__2);
/* Computing MAX */
	i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1,i__2), i__2 = *n + 
		hsdwor;
	maxwrk = max(i__1,i__2);
/* Computing MAX */
	i__1 = maxwrk, i__2 = *n << 2;
	maxwrk = max(i__1,i__2);
	dwork[1] = (doublereal) maxwrk;
    }
    if (*ldwork < minwrk) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB05MY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants. */

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1);
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM]. */

    anrm = dlange_("M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
    scalea = FALSE_;
    if (anrm > 0. && anrm < smlnum) {
	scalea = TRUE_;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = TRUE_;
	cscale = bignum;
    }
    if (scalea) {
	dlascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &
		ierr, (ftnlen)1);
    }

/*     Balance the matrix, if requested. (Permutation is not possible.) */
/*     (Workspace: need N) */

    ibal = 1;
    dgebal_(balanc, n, &a[a_offset], lda, &ilo, &ihi, &dwork[ibal], &ierr, (
	    ftnlen)1);

/*     Reduce to upper Hessenberg form. */
/*     (Workspace: need 3*N, prefer 2*N+N*NB) */

    itau = ibal + *n;
    jwork = itau + *n;
    i__1 = *ldwork - jwork + 1;
    dgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &dwork[itau], &dwork[jwork], &
	    i__1, &ierr);

/*     Compute right eigenvectors of T. */
/*     Copy Householder vectors to Q. */

    dlacpy_("Lower", n, n, &a[a_offset], lda, &q[q_offset], ldq, (ftnlen)5);

/*     Generate orthogonal matrix in Q. */
/*     (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

    i__1 = *ldwork - jwork + 1;
    dorghr_(n, &ilo, &ihi, &q[q_offset], ldq, &dwork[itau], &dwork[jwork], &
	    i__1, &ierr);

/*     Perform QR iteration, accumulating Schur vectors in Q. */
/*     (Workspace: need N+1, prefer N+HSDWOR (see comments) ) */

    jwork = itau;
    i__1 = *ldwork - jwork + 1;
    dhseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &q[
	    q_offset], ldq, &dwork[jwork], &i__1, info, (ftnlen)1, (ftnlen)1);

/*     If INFO > 0 from DHSEQR, then quit. */

    if (*info > 0) {
	goto L10;
    }

/*     Compute right eigenvectors of T in R. */
/*     (Workspace: need 4*N) */

    dtrevc_("Right", "All", select, n, &a[a_offset], lda, dum, &c__1, &r__[
	    r_offset], ldr, n, &nout, &dwork[jwork], &ierr, (ftnlen)5, (
	    ftnlen)3);

/*     Undo scaling if necessary. */

L10:
    if (scalea) {
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 
		1], &i__2, &ierr, (ftnlen)1);
	i__1 = *n - *info;
/* Computing MAX */
	i__3 = *n - *info;
	i__2 = max(i__3,1);
	dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 
		1], &i__2, &ierr, (ftnlen)1);
	if (*info > 0) {
	    i__1 = ilo - 1;
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], 
		    n, &ierr, (ftnlen)1);
	    i__1 = ilo - 1;
	    dlascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], 
		    n, &ierr, (ftnlen)1);
	}
    }

    if (scale) {
	for (k = *n; k >= 1; --k) {
	    dwork[k + 1] = dwork[k];
/* L20: */
	}
    }
    dwork[1] = (doublereal) maxwrk;

    return 0;
/* *** Last line of MB05MY *** */
} /* mb05my_ */

