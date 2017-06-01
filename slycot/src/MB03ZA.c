/* MB03ZA.f -- translated by f2c (version 20100827).
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
static doublereal c_b21 = 1.;
static integer c__4 = 4;
static logical c_true = TRUE_;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b479 = -1.;
static integer c__12 = 12;

/* Subroutine */ int mb03za_(char *compc, char *compu, char *compv, char *
	compw, char *which, logical *select, integer *n, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *u1, integer *ldu1, doublereal *u2, integer *ldu2, 
	doublereal *v1, integer *ldv1, doublereal *v2, integer *ldv2, 
	doublereal *w, integer *ldw, doublereal *wr, doublereal *wi, integer *
	m, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	compc_len, ftnlen compu_len, ftnlen compv_len, ftnlen compw_len, 
	ftnlen which_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, v1_dim1, v1_offset, v2_dim1, 
	    v2_offset, w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, l;
    static doublereal q[16]	/* was [4][4] */, t[16]	/* was [4][4] */, z__[
	    16]	/* was [4][4] */;
    static integer nb, mm, ks, pw, nbf, nbl;
    static doublereal dw12[12];
    static integer len, pwc, pwd, pos, here;
    static logical pair;
    static integer idum[1], ierr;
    static logical ldum[1];
    static integer pwck, ifst, pwdl;
    static doublereal temp;
    static logical swap;
    static integer ilst;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgees_(char *, char *, L_fp, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen), dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), mb03wa_(
	    logical *, logical *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen), lfdum_();
    static logical wantc;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal winew[4];
    static logical initw, wantu, wantv, wantw;
    static doublereal wrnew[4];
    static logical cmpall;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical selnew[4];
    static integer nbnext;
    extern /* Subroutine */ int dtrsen_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer wrkmin;


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

/*     1. To compute, for a given matrix pair (A,B) in periodic Schur */
/*        form, orthogonal matrices Ur and Vr so that */

/*            T           [ A11  A12 ]     T           [ B11  B12 ] */
/*          Vr * A * Ur = [          ],  Ur * B * Vr = [          ], (1) */
/*                        [  0   A22 ]                 [  0   B22 ] */

/*        is in periodic Schur form, and the eigenvalues of A11*B11 */
/*        form a selected cluster of eigenvalues. */

/*     2. To compute an orthogonal matrix W so that */

/*                   T  [  0  -A11 ]       [  R11   R12 ] */
/*                  W * [          ] * W = [            ],           (2) */
/*                      [ B11   0  ]       [   0    R22 ] */

/*        where the eigenvalues of R11 and -R22 coincide and have */
/*        positive real part. */

/*     Optionally, the matrix C is overwritten by Ur'*C*Vr. */

/*     All eigenvalues of A11*B11 must either be complex or real and */
/*     negative. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPC   CHARACTER*1 */
/*             = 'U':  update the matrix C; */
/*             = 'N':  do not update C. */

/*     COMPU   CHARACTER*1 */
/*             = 'U':  update the matrices U1 and U2; */
/*             = 'N':  do not update U1 and U2. */
/*             See the description of U1 and U2. */

/*     COMPV   CHARACTER*1 */
/*             = 'U':  update the matrices V1 and V2; */
/*             = 'N':  do not update V1 and V2. */
/*             See the description of V1 and V2. */

/*     COMPW   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix W as follows: */
/*             = 'N':  the matrix W is not required; */
/*             = 'I':  W is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix W is returned; */
/*             = 'V':  W must contain an orthogonal matrix Q on entry, */
/*                     and the product Q*W is returned. */

/*     WHICH   CHARACTER*1 */
/*             = 'A':  select all eigenvalues, this effectively means */
/*                     that Ur and Vr are identity matrices and A11 = A, */
/*                     B11 = B; */
/*             = 'S':  select a cluster of eigenvalues specified by */
/*                     SELECT. */

/*     SELECT  LOGICAL array, dimension (N) */
/*             If WHICH = 'S', then SELECT specifies the eigenvalues of */
/*             A*B in the selected cluster. To select a real eigenvalue */
/*             w(j), SELECT(j) must be set to .TRUE.. To select a complex */
/*             conjugate pair of eigenvalues w(j) and w(j+1), */
/*             corresponding to a 2-by-2 diagonal block in A, both */
/*             SELECT(j) and SELECT(j+1) must be set to .TRUE.; a complex */
/*             conjugate pair of eigenvalues must be either both included */
/*             in the cluster or both excluded. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A of the matrix */
/*             pair (A,B) in periodic Schur form. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the matrix R22 in (2). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix B of the matrix pair */
/*             (A,B) in periodic Schur form. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, if COMPC = 'U', the leading N-by-N part of this */
/*             array must contain a general matrix C. */
/*             On exit, if COMPC = 'U', the leading N-by-N part of this */
/*             array contains the updated matrix Ur'*C*Vr. */
/*             If COMPC = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= 1. */
/*             LDC >= N,  if COMPC = 'U' and WHICH = 'S'. */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain U1, the (1,1) */
/*             block of an orthogonal symplectic matrix */
/*             U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains U1*Ur. */
/*             If COMPU = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= 1. */
/*             LDU1 >= N,  if COMPU = 'U' and WHICH = 'S'. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain U2, the (1,2) */
/*             block of an orthogonal symplectic matrix */
/*             U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains U2*Ur. */
/*             If COMPU = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= 1. */
/*             LDU2 >= N,  if COMPU = 'U' and WHICH = 'S'. */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain V1, the (1,1) */
/*             block of an orthogonal symplectic matrix */
/*             V = [ V1, V2; -V2, V1 ]. */
/*             On exit, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains V1*Vr. */
/*             If COMPV = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= 1. */
/*             LDV1 >= N,  if COMPV = 'U' and WHICH = 'S'. */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,N) */
/*             On entry, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain V2, the (1,2) */
/*             block of an orthogonal symplectic matrix */
/*             V = [ V1, V2; -V2, V1 ]. */
/*             On exit, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains V2*Vr. */
/*             If COMPV = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= 1. */
/*             LDV2 >= N,  if COMPV = 'U' and WHICH = 'S'. */

/*     W       (input/output) DOUBLE PRECISION array, dimension (LDW,2*M) */
/*             On entry, if COMPW = 'V', then the leading 2*M-by-2*M part */
/*             of this array must contain a matrix W. */
/*             If COMPW = 'I', then W need not be set on entry, W is set */
/*             to the identity matrix. */
/*             On exit, if COMPW = 'I' or 'V' the leading 2*M-by-2*M part */
/*             of this array is post-multiplied by the transformation */
/*             matrix that produced (2). */
/*             If COMPW = 'N', this array is not referenced. */

/*     LDW     INTEGER */
/*             The leading dimension of the array W.  LDW >= 1. */
/*             LDW >= 2*M,  if COMPW = 'I' or COMPW = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (M) */
/*     WI      (output) DOUBLE PRECISION array, dimension (M) */
/*             The real and imaginary parts, respectively, of the */
/*             eigenvalues of R22. The eigenvalues are stored in the same */
/*             order as on the diagonal of R22, with */
/*             WR(i) = R22(i,i) and, if R22(i:i+1,i:i+1) is a 2-by-2 */
/*             diagonal block, WI(i) > 0 and WI(i+1) = -WI(i). */
/*             In exact arithmetic, these eigenvalue are the positive */
/*             square roots of the selected eigenvalues of the product */
/*             A*B. However, if an eigenvalue is sufficiently */
/*             ill-conditioned, then its value may differ significantly. */

/*     M       (output) INTEGER */
/*             The number of selected eigenvalues. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = -28,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, 4*N, 8*M ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  reordering of the product A*B in Step 1 failed */
/*                   because some eigenvalues are too close to separate; */
/*             = 2:  reordering of some submatrix in Step 2 failed */
/*                   because some eigenvalues are too close to separate; */
/*             = 3:  the QR algorithm failed to compute the Schur form */
/*                   of some submatrix in Step 2; */
/*             = 4:  the condition that all eigenvalues of A11*B11 must */
/*                   either be complex or real and negative is */
/*                   numerically violated. */

/*     METHOD */

/*     Step 1 is performed using a reordering technique analogous to the */
/*     LAPACK routine DTGSEN for reordering matrix pencils [1,2]. Step 2 */
/*     is an implementation of Algorithm 2 in [3]. It requires O(M*N*N) */
/*     floating point operations. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A direct method for reordering eigenvalues in the generalized */
/*         real Schur form of a regular matrix pair (A,B), in M.S. Moonen */
/*         et al (eds), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., 1993, pp. 195-218. */

/*     [2] Kagstrom, B. and Poromaa P.: */
/*         Computing eigenspaces with specified eigenvalues of a regular */
/*         matrix pair (A, B) and condition estimation: Theory, */
/*         algorithms and software, Numer. Algorithms, 1996, vol. 12, */
/*         pp. 369-407. */

/*     [3] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix,  J. Comput. Appl. Math., 86, */
/*         pp. 17-43, 1997. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLABMX). */

/*     KEYWORDS */

/*     Hamiltonian matrix, invariant subspace. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and check input parameters */

    /* Parameter adjustments */
    --select;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    v1_dim1 = *ldv1;
    v1_offset = 1 + v1_dim1;
    v1 -= v1_offset;
    v2_dim1 = *ldv2;
    v2_offset = 1 + v2_dim1;
    v2 -= v2_offset;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    wantc = lsame_(compc, "U", (ftnlen)1, (ftnlen)1);
    wantu = lsame_(compu, "U", (ftnlen)1, (ftnlen)1);
    wantv = lsame_(compv, "U", (ftnlen)1, (ftnlen)1);
    initw = lsame_(compw, "I", (ftnlen)1, (ftnlen)1);
    wantw = initw || lsame_(compw, "V", (ftnlen)1, (ftnlen)1);
    cmpall = lsame_(which, "A", (ftnlen)1, (ftnlen)1);
/* Computing MAX */
    i__1 = 1, i__2 = *n << 2;
    wrkmin = max(i__1,i__2);

    *info = 0;
    if (! wantc && ! lsame_(compc, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! wantu && ! lsame_(compu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! wantv && ! lsame_(compv, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (! wantw && ! lsame_(compw, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -4;
    } else if (! cmpall && ! lsame_(which, "S", (ftnlen)1, (ftnlen)1)) {
	*info = -5;
    } else {
	if (cmpall) {
	    *m = *n;
	} else {

/*           Set M to the dimension of the specified invariant subspace. */

	    *m = 0;
	    pair = FALSE_;
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (pair) {
		    pair = FALSE_;
		} else {
		    if (k < *n) {
			if (a[k + 1 + k * a_dim1] == 0.) {
			    if (select[k]) {
				++(*m);
			    }
			} else {
			    pair = TRUE_;
			    if (select[k] || select[k + 1]) {
				*m += 2;
			    }
			}
		    } else {
			if (select[*n]) {
			    ++(*m);
			}
		    }
		}
/* L10: */
	    }
	}

/*        Compute workspace requirements. */

/* Computing MAX */
	i__1 = wrkmin, i__2 = *m << 3;
	wrkmin = max(i__1,i__2);

	if (*n < 0) {
	    *info = -7;
	} else if (*lda < max(1,*n)) {
	    *info = -9;
	} else if (*ldb < max(1,*n)) {
	    *info = -11;
	} else if (*ldc < 1 || wantc && ! cmpall && *ldc < *n) {
	    *info = -13;
	} else if (*ldu1 < 1 || wantu && ! cmpall && *ldu1 < *n) {
	    *info = -15;
	} else if (*ldu2 < 1 || wantu && ! cmpall && *ldu2 < *n) {
	    *info = -17;
	} else if (*ldv1 < 1 || wantv && ! cmpall && *ldv1 < *n) {
	    *info = -19;
	} else if (*ldv2 < 1 || wantv && ! cmpall && *ldv2 < *n) {
	    *info = -21;
	} else if (*ldw < 1 || wantw && *ldw < *m << 1) {
	    *info = -23;
	} else if (*ldwork < wrkmin) {
	    *info = -28;
	    dwork[1] = (doublereal) wrkmin;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03ZA", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Jump immediately to Step 2, if all eigenvalues are requested. */

    if (cmpall) {
	goto L50;
    }

/*     Step 1: Collect the selected blocks at the top-left corner of A*B. */

    ks = 0;
    pair = FALSE_;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (pair) {
	    pair = FALSE_;
	} else {
	    swap = select[k];
	    if (k < *n) {
		if (a[k + 1 + k * a_dim1] != 0.) {
		    pair = TRUE_;
		    swap = swap || select[k + 1];
		}
	    }

	    if (pair) {
		nbf = 2;
	    } else {
		nbf = 1;
	    }

	    if (swap) {
		++ks;
		ifst = k;

/*              Swap the K-th block to position KS. */

		ilst = ks;
		nbl = 1;
		if (ilst > 1) {
		    if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
			--ilst;
			nbl = 2;
		    }
		}

		if (ilst == ifst) {
		    goto L30;
		}

		here = ifst;
L20:

/*              Swap block with next one above. */

		if (nbf == 1 || nbf == 2) {

/*                 Current block either 1-by-1 or 2-by-2. */

		    nbnext = 1;
		    if (here >= 3) {
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
			    nbnext = 2;
			}
		    }
		    pos = here - nbnext;
		    nb = nbnext + nbf;
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
			    ftnlen)3);
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
			    ftnlen)3);

		    mb03wa_(&c_true, &c_true, &nbnext, &nbf, &a[pos + pos * 
			    a_dim1], lda, &b[pos + pos * b_dim1], ldb, q, &
			    c__4, z__, &c__4, &ierr);

		    if (ierr != 0) {
			dwork[1] = (doublereal) wrkmin;
			*info = 1;
			return 0;
		    }

/*                 Update rest of A. */

		    if (pos > 1) {
			i__2 = pos - 1;
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &a[pos * a_dim1 + 1], lda, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			i__2 = pos - 1;
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				a_dim1 + 1], lda, (ftnlen)3);
		    }
		    if (pos + nb <= *n) {
			i__2 = *n - pos - nb + 1;
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, q, &c__4, &a[pos + (pos + nb) * a_dim1]
				, lda, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				ftnlen)12);
			i__2 = *n - pos - nb + 1;
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos + (
				pos + nb) * a_dim1], lda, (ftnlen)3);
		    }

/*                 Update rest of B. */

		    if (pos > 1) {
			i__2 = pos - 1;
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &b[pos * b_dim1 + 1], ldb, q, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			i__2 = pos - 1;
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				b_dim1 + 1], ldb, (ftnlen)3);
		    }
		    if (pos + nb <= *n) {
			i__2 = *n - pos - nb + 1;
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, z__, &c__4, &b[pos + (pos + nb) * 
				b_dim1], ldb, &c_b20, &dwork[1], &nb, (ftnlen)
				9, (ftnlen)12);
			i__2 = *n - pos - nb + 1;
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos + (
				pos + nb) * b_dim1], ldb, (ftnlen)3);
		    }

/*                 Update C. */

		    if (wantc) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &c__[pos * c_dim1 + 1], ldc, q, &c__4, 
				&c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				c_dim1 + 1], ldc, (ftnlen)3);
			dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				c_b21, z__, &c__4, &c__[pos + c_dim1], ldc, &
				c_b20, &dwork[1], &nb, (ftnlen)9, (ftnlen)12);
			dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				c_dim1], ldc, (ftnlen)3);
		    }

/*                 Update U. */

		    if (wantu) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u1[pos * u1_dim1 + 1], ldu1, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				u1_dim1 + 1], ldu1, (ftnlen)3);
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u2[pos * u2_dim1 + 1], ldu2, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				u2_dim1 + 1], ldu2, (ftnlen)3);
		    }

/*                 Update V. */

		    if (wantv) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v1[pos * v1_dim1 + 1], ldv1, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
			dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				v1_dim1 + 1], ldv1, (ftnlen)3);
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v2[pos * v2_dim1 + 1], ldv2, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
			dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				v2_dim1 + 1], ldv2, (ftnlen)3);
		    }

		    here -= nbnext;

/*                 Test if 2-by-2 block breaks into two 1-by-1 blocks. */

		    if (nbf == 2) {
			if (a[here + 1 + here * a_dim1] == 0.) {
			    nbf = 3;
			}
		    }

		} else {

/*                 Current block consists of two 1 by 1 blocks each of */
/*                 which must be swapped individually. */

		    nbnext = 1;
		    if (here >= 3) {
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
			    nbnext = 2;
			}
		    }
		    pos = here - nbnext;
		    nb = nbnext + 1;
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
			    ftnlen)3);
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
			    ftnlen)3);

		    mb03wa_(&c_true, &c_true, &nbnext, &c__1, &a[pos + pos * 
			    a_dim1], lda, &b[pos + pos * b_dim1], ldb, q, &
			    c__4, z__, &c__4, &ierr);

		    if (ierr != 0) {
			dwork[1] = (doublereal) wrkmin;
			*info = 1;
			return 0;
		    }

/*                 Update rest of A. */

		    if (pos > 1) {
			i__2 = pos - 1;
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &a[pos * a_dim1 + 1], lda, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			i__2 = pos - 1;
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				a_dim1 + 1], lda, (ftnlen)3);
		    }
		    if (pos + nb <= *n) {
			i__2 = *n - pos - nb + 1;
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, q, &c__4, &a[pos + (pos + nb) * a_dim1]
				, lda, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				ftnlen)12);
			i__2 = *n - pos - nb + 1;
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos + (
				pos + nb) * a_dim1], lda, (ftnlen)3);
		    }

/*                 Update rest of B. */

		    if (pos > 1) {
			i__2 = pos - 1;
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &b[pos * b_dim1 + 1], ldb, q, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			i__2 = pos - 1;
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				b_dim1 + 1], ldb, (ftnlen)3);
		    }
		    if (pos + nb <= *n) {
			i__2 = *n - pos - nb + 1;
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, z__, &c__4, &b[pos + (pos + nb) * 
				b_dim1], ldb, &c_b20, &dwork[1], &nb, (ftnlen)
				9, (ftnlen)12);
			i__2 = *n - pos - nb + 1;
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos + (
				pos + nb) * b_dim1], ldb, (ftnlen)3);
		    }

/*                 Update C. */

		    if (wantc) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &c__[pos * c_dim1 + 1], ldc, q, &c__4, 
				&c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				c_dim1 + 1], ldc, (ftnlen)3);
			dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				c_b21, z__, &c__4, &c__[pos + c_dim1], ldc, &
				c_b20, &dwork[1], &nb, (ftnlen)9, (ftnlen)12);
			dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				c_dim1], ldc, (ftnlen)3);
		    }

/*                 Update U. */

		    if (wantu) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u1[pos * u1_dim1 + 1], ldu1, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				u1_dim1 + 1], ldu1, (ftnlen)3);
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u2[pos * u2_dim1 + 1], ldu2, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
			dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				u2_dim1 + 1], ldu2, (ftnlen)3);
		    }

/*                 Update V. */

		    if (wantv) {
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v1[pos * v1_dim1 + 1], ldv1, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
			dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				v1_dim1 + 1], ldv1, (ftnlen)3);
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v2[pos * v2_dim1 + 1], ldv2, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
			dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				v2_dim1 + 1], ldv2, (ftnlen)3);
		    }

		    if (nbnext == 1) {

/*                    Swap two 1-by-1 blocks. */

			pos = here;
			nb = nbnext + 1;
			dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
				ftnlen)3);
			dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
				ftnlen)3);

			mb03wa_(&c_true, &c_true, &nbnext, &c__1, &a[pos + 
				pos * a_dim1], lda, &b[pos + pos * b_dim1], 
				ldb, q, &c__4, z__, &c__4, &ierr);

			if (ierr != 0) {
			    dwork[1] = (doublereal) wrkmin;
			    *info = 1;
			    return 0;
			}

/*                    Update rest of A. */

			if (pos > 1) {
			    i__2 = pos - 1;
			    dgemm_("No Transpose", "No Transpose", &i__2, &nb,
				     &nb, &c_b21, &a[pos * a_dim1 + 1], lda, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    i__2 = pos - 1;
			    dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				    a_dim1 + 1], lda, (ftnlen)3);
			}
			if (pos + nb <= *n) {
			    i__2 = *n - pos - nb + 1;
			    dgemm_("Transpose", "No Transpose", &nb, &i__2, &
				    nb, &c_b21, q, &c__4, &a[pos + (pos + nb) 
				    * a_dim1], lda, &c_b20, &dwork[1], &nb, (
				    ftnlen)9, (ftnlen)12);
			    i__2 = *n - pos - nb + 1;
			    dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos 
				    + (pos + nb) * a_dim1], lda, (ftnlen)3);
			}

/*                    Update rest of B. */

			if (pos > 1) {
			    i__2 = pos - 1;
			    dgemm_("No Transpose", "No Transpose", &i__2, &nb,
				     &nb, &c_b21, &b[pos * b_dim1 + 1], ldb, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    i__2 = pos - 1;
			    dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				    b_dim1 + 1], ldb, (ftnlen)3);
			}
			if (pos + nb <= *n) {
			    i__2 = *n - pos - nb + 1;
			    dgemm_("Transpose", "No Transpose", &nb, &i__2, &
				    nb, &c_b21, z__, &c__4, &b[pos + (pos + 
				    nb) * b_dim1], ldb, &c_b20, &dwork[1], &
				    nb, (ftnlen)9, (ftnlen)12);
			    i__2 = *n - pos - nb + 1;
			    dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos 
				    + (pos + nb) * b_dim1], ldb, (ftnlen)3);
			}

/*                    Update C. */

			if (wantc) {
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &c__[pos * c_dim1 + 1], ldc, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				    c_dim1 + 1], ldc, (ftnlen)3);
			    dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				    c_b21, z__, &c__4, &c__[pos + c_dim1], 
				    ldc, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				    ftnlen)12);
			    dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				    c_dim1], ldc, (ftnlen)3);
			}

/*                    Update U. */

			if (wantu) {
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &u1[pos * u1_dim1 + 1], ldu1, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				    u1_dim1 + 1], ldu1, (ftnlen)3);
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &u2[pos * u2_dim1 + 1], ldu2, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				    u2_dim1 + 1], ldu2, (ftnlen)3);
			}

/*                    Update V. */

			if (wantv) {
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &v1[pos * v1_dim1 + 1], ldv1, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				    v1_dim1 + 1], ldv1, (ftnlen)3);
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &v2[pos * v2_dim1 + 1], ldv2, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
			    dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				    v2_dim1 + 1], ldv2, (ftnlen)3);
			}

			--here;
		    } else {

/*                    Recompute NBNEXT in case 2-by-2 split. */

			if (a[here + (here - 1) * a_dim1] == 0.) {
			    nbnext = 1;
			}

			if (nbnext == 2) {

/*                       2-by-2 block did not split. */

			    pos = here - 1;
			    nb = 3;
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

			    if (ierr != 0) {
				dwork[1] = (doublereal) wrkmin;
				*info = 1;
				return 0;
			    }

/*                       Update rest of A. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
			    }

/*                       Update rest of B. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
			    }

/*                       Update C. */

			    if (wantc) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
			    }

/*                       Update U. */

			    if (wantu) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
			    }

/*                       Update V. */

			    if (wantv) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
			    }

			    here += -2;
			} else {

/*                       2-by-2 block did split. */

			    pos = here;
			    nb = 2;
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

			    if (ierr != 0) {
				dwork[1] = (doublereal) wrkmin;
				*info = 1;
				return 0;
			    }

/*                       Update rest of A. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
			    }

/*                       Update rest of B. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
			    }

/*                       Update C. */

			    if (wantc) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
			    }

/*                       Update U. */

			    if (wantu) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
			    }

/*                       Update V. */

			    if (wantv) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
			    }

			    pos = here - 1;
			    nb = 2;
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

			    if (ierr != 0) {
				dwork[1] = (doublereal) wrkmin;
				*info = 1;
				return 0;
			    }

/*                       Update rest of A. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
			    }

/*                       Update rest of B. */

			    if (pos > 1) {
				i__2 = pos - 1;
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
				i__2 = pos - 1;
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
			    }
			    if (pos + nb <= *n) {
				i__2 = *n - pos - nb + 1;
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
				i__2 = *n - pos - nb + 1;
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
			    }

/*                       Update C. */

			    if (wantc) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
			    }

/*                       Update U. */

			    if (wantu) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
			    }

/*                       Update V. */

			    if (wantv) {
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
			    }

			    here += -2;
			}
		    }
		}

		if (here > ilst) {
		    goto L20;
		}

L30:
		if (pair) {
		    ++ks;
		}
	    }
	}
/* L40: */
    }

L50:

/*     Step 2: Compute an ordered Schur decomposition of */
/*             [ 0, -A11; B11, 0 ]. */

    if (initw) {
	i__1 = *m << 1;
	i__2 = *m << 1;
	dlaset_("All", &i__1, &i__2, &c_b20, &c_b21, &w[w_offset], ldw, (
		ftnlen)3);
    }
    pwc = 1;
    pwd = pwc + (*m << 1);
    pw = pwd + (*m << 1);
    pair = FALSE_;
    nb = 1;

    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	if (pair) {
	    pair = FALSE_;
	    nb = 1;
	} else {
	    if (k < *n) {
		if (a[k + 1 + k * a_dim1] != 0.) {
		    pair = TRUE_;
		    nb = 2;
		}
	    }
	    pwck = pwc + (k - 1 << 1);
	    pwdl = pwd + (k - 1 << 1);
	    i__2 = *m - k + 1;
	    dlaset_("All", &nb, &i__2, &c_b20, &c_b20, &dwork[pwck], &c__2, (
		    ftnlen)3);
	    i__2 = *m - k + 1;
	    dlacpy_("All", &nb, &i__2, &a[k + k * a_dim1], lda, &dwork[pwdl], 
		    &c__2, (ftnlen)3);
	    i__2 = *m - k + 1;
	    dlaset_("All", &nb, &i__2, &c_b20, &c_b20, &a[k + k * a_dim1], 
		    lda, (ftnlen)3);

	    l = k;

/*           WHILE L >= 1 DO */

L60:

	    if (k == l) {

/*                 Annihilate B(k,k). */

		nbl = nb;
		i__2 = nb + nbl;
		i__3 = nb + nbl;
		dlaset_("All", &i__2, &i__3, &c_b20, &c_b20, t, &c__4, (
			ftnlen)3);
		dlacpy_("Upper", &nbl, &nbl, &b[l + l * b_dim1], ldb, &t[nb], 
			&c__4, (ftnlen)5);
		if (nb == 1) {
		    dwork[pwdl] = -dwork[pwdl];
		} else {
		    i__2 = nb << 1;
		    dscal_(&i__2, &c_b479, &dwork[pwdl], &c__1);
		}
		dlacpy_("All", &nb, &nb, &dwork[pwdl], &c__2, &t[(nb + 1 << 2)
			 - 4], &c__4, (ftnlen)3);
	    } else {

/*                 Annihilate B(l,k). */

		i__2 = nbl + nb;
		i__3 = nbl + nb;
		dlaset_("All", &i__2, &i__3, &c_b20, &c_b20, t, &c__4, (
			ftnlen)3);
		dlacpy_("All", &nbl, &nbl, &a[l + l * a_dim1], lda, t, &c__4, 
			(ftnlen)3);
		dlacpy_("All", &nbl, &nb, &b[l + k * b_dim1], ldb, &t[(nbl + 
			1 << 2) - 4], &c__4, (ftnlen)3);
		dlacpy_("All", &nb, &nb, &dwork[pwck], &c__2, &t[nbl + 1 + (
			nbl + 1 << 2) - 5], &c__4, (ftnlen)3);
		pwdl = pwd + (l - 1 << 1);
	    }

	    i__2 = nb + nbl;
	    dgees_("V", "Not Sorted", (L_fp)lfdum_, &i__2, t, &c__4, &mm, 
		    wrnew, winew, q, &c__4, dw12, &c__12, ldum, &ierr, (
		    ftnlen)1, (ftnlen)10);
	    if (ierr != 0) {
		dwork[1] = (doublereal) wrkmin;
		*info = 3;
		return 0;
	    }

/*              Reorder Schur form. */

	    mm = 0;
	    i__2 = nb + nbl;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (wrnew[i__ - 1] > 0.) {
		    ++mm;
		    selnew[i__ - 1] = TRUE_;
		} else {
		    selnew[i__ - 1] = FALSE_;
		}
/* L70: */
	    }
	    if (mm < nb) {
		dwork[1] = (doublereal) wrkmin;
		*info = 4;
		return 0;
	    }
	    i__2 = nb + nbl;
	    dtrsen_("None", "V", selnew, &i__2, t, &c__4, q, &c__4, wrnew, 
		    winew, &mm, &temp, &temp, dw12, &c__4, idum, &c__1, &ierr,
		     (ftnlen)4, (ftnlen)1);
	    if (ierr != 0) {
		dwork[1] = (doublereal) wrkmin;
		*info = 2;
		return 0;
	    }

/*              Permute Q if necessary. */

	    if (k != l) {
		i__2 = nb + nbl;
		dlacpy_("All", &nbl, &i__2, q, &c__4, &z__[nb], &c__4, (
			ftnlen)3);
		i__2 = nb + nbl;
		dlacpy_("All", &nb, &i__2, &q[nbl], &c__4, z__, &c__4, (
			ftnlen)3);
		i__2 = nb + nbl;
		i__3 = nb + nbl;
		dlacpy_("All", &i__2, &i__3, z__, &c__4, q, &c__4, (ftnlen)3);
	    }

/*              Update "diagonal" blocks. */

	    dlacpy_("All", &nb, &nb, t, &c__4, &dwork[pwck], &c__2, (ftnlen)3)
		    ;
	    dlacpy_("All", &nb, &nbl, &t[(nb + 1 << 2) - 4], &c__4, &dwork[
		    pwdl], &c__2, (ftnlen)3);
	    if (nb == 1) {
		dscal_(&nbl, &c_b479, &dwork[pwdl], &c__2);
	    } else {
		i__2 = nbl << 1;
		dscal_(&i__2, &c_b479, &dwork[pwdl], &c__1);
	    }
	    dlacpy_("All", &nbl, &nbl, &t[nb + 1 + (nb + 1 << 2) - 5], &c__4, 
		    &a[l + l * a_dim1], lda, (ftnlen)3);

/*              Update block columns of A and B. */

	    len = l - 1;
	    if (len > 0) {
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &b[k * b_dim1 + 1], ldb, q, &c__4, &c_b20, &dwork[pw]
			, m, (ftnlen)12, (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &b[k * b_dim1 + 1], ldb, &q[(nb + 1 << 2) - 4],
			 &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (ftnlen)12,
			 (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &a[l * a_dim1 + 1], lda, &q[nb], &c__4, &c_b21,
			 &dwork[pw], m, (ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nb, &dwork[pw], m, &b[k * b_dim1 + 1], 
			ldb, (ftnlen)3);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &a[l * a_dim1 + 1], lda, &q[nb + 1 + (nb + 1 <<
			 2) - 5], &c__4, &c_b21, &dwork[pw + (*m << 1)], m, (
			ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &a[l * 
			a_dim1 + 1], lda, (ftnlen)3);
	    }

/*              Update block column of A. */

	    len = *m - l - nbl + 1;
	    if (len > 0) {
		dgemm_("Transpose", "No Transpose", &nb, &len, &nb, &c_b21, q,
			 &c__4, &dwork[pwdl + (nbl << 1)], &c__2, &c_b20, &
			dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nb, &c_b479, 
			&q[(nb + 1 << 2) - 4], &c__4, &dwork[pwdl + (nbl << 1)
			], &c__2, &c_b20, &dwork[pw + (*m << 1)], &c__2, (
			ftnlen)9, (ftnlen)12);
		dgemm_("Transpose", "No Transpose", &nb, &len, &nbl, &c_b479, 
			&q[nb], &c__4, &a[l + (l + nbl) * a_dim1], lda, &
			c_b21, &dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
		dlacpy_("All", &nb, &len, &dwork[pw], &c__2, &dwork[pwdl + (
			nbl << 1)], &c__2, (ftnlen)3);
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nbl, &c_b21, 
			&q[nb + 1 + (nb + 1 << 2) - 5], &c__4, &a[l + (l + 
			nbl) * a_dim1], lda, &c_b21, &dwork[pw + (*m << 1)], &
			c__2, (ftnlen)9, (ftnlen)12);
		dlacpy_("All", &nbl, &len, &dwork[pw + (*m << 1)], &c__2, &a[
			l + (l + nbl) * a_dim1], lda, (ftnlen)3);
	    }

/*              Update block row of B. */

	    len = *m - k - nb + 1;
	    if (len > 0) {
		dgemm_("Transpose", "No Transpose", &nb, &len, &nb, &c_b21, q,
			 &c__4, &dwork[pwck + (nb << 1)], &c__2, &c_b20, &
			dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nb, &c_b21, &
			q[(nb + 1 << 2) - 4], &c__4, &dwork[pwck + (nb << 1)],
			 &c__2, &c_b20, &dwork[pw + (*m << 1)], &c__2, (
			ftnlen)9, (ftnlen)12);
		dgemm_("Transpose", "No Transpose", &nb, &len, &nbl, &c_b21, &
			q[nb], &c__4, &b[l + (k + nb) * b_dim1], ldb, &c_b21, 
			&dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
		dlacpy_("All", &nb, &len, &dwork[pw], &c__2, &dwork[pwck + (
			nb << 1)], &c__2, (ftnlen)3);
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nbl, &c_b21, 
			&q[nb + 1 + (nb + 1 << 2) - 5], &c__4, &b[l + (k + nb)
			 * b_dim1], ldb, &c_b21, &dwork[pw + (*m << 1)], &
			c__2, (ftnlen)9, (ftnlen)12);
		dlacpy_("All", &nbl, &len, &dwork[pw + (*m << 1)], &c__2, &b[
			l + (k + nb) * b_dim1], ldb, (ftnlen)3);
	    }

/*              Update W. */

	    if (wantw) {
		if (initw) {
		    pos = l;
		    len = k + nb - l;
		} else {
		    pos = 1;
		    len = *m;
		}
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &w[pos + k * w_dim1], ldw, q, &c__4, &c_b20, &dwork[
			pw], m, (ftnlen)12, (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &w[pos + k * w_dim1], ldw, &q[(nb + 1 << 2) - 
			4], &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (ftnlen)
			12, (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &w[pos + (*m + l) * w_dim1], ldw, &q[nb], &
			c__4, &c_b21, &dwork[pw], m, (ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nb, &dwork[pw], m, &w[pos + k * w_dim1],
			 ldw, (ftnlen)3);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &w[pos + (*m + l) * w_dim1], ldw, &q[nb + 1 + (
			nb + 1 << 2) - 5], &c__4, &c_b21, &dwork[pw + (*m << 
			1)], m, (ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &w[pos 
			+ (*m + l) * w_dim1], ldw, (ftnlen)3);

		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &w[*m + pos + k * w_dim1], ldw, q, &c__4, &c_b20, &
			dwork[pw], m, (ftnlen)12, (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &w[*m + pos + k * w_dim1], ldw, &q[(nb + 1 << 
			2) - 4], &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (
			ftnlen)12, (ftnlen)12);
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &w[*m + pos + (*m + l) * w_dim1], ldw, &q[nb], 
			&c__4, &c_b21, &dwork[pw], m, (ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nb, &dwork[pw], m, &w[*m + pos + k * 
			w_dim1], ldw, (ftnlen)3);
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &w[*m + pos + (*m + l) * w_dim1], ldw, &q[nb + 
			1 + (nb + 1 << 2) - 5], &c__4, &c_b21, &dwork[pw + (*
			m << 1)], m, (ftnlen)12, (ftnlen)12);
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &w[*m + 
			pos + (*m + l) * w_dim1], ldw, (ftnlen)3);
	    }

	    --l;
	    nbl = 1;
	    if (l > 1) {
		if (a[l + (l - 1) * a_dim1] != 0.) {
		    nbl = 2;
		    --l;
		}
	    }

/*           END WHILE L >= 1 DO */

	    if (l >= 1) {
		goto L60;
	    }

/*           Copy recomputed eigenvalues. */

	    dcopy_(&nb, wrnew, &c__1, &wr[k], &c__1);
	    dcopy_(&nb, winew, &c__1, &wi[k], &c__1);
	}
/* L80: */
    }
    dwork[1] = (doublereal) wrkmin;
    return 0;
/* *** Last line of MB03ZA *** */
} /* mb03za_ */


logical lfdum_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    logical ret_val;


/*     Void logical function for DGEES. */

    ret_val = FALSE_;
    return ret_val;
/* *** Last line of LFDUM *** */
} /* lfdum_ */

