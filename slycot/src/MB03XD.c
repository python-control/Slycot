/* MB03XD.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static doublereal c_b30 = 1.;
static doublereal c_b31 = -1.;
static doublereal c_b56 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03xd_(char *balanc, char *job, char *jobu, char *jobv, 
	integer *n, doublereal *a, integer *lda, doublereal *qg, integer *
	ldqg, doublereal *t, integer *ldt, doublereal *u1, integer *ldu1, 
	doublereal *u2, integer *ldu2, doublereal *v1, integer *ldv1, 
	doublereal *v2, integer *ldv2, doublereal *wr, doublereal *wi, 
	integer *ilo, doublereal *scale, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen balanc_len, ftnlen job_len, ftnlen jobu_len, 
	ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, t_dim1, t_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, v1_dim1, v1_offset, v2_dim1, 
	    v2_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, pq, pz;
    static doublereal eps;
    static integer pdw, ilo1, ierr, pcsl;
    static doublereal hnrm, temp;
    static integer pcsr;
    extern /* Subroutine */ int ma01ad_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), mb04dd_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    extern doublereal ma02id_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04qb_(char *, char *, char *, char *, char *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *), mb04tb_(char *,
	     char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen), dgemm_(char 
	    *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer pbeta;
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char uchar[1], vchar[1];
    extern /* Subroutine */ int mb03xp_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);
    static doublereal tempi;
    static logical lperm, wantg;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ptaul;
    static doublereal tempr;
    static integer ptaur;
    static logical wants, wantu, wantv;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    static logical scaleh;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum;
    static integer wrkmin;
    static doublereal smlnum;
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

/*     To compute the eigenvalues of a Hamiltonian matrix, */

/*                   [  A   G  ]         T        T */
/*             H  =  [       T ],   G = G,   Q = Q,                  (1) */
/*                   [  Q  -A  ] */

/*     where A, G and Q are real n-by-n matrices. */

/*     Due to the structure of H all eigenvalues appear in pairs */
/*     (lambda,-lambda). This routine computes the eigenvalues of H */
/*     using an algorithm based on the symplectic URV and the periodic */
/*     Schur decompositions as described in [1], */

/*           T       [  T   G  ] */
/*          U H V =  [       T ],                                    (2) */
/*                   [  0  -S  ] */

/*     where U and V are 2n-by-2n orthogonal symplectic matrices, */
/*     S is in real Schur form and T is upper triangular. */

/*     The algorithm is backward stable and preserves the eigenvalue */
/*     pairings in finite precision arithmetic. */

/*     Optionally, a symplectic balancing transformation to improve the */
/*     conditioning of eigenvalues is computed (see MB04DD). In this */
/*     case, the matrix H in decomposition (2) must be replaced by the */
/*     balanced matrix. */

/*     The SLICOT Library routine MB03ZD can be used to compute invariant */
/*     subspaces of H from the output of this routine. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Indicates how H should be diagonally scaled and/or */
/*             permuted to reduce its norm. */
/*             = 'N': Do not diagonally scale or permute; */
/*             = 'P': Perform symplectic permutations to make the matrix */
/*                    closer to Hamiltonian Schur form. Do not diagonally */
/*                    scale; */
/*             = 'S': Diagonally scale the matrix, i.e., replace A, G and */
/*                    Q by D*A*D**(-1), D*G*D and D**(-1)*Q*D**(-1) where */
/*                    D is a diagonal matrix chosen to make the rows and */
/*                    columns of H more equal in norm. Do not permute; */
/*             = 'B': Both diagonally scale and permute A, G and Q. */
/*             Permuting does not change the norm of H, but scaling does. */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             decomposition (2) or the eigenvalues only, as follows: */
/*             = 'E': compute the eigenvalues only; */
/*             = 'S': compute matrices T and S of (2); */
/*             = 'G': compute matrices T, S and G of (2). */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether or not the user wishes to compute the */
/*             orthogonal symplectic matrix U of (2) as follows: */
/*             = 'N': the matrix U is not computed; */
/*             = 'U': the matrix U is computed. */

/*     JOBV    CHARACTER*1 */
/*             Indicates whether or not the user wishes to compute the */
/*             orthogonal symplectic matrix V of (2) as follows: */
/*             = 'N': the matrix V is not computed; */
/*             = 'V': the matrix V is computed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, this array is overwritten. If JOB = 'S' or */
/*             JOB = 'G', the leading N-by-N part of this array contains */
/*             the matrix S in real Schur form of decomposition (2). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the lower triangular part of the */
/*             matrix Q and in columns 2:N+1 the upper triangular part */
/*             of the matrix G. */
/*             On exit, this array is overwritten. If JOB = 'G', the */
/*             leading N-by-N+1 part of this array contains in columns */
/*             2:N+1 the matrix G of decomposition (2). */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= max(1,N). */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             On exit, if JOB = 'S' or JOB = 'G', the leading N-by-N */
/*             part of this array contains the upper triangular matrix T */
/*             of the decomposition (2). Otherwise, this array is used as */
/*             workspace. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,N). */

/*     U1      (output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On exit, if JOBU = 'U', the leading N-by-N part of this */
/*             array contains the (1,1) block of the orthogonal */
/*             symplectic matrix U of decomposition (2). */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= 1. */
/*             LDU1 >= N,    if JOBU = 'U'. */

/*     U2      (output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On exit, if JOBU = 'U', the leading N-by-N part of this */
/*             array contains the (2,1) block of the orthogonal */
/*             symplectic matrix U of decomposition (2). */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= 1. */
/*             LDU2 >= N,    if JOBU = 'U'. */

/*     V1      (output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On exit, if JOBV = 'V', the leading N-by-N part of this */
/*             array contains the (1,1) block of the orthogonal */
/*             symplectic matrix V of decomposition (2). */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= 1. */
/*             LDV1 >= N,    if JOBV = 'V'. */

/*     V2      (output) DOUBLE PRECISION array, dimension (LDV2,N) */
/*             On exit, if JOBV = 'V', the leading N-by-N part of this */
/*             array contains the (2,1) block of the orthogonal */
/*             symplectic matrix V of decomposition (2). */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= 1. */
/*             LDV2 >= N,    if JOBV = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, the leading N elements of WR and WI contain the */
/*             real and imaginary parts, respectively, of N eigenvalues */
/*             that have nonpositive real part. Complex conjugate pairs */
/*             of eigenvalues with real part not equal to zero will */
/*             appear consecutively with the eigenvalue having the */
/*             positive imaginary part first. For complex conjugate pairs */
/*             of eigenvalues on the imaginary axis only the eigenvalue */
/*             having nonnegative imaginary part will be returned. */

/*     ILO     (output) INTEGER */
/*             ILO is an integer value determined when H was balanced. */
/*             The balanced A(i,j) = 0 if I > J and J = 1,...,ILO-1. */
/*             The balanced Q(i,j) = 0 if J = 1,...,ILO-1 or */
/*             I = 1,...,ILO-1. */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, if SCALE = 'S', the leading N elements of this */
/*             array contain details of the permutation and scaling */
/*             factors applied when balancing H, see MB04DD. */
/*             This array is not referenced if BALANC = 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -25,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  (input) INTEGER */
/*             The dimension of the array DWORK. LDWORK >= max( 1, 8*N ). */
/*             Moreover: */
/*             If JOB = 'E' or 'S' and JOBU = 'N' and JOBV = 'N', */
/*                LDWORK >= 7*N+N*N. */
/*             If JOB = 'G' and JOBU = 'N' and JOBV = 'N', */
/*                LDWORK >= max( 7*N+N*N, 2*N+3*N*N ). */
/*             If JOB = 'G' and JOBU = 'U' and JOBV = 'N', */
/*                LDWORK >= 7*N+2*N*N. */
/*             If JOB = 'G' and JOBU = 'N' and JOBV = 'V', */
/*                LDWORK >= 7*N+2*N*N. */
/*             If JOB = 'G' and JOBU = 'U' and JOBV = 'V', */
/*                LDWORK >= 7*N+N*N. */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO     (output) INTEGER */
/*              = 0:  successful exit; */
/*              < 0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*              > 0:  if INFO = i, the periodic QR algorithm failed to */
/*                    compute all the eigenvalues, elements i+1:N of WR */
/*                    and WI contain eigenvalues which have converged. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. */
/*         Numer. Math., Vol. 78(3), pp. 329-358, 1998. */

/*     [2] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix,  J. Comput. Appl. Math., vol. 86, */
/*         pp. 17-43, 1997. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAESU). */

/*     KEYWORDS */

/*     Eigenvalues, invariant subspace, Hamiltonian matrix. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
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
    --wr;
    --wi;
    --scale;
    --dwork;

    /* Function Body */
    *info = 0;
    lperm = lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (
	    ftnlen)1, (ftnlen)1);
    lscal = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (
	    ftnlen)1, (ftnlen)1);
    wants = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "G", (
	    ftnlen)1, (ftnlen)1);
    wantg = lsame_(job, "G", (ftnlen)1, (ftnlen)1);
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);

    if (wantg) {
	if (wantu) {
	    if (wantv) {
/* Computing MAX */
		i__1 = 1, i__2 = *n * 7 + *n * *n;
		wrkmin = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = 1, i__2 = *n * 7 + (*n << 1) * *n;
		wrkmin = max(i__1,i__2);
	    }
	} else {
	    if (wantv) {
/* Computing MAX */
		i__1 = 1, i__2 = *n * 7 + (*n << 1) * *n;
		wrkmin = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = 1, i__2 = *n * 7 + *n * *n, i__1 = max(i__1,i__2), 
			i__2 = (*n << 1) + *n * 3 * *n;
		wrkmin = max(i__1,i__2);
	    }
	}
    } else {
	if (wantu) {
	    if (wantv) {
/* Computing MAX */
		i__1 = 1, i__2 = *n << 3;
		wrkmin = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = 1, i__2 = *n << 3;
		wrkmin = max(i__1,i__2);
	    }
	} else {
	    if (wantv) {
/* Computing MAX */
		i__1 = 1, i__2 = *n << 3;
		wrkmin = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = 1, i__2 = *n * 7 + *n * *n;
		wrkmin = max(i__1,i__2);
	    }
	}
    }

    wrkopt = wrkmin;

/*     Test the scalar input parameters. */

    if (! lperm && ! lscal && ! lsame_(balanc, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! wants && ! lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! wantu && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (! wantv && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldqg < max(1,*n)) {
	*info = -9;
    } else if (*ldt < max(1,*n)) {
	*info = -11;
    } else if (*ldu1 < 1 || wantu && *ldu1 < *n) {
	*info = -13;
    } else if (*ldu2 < 1 || wantu && *ldu2 < *n) {
	*info = -15;
    } else if (*ldv1 < 1 || wantv && *ldv1 < *n) {
	*info = -17;
    } else if (*ldv2 < 1 || wantv && *ldv2 < *n) {
	*info = -19;
    } else if (*ldwork < wrkmin) {
	*info = -25;
	dwork[1] = (doublereal) wrkmin;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03XD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *ilo = 0;
    if (*n == 0) {
	return 0;
    }

    eps = dlamch_("P", (ftnlen)1);
    smlnum = dlamch_("S", (ftnlen)1);
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;

/*     Scale H if maximal element is outside range [SMLNUM,BIGNUM]. */

    hnrm = ma02id_("Hamiltonian", "MaxElement", n, &a[a_offset], lda, &qg[
	    qg_offset], ldqg, &dwork[1], (ftnlen)11, (ftnlen)10);
    scaleh = FALSE_;
    if (hnrm > 0. && hnrm < smlnum) {
	scaleh = TRUE_;
	cscale = smlnum;
    } else if (hnrm > bignum) {
	scaleh = TRUE_;
	cscale = bignum;
    }
    if (scaleh) {
	dlascl_("General", &c__0, &c__0, &hnrm, &cscale, n, n, &a[a_offset], 
		lda, &ierr, (ftnlen)7);
	i__1 = *n + 1;
	dlascl_("General", &c__0, &c__0, &hnrm, &cscale, n, &i__1, &qg[
		qg_offset], ldqg, &ierr, (ftnlen)7);
    }

/*     Balance the matrix. */

    mb04dd_(balanc, n, &a[a_offset], lda, &qg[qg_offset], ldqg, ilo, &scale[1]
	    , &ierr, (ftnlen)1);

/*     Copy A to T and multiply A by -1. */

    dlacpy_("All", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)3);
    dlascl_("General", &c__0, &c__0, &c_b30, &c_b31, n, n, &a[a_offset], lda, 
	    &ierr, (ftnlen)7);

/*     --------------------------------------------- */
/*     Step 1: Compute symplectic URV decomposition. */
/*     --------------------------------------------- */

    pcsl = 1;
    pcsr = pcsl + (*n << 1);
    ptaul = pcsr + (*n << 1);
    ptaur = ptaul + *n;
    pdw = ptaur + *n;
    if (! wantu && ! wantv) {

/*         Copy Q and Q' to workspace. */

	pq = pdw;
	pdw += *n * *n;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    k = pq + (*n + 1) * (j - 1);
	    l = k;
	    dwork[k] = qg[j + j * qg_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		++k;
		l += *n;
		temp = qg[i__ + j * qg_dim1];
		dwork[k] = temp;
		dwork[l] = temp;
/* L10: */
	    }
/* L20: */
	}
    } else if (wantu) {

/*         Copy Q and Q' to U2. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    u2[j + j * u2_dim1] = qg[j + j * qg_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		temp = qg[i__ + j * qg_dim1];
		u2[i__ + j * u2_dim1] = temp;
		u2[j + i__ * u2_dim1] = temp;
/* L30: */
	    }
/* L40: */
	}
    } else {

/*         Copy Q and Q' to V2. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    v2[j + j * v2_dim1] = qg[j + j * qg_dim1];
	    i__2 = *n;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		temp = qg[i__ + j * qg_dim1];
		v2[i__ + j * v2_dim1] = temp;
		v2[j + i__ * v2_dim1] = temp;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Transpose G. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    qg[i__ + (j + 1) * qg_dim1] = qg[j + (i__ + 1) * qg_dim1];
/* L70: */
	}
/* L80: */
    }

    if (! wantu && ! wantv) {
	i__1 = *ldwork - pdw + 1;
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &dwork[pq], n, 
		&dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[ptaur], &
		dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
    } else if (wantu) {
	i__1 = *ldwork - pdw + 1;
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &u2[u2_offset],
		 ldu2, &dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[
		ptaur], &dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
    } else {
	i__1 = *ldwork - pdw + 1;
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &v2[v2_offset],
		 ldv2, &dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[
		ptaur], &dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
    wrkopt = max(i__1,i__2);

    if (wantu && ! wantv && ! wantg) {
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
	}
    } else if (! wantu && wantv && ! wantg) {
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlacpy_("Upper", &i__1, &i__2, &v2[(v2_dim1 << 1) + 1], ldv2, &qg[
		    (qg_dim1 << 1) + 1], ldqg, (ftnlen)5);
	}
    } else if (wantu && wantv && ! wantg) {
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &v2[v2_dim1 + 
		    2], ldv2, (ftnlen)5);
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
	}
    } else if (wantu && ! wantv && wantg) {
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &dwork[pdw + *
		    n * *n + *n], &i__3, (ftnlen)5);
	}
    } else if (! wantu && wantv && wantg) {
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    i__3 = *n - 2;
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 3], lda, &dwork[pdw + *
		    n * *n + *n], &i__3, (ftnlen)5);
	}
    } else if (wantu && wantv && wantg) {
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &dwork[pdw + *
		    n], &i__3, (ftnlen)5);
	}
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 3], lda, &v2[v2_dim1 + 
		    3], ldv2, (ftnlen)5);
	}
    }

/*     ---------------------------------------------- */
/*     Step 2:  Compute periodic Schur decomposition. */
/*     ---------------------------------------------- */

    if (*n > 2) {
	i__1 = *n - 2;
	i__2 = *n - 2;
	dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], lda, (
		ftnlen)5);
    }
    if (*n > 1) {
	i__1 = *n - 1;
	i__2 = *n - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], ldt, (
		ftnlen)5);
    }
    if (! wantu && ! wantv) {
	pbeta = 1;
    } else {
	pbeta = pdw;
    }

    if (! wantg) {

/*        Workspace requirements: 2*N (8*N with U or V). */

	pdw = pbeta + *n;
	if (wantu) {
	    *(unsigned char *)uchar = 'I';
	} else {
	    *(unsigned char *)uchar = 'N';
	}
	if (wantv) {
	    *(unsigned char *)vchar = 'I';
	} else {
	    *(unsigned char *)vchar = 'N';
	}
	i__1 = *ldwork - pdw + 1;
	mb03xp_(job, vchar, uchar, n, ilo, n, &a[a_offset], lda, &t[t_offset],
		 ldt, &v1[v1_offset], ldv1, &u1[u1_offset], ldu1, &wr[1], &wi[
		1], &dwork[pbeta], &dwork[pdw], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	if (*info != 0) {
	    goto L90;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);

    } else if (! wantu && ! wantv && wantg) {

/*        Workspace requirements: 3*N*N + 2*N. */

	pq = pbeta + *n;
	pz = pq + *n * *n;
	pdw = pz + *n * *n;
	i__1 = *ldwork - pdw + 1;
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &dwork[pq], n, &dwork[pz], n, &wr[1], &wi[1], 
		&dwork[pbeta], &dwork[pdw], &i__1, info, (ftnlen)5, (ftnlen)4,
		 (ftnlen)4);
	if (*info != 0) {
	    goto L90;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pz], n, &
		qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (ftnlen)
		9, (ftnlen)12);
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &dwork[pq], n, &c_b56, &qg[(qg_dim1 << 1) + 1], ldqg, (
		ftnlen)12, (ftnlen)12);
    } else if (wantu && ! wantv && wantg) {

/*        Workspace requirements: 2*N*N + 7*N. */

	pq = pbeta + *n;
	pdw = pq + *n * *n;
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &dwork[pq], n, &u1[u1_offset], ldu1, &wr[1], &
		wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)], &
		i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
	if (*info != 0) {
	    goto L90;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &t[t_dim1 + 2],
		     ldt, (ftnlen)5);
	}
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &u1[u1_offset], 
		ldu1, &qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (
		ftnlen)9, (ftnlen)12);
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &dwork[pq], n, &c_b56, &qg[(qg_dim1 << 1) + 1], ldqg, (
		ftnlen)12, (ftnlen)12);

    } else if (! wantu && wantv && wantg) {

/*        Workspace requirements: 2*N*N + 7*N */

	pz = pbeta + *n;
	pdw = pz + *n * *n;
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &v1[v1_offset], ldv1, &dwork[pz], n, &wr[1], &
		wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)], &
		i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
	if (*info != 0) {
	    goto L90;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    i__3 = *n - 2;
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &a[a_dim1 + 3],
		     lda, (ftnlen)5);
	}
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pz], n, &
		qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (ftnlen)
		9, (ftnlen)12);
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &v1[v1_offset], ldv1, &c_b56, &qg[(qg_dim1 << 1) + 1], 
		ldqg, (ftnlen)12, (ftnlen)12);

    } else if (wantu && wantv && wantg) {

/*        Workspace requirements: N*N + 7*N. */

	pdw = pbeta + *n;
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &v1[v1_offset], ldv1, &u1[u1_offset], ldu1, &
		wr[1], &wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)
		], &i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
	if (*info != 0) {
	    goto L90;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &t[t_dim1 + 2],
		     ldt, (ftnlen)5);
	}
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    dlacpy_("Lower", &i__1, &i__2, &v2[v2_dim1 + 3], ldv2, &a[a_dim1 
		    + 3], lda, (ftnlen)5);
	}
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &u1[u1_offset], 
		ldu1, &qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (
		ftnlen)9, (ftnlen)12);
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &v1[v1_offset], ldv1, &c_b56, &qg[(qg_dim1 << 1) + 1], 
		ldqg, (ftnlen)12, (ftnlen)12);
    }

L90:

/*     Compute square roots of eigenvalues and rescale. */

    i__1 = *n;
    for (i__ = *info + 1; i__ <= i__1; ++i__) {
	tempr = wr[i__];
	tempi = wi[i__];
	temp = dwork[pbeta + i__ - 1];
	if (temp > 0.) {
	    tempr = -tempr;
	}
	temp = abs(temp);
	if (tempi == 0.) {
	    if (tempr < 0.) {
		wr[i__] = 0.;
		wi[i__] = sqrt(temp) * sqrt(-tempr);
	    } else {
		wr[i__] = -sqrt(temp) * sqrt(tempr);
		wi[i__] = 0.;
	    }
	} else {
	    ma01ad_(&tempr, &tempi, &wr[i__], &wi[i__]);
	    wr[i__] = -wr[i__] * sqrt(temp);
	    if (temp > 0.) {
		wi[i__] *= sqrt(temp);
	    } else {
		wi[i__] = 0.;
	    }
	}
/* L100: */
    }

    if (scaleh) {

/*        Undo scaling. */

	dlascl_("Hessenberg", &c__0, &c__0, &cscale, &hnrm, n, n, &a[a_offset]
		, lda, &ierr, (ftnlen)10);
	dlascl_("Upper", &c__0, &c__0, &cscale, &hnrm, n, n, &t[t_offset], 
		ldt, &ierr, (ftnlen)5);
	if (wantg) {
	    dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, n, &qg[(
		    qg_dim1 << 1) + 1], ldqg, &ierr, (ftnlen)7);
	}
	dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, &c__1, &wr[1], n, 
		&ierr, (ftnlen)7);
	dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, &c__1, &wi[1], n, 
		&ierr, (ftnlen)7);
    }

    if (*info != 0) {
	return 0;
    }

/*     ----------------------------------------------- */
/*     Step 3:  Compute orthogonal symplectic factors. */
/*     ----------------------------------------------- */

/*     Fix CSL and CSR for MB04QB. */

    if (wantu) {
	dscal_(n, &c_b31, &dwork[pcsl + 1], &c__2);
    }
    if (wantv) {
	i__1 = *n - 1;
	dscal_(&i__1, &c_b31, &dwork[pcsr + 1], &c__2);
    }
/* Computing MIN */
    i__1 = *n, i__2 = *ilo + 1;
    ilo1 = min(i__1,i__2);

    if (wantu && ! wantv && ! wantg) {

/*        Workspace requirements: 7*N. */

	pdw = ptaur;
	i__1 = *ldt + 1;
	dcopy_(n, &t[t_dim1 + 1], &i__1, &dwork[pdw], &c__1);
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &t[t_offset], ldt, (
		ftnlen)5);
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
	i__1 = *n - *ilo + 1;
	i__2 = *n - *ilo + 1;
	i__3 = *ldwork - pdw - *n + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &qg[*ilo + *ilo * qg_dim1], 
		ldqg, &t[*ilo + *ilo * t_dim1], ldt, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw + *n], &i__3, &ierr, 
		(ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
	wrkopt = max(i__1,i__2);
	i__1 = *ldt + 1;
	dcopy_(n, &dwork[pdw], &c__1, &t[t_dim1 + 1], &i__1);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
	}

    } else if (! wantu && wantv && ! wantg) {

/*        Workspace requirements: 7*N. */

	pdw = ptaur + *n;
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
	i__2 = 0, i__3 = *n - *ilo;
	i__1 = max(i__2,i__3);
/* Computing MAX */
	i__5 = 0, i__6 = *n - *ilo;
	i__4 = max(i__5,i__6);
	i__7 = *ldwork - pdw + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &qg[ilo1 + *ilo * qg_dim1], ldqg, 
		&qg[*ilo + ilo1 * qg_dim1], ldqg, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw], &i__7, &ierr, (ftnlen)
		12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);

    } else if (wantu && wantv && ! wantg) {

/*        Workspace requirements: 8*N. */

	pdw = ptaur + *n;
	i__1 = *ldt + 1;
	dcopy_(n, &t[t_dim1 + 1], &i__1, &dwork[pdw], &c__1);
	dlacpy_("Lower", n, n, &v2[v2_offset], ldv2, &t[t_offset], ldt, (
		ftnlen)5);
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
	i__2 = 0, i__3 = *n - *ilo;
	i__1 = max(i__2,i__3);
/* Computing MAX */
	i__5 = 0, i__6 = *n - *ilo;
	i__4 = max(i__5,i__6);
	i__7 = *ldwork - pdw - *n + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &qg[ilo1 + *ilo * qg_dim1], ldqg, 
		&u2[*ilo + ilo1 * u2_dim1], ldu2, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw + *n], &i__7, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
	wrkopt = max(i__1,i__2);

	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &qg[qg_offset], ldqg, (
		ftnlen)5);
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
	i__1 = *n - *ilo + 1;
	i__2 = *n - *ilo + 1;
	i__3 = *ldwork - pdw - *n + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&qg[*ilo + *ilo * qg_dim1], ldqg, &u1[*ilo + u1_dim1], ldu1, &
		u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 2], &
		dwork[ptaul + *ilo - 1], &dwork[pdw + *n], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
	wrkopt = max(i__1,i__2);
	i__1 = *ldt + 1;
	dcopy_(n, &dwork[pdw], &c__1, &t[t_dim1 + 1], &i__1);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
	}

    } else if (wantu && ! wantv && wantg) {

/*        Workspace requirements: 6*N + N*N. */

	pq = ptaur;
	pdw = pq + *n * *n;
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &dwork[pq], n, (ftnlen)5)
		;
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
	i__1 = *n - *ilo + 1;
	i__2 = *n - *ilo + 1;
	i__3 = *ldwork - pdw + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&dwork[pq + (*ilo - 1) * (*n + 1)], n, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
	}

    } else if (! wantu && wantv && wantg) {

/*        Workspace requirements: 7*N + N*N. */

	pq = ptaur + *n;
	pdw = pq + *n * *n;
	dlacpy_("Upper", n, n, &v2[v2_offset], ldv2, &dwork[pq], n, (ftnlen)5)
		;
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
	i__2 = 0, i__3 = *n - *ilo;
	i__1 = max(i__2,i__3);
/* Computing MAX */
	i__5 = 0, i__6 = *n - *ilo;
	i__4 = max(i__5,i__6);
	i__7 = *ldwork - pdw - *n + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &a[ilo1 + *ilo * a_dim1], lda, &
		dwork[pq + *ilo * *n + *ilo - 1], n, &v1[ilo1 + v1_dim1], 
		ldv1, &v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 
		2], &dwork[ptaur + *ilo - 1], &dwork[pdw + *n], &i__7, &ierr, 
		(ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], 
		    lda, (ftnlen)5);
	}

    } else if (wantu && wantv && wantg) {

/*        Workspace requirements: 6*N + N*N. */

	pdw = ptaur + *n;
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
	i__2 = 0, i__3 = *n - *ilo;
	i__1 = max(i__2,i__3);
/* Computing MAX */
	i__5 = 0, i__6 = *n - *ilo;
	i__4 = max(i__5,i__6);
	i__7 = *ldwork - pdw + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &a[ilo1 + *ilo * a_dim1], lda, &
		u2[*ilo + ilo1 * u2_dim1], ldu2, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw], &i__7, &ierr, (ftnlen)
		12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);

	pq = ptaur;
	pdw = pq + *n * *n;
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &dwork[pq], n, (ftnlen)5)
		;
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
	i__1 = *n - *ilo + 1;
	i__2 = *n - *ilo + 1;
	i__3 = *ldwork - pdw + 1;
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&dwork[pq + (*ilo - 1) * (*n + 1)], n, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);
	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], 
		    lda, (ftnlen)5);
	}
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
	}
    }

    dwork[1] = (doublereal) wrkopt;
    return 0;
/* *** Last line of MB03XD *** */
} /* mb03xd_ */

