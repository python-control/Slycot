/* TG01JD.f -- translated by f2c (version 20100827).
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

static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int tg01jd_(char *job, char *systyp, char *equil, integer *n,
	 integer *m, integer *p, doublereal *a, integer *lda, doublereal *e, 
	integer *lde, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, integer *nr, integer *infred, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen systyp_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer m1, n1, p1, nc, lba, lbe, ldm, ldp, ldq, kwa, kwb, kwc;
    static doublereal dum[1];
    static integer kwe, ldz;
    static char jobq[1], jobz[1];
    extern /* Subroutine */ int ma02cd_(integer *, integer *, integer *, 
	    doublereal *, integer *), tg01ad_(char *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ljobc;
    static integer nblck;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical ljobo;
    extern /* Subroutine */ int tg01hx_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer maxmp;
    static logical lsysp, lsysr, lsyss, lspace, fincon, infcon;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical finobs, infobs, ljobir, lequil;


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

/*     To find a reduced (controllable, observable, or irreducible) */
/*     descriptor representation (Ar-lambda*Er,Br,Cr) for an original */
/*     descriptor representation (A-lambda*E,B,C). */
/*     The pencil Ar-lambda*Er is in an upper block Hessenberg form, with */
/*     either Ar or Er upper triangular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to remove the */
/*             uncontrollable and/or unobservable parts as follows: */
/*             = 'I':  Remove both the uncontrollable and unobservable */
/*                     parts to get an irreducible descriptor */
/*                     representation; */
/*             = 'C':  Remove the uncontrollable part only to get a */
/*                     controllable descriptor representation; */
/*             = 'O':  Remove the unobservable part only to get an */
/*                     observable descriptor representation. */

/*     SYSTYP  CHARACTER*1 */
/*             Indicates the type of descriptor system algorithm */
/*             to be applied according to the assumed */
/*             transfer-function matrix as follows: */
/*             = 'R':  Rational transfer-function matrix; */
/*             = 'S':  Proper (standard) transfer-function matrix; */
/*             = 'P':  Polynomial transfer-function matrix. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily scale */
/*             the system (A-lambda*E,B,C) as follows: */
/*             = 'S':  Perform scaling; */
/*             = 'N':  Do not perform scaling. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the descriptor state vector; also the */
/*             order of square matrices A and E, the number of rows of */
/*             matrix B, and the number of columns of matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of descriptor system input vector; also the */
/*             number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of descriptor system output vector; also the */
/*             number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state matrix A. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the reduced order state matrix Ar of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             The matrix Ar is upper triangular if SYSTYP = 'R' or 'P'. */
/*             If SYSTYP = 'S' and JOB = 'C', the matrix [Br Ar] */
/*             is in a controllable staircase form (see TG01HD). */
/*             If SYSTYP = 'S' and JOB = 'I' or 'O', the matrix ( Ar ) */
/*                                                              ( Cr ) */
/*             is in an observable staircase form (see TG01HD). */
/*             The block structure of staircase forms is contained */
/*             in the leading INFRED(7) elements of IWORK. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original descriptor matrix E. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the reduced order descriptor matrix Er of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             The resulting Er has INFRED(6) nonzero sub-diagonals. */
/*             If at least for one k = 1,...,4, INFRED(k) >= 0, then the */
/*             resulting Er is structured being either upper triangular */
/*             or block Hessenberg, in accordance to the last */
/*             performed order reduction phase (see METHOD). */
/*             The block structure of staircase forms is contained */
/*             in the leading INFRED(7) elements of IWORK. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M), */
/*             if JOB = 'C', or (LDB,MAX(M,P)), otherwise. */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input matrix B; if JOB = 'I', */
/*             or JOB = 'O', the remainder of the leading N-by-MAX(M,P) */
/*             part is used as internal workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the reduced input matrix Br of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'C', only the first IWORK(1) rows of B are */
/*             nonzero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original output matrix C; if JOB = 'I', */
/*             or JOB = 'O', the remainder of the leading MAX(M,P)-by-N */
/*             part is used as internal workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix Cr of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'I', or JOB = 'O', only the last IWORK(1) columns */
/*             (in the first NR columns) of C are nonzero. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1          if N = 0. */

/*     NR      (output) INTEGER */
/*             The order of the reduced descriptor representation */
/*             (Ar-lambda*Er,Br,Cr) of an irreducible, controllable, */
/*             or observable realization for the original system, */
/*             depending on JOB = 'I', JOB = 'C', or JOB = 'O', */
/*             respectively. */

/*     INFRED  (output) INTEGER array, dimension 7 */
/*             This array contains information on performed reduction */
/*             and on structure of resulting system matrices as follows: */
/*             INFRED(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction */
/*                            (see METHOD) has been performed. In this */
/*                            case, INFRED(k) is the achieved order */
/*                            reduction in Phase k. */
/*             INFRED(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not */
/*                            performed. */
/*             INFRED(5)  -   the number of nonzero sub-diagonals of A. */
/*             INFRED(6)  -   the number of nonzero sub-diagonals of E. */
/*             INFRED(7)  -   the number of blocks in the resulting */
/*                            staircase form at last performed reduction */
/*                            phase. The block dimensions are contained */
/*                            in the first INFRED(7) elements of IWORK. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determinations when */
/*             transforming (A-lambda*E,B,C). If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             reciprocal condition numbers in rank determinations; a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N+MAX(M,P) */
/*             On exit, if INFO = 0, the leading INFRED(7) elements of */
/*             IWORK contain the orders of the diagonal blocks of */
/*             Ar-lambda*Er. */

/*     DWORK   DOUBLE PRECISION array, dimension LDWORK */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(8*N,2*M,2*P), if EQUIL = 'S'; */
/*             LDWORK >= MAX(N,2*M,2*P),   if EQUIL = 'N'. */
/*             If LDWORK >= MAX(2*N*N+N*M+N*P)+MAX(N,2*M,2*P) then more */
/*             accurate results are to be expected by performing only */
/*             those reductions phases (see METHOD), where effective */
/*             order reduction occurs. This is achieved by saving the */
/*             system matrices before each phase and restoring them if no */
/*             order reduction took place. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the reduction algorithms of [1]. */
/*     The order reduction is performed in 4 phases: */
/*     Phase 1: Eliminate all finite uncontrolable eigenvalues. */
/*              The resulting matrix ( Br Ar ) is in a controllable */
/*              staircase form (see SLICOT Library routine TG01HD), and */
/*              Er is upper triangular. */
/*              This phase is performed if JOB = 'I' or 'C' and */
/*              SYSTYP = 'R' or 'S'. */
/*     Phase 2: Eliminate all infinite and finite nonzero uncontrollable */
/*              eigenvalues. The resulting matrix ( Br Er ) is in a */
/*              controllable staircase form (see TG01HD), and Ar is */
/*              upper triangular. */
/*              This phase is performed if JOB = 'I' or 'C' and */
/*              SYSTYP = 'R' or 'P'. */
/*     Phase 3: Eliminate all finite unobservable eigenvalues. */
/*              The resulting matrix ( Ar ) is in an observable */
/*                                   ( Cr ) */
/*              staircase form (see SLICOT Library routine TG01ID), and */
/*              Er is upper triangular. */
/*              This phase is performed if JOB = 'I' or 'O' and */
/*              SYSTYP = 'R' or 'S'. */
/*     Phase 4: Eliminate all infinite and finite nonzero unobservable */
/*              eigenvalues. The resulting matrix ( Er ) is in an */
/*                                                ( Cr ) */
/*              observable staircase form (see TG01ID), and Ar is */
/*              upper triangular. */
/*              This phase is performed if JOB = 'I' or 'O' and */
/*              SYSTYP = 'R' or 'P'. */

/*     REFERENCES */

/*     [1] A. Varga */
/*         Computation of Irreducible Generalized State-Space */
/*         Realizations. */
/*         Kybernetika, vol. 26, pp. 89-106, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( N**3 )  floating point operations. */

/*     FURTHER COMMENTS */

/*     If the pencil (A-lambda*E) has no zero eigenvalues, then an */
/*     irreducible realization can be computed skipping Phases 1 and 3 */
/*     by using the setting: JOB = 'I' and SYSTYP = 'P'. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     April 1999. Based on the RASP routine RPDSIR. */

/*     REVISIONS */

/*     July 1999, V. Sima, Research Institute for Informatics, Bucharest. */
/*     May 2003, A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     May 2003, March 2004, V. Sima. */

/*     KEYWORDS */

/*     Controllability, irreducible realization, observability, */
/*     orthogonal canonical form, orthogonal transformation. */

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

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --infred;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    maxmp = max(*m,*p);
    n1 = max(1,*n);

/*     Decode JOB. */

    ljobir = lsame_(job, "I", (ftnlen)1, (ftnlen)1);
    ljobc = ljobir || lsame_(job, "C", (ftnlen)1, (ftnlen)1);
    ljobo = ljobir || lsame_(job, "O", (ftnlen)1, (ftnlen)1);

/*     Decode SYSTYP. */

    lsysr = lsame_(systyp, "R", (ftnlen)1, (ftnlen)1);
    lsyss = lsysr || lsame_(systyp, "S", (ftnlen)1, (ftnlen)1);
    lsysp = lsysr || lsame_(systyp, "P", (ftnlen)1, (ftnlen)1);

    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljobc && ! ljobo) {
	*info = -1;
    } else if (! lsyss && ! lsysp) {
	*info = -2;
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*lda < n1) {
	*info = -8;
    } else if (*lde < n1) {
	*info = -10;
    } else if (*ldb < n1) {
	*info = -12;
    } else if (*ldc < 1 || *n > 0 && *ldc < maxmp) {
	*info = -14;
    } else if (*tol >= 1.) {
	*info = -17;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = *n, i__2 = maxmp << 1;
/* Computing MAX */
	i__3 = *n << 3, i__4 = maxmp << 1;
	if (! lequil && *ldwork < max(i__1,i__2) || lequil && *ldwork < max(
		i__3,i__4)) {
	    *info = -20;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TG01JD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    infred[1] = -1;
    infred[2] = -1;
    infred[3] = -1;
    infred[4] = -1;
    infred[5] = 0;
    infred[6] = 0;
    infred[7] = 0;

    if (max(*n,maxmp) == 0) {
	*nr = 0;
	return 0;
    }

    m1 = max(1,*m);
    p1 = max(1,*p);
    ldm = max(*ldc,*m);
    ldp = max(*ldc,*p);

/*     Set controllability/observability determination options. */

    fincon = ljobc && lsyss;
    infcon = ljobc && lsysp;
    finobs = ljobo && lsyss;
    infobs = ljobo && lsysp;

/*     Set large workspace option and determine offsets. */

/* Computing MAX */
    i__1 = *n, i__2 = maxmp << 1;
    lspace = *ldwork >= *n * ((*n << 1) + *m + *p) + max(i__1,i__2);
/* Computing MAX */
    i__1 = *n, i__2 = maxmp << 1;
    kwa = max(i__1,i__2) + 1;
    kwe = kwa + *n * *n;
    kwb = kwe + *n * *n;
    kwc = kwb + *n * *m;

/*     If required, scale the system (A-lambda*E,B,C). */
/*     Workspace: need 8*N. */

    if (lequil) {
	tg01ad_("All", n, n, m, p, &c_b12, &a[a_offset], lda, &e[e_offset], 
		lde, &b[b_offset], ldb, &c__[c_offset], &ldp, &dwork[1], &
		dwork[*n + 1], &dwork[(*n << 1) + 1], info, (ftnlen)3);
    }

    *(unsigned char *)jobq = 'N';
    *(unsigned char *)jobz = 'N';
    ldq = 1;
    ldz = 1;
/* Computing MAX */
    i__1 = 0, i__2 = *n - 1;
    lba = max(i__1,i__2);
    lbe = lba;
    nc = *n;
    *nr = *n;

    if (fincon) {

/*        Phase 1: Eliminate all finite uncontrolable eigenvalues. */

	if (lspace) {

/*           Save system matrices. */

	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, m, &b[b_offset], ldb, &dwork[kwb], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", p, &nc, &c__[c_offset], ldc, &dwork[kwc], &p1, (
		    ftnlen)4);
	}

/*        Perform finite controllability form reduction. */
/*        Workspace: need   MAX(N,2*M). */

	tg01hx_(jobq, jobz, &nc, &nc, m, p, &nc, &lbe, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], &ldp, dum, 
		&ldq, dum, &ldz, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (*nr < nc || ! lspace) {
	    if (nblck > 1) {
		lba = iwork[1] + iwork[2] - 1;
	    } else if (nblck == 1) {
		lba = iwork[1] - 1;
	    } else {
		lba = 0;
	    }
	    lbe = 0;
	    infred[1] = nc - *nr;
	    infred[7] = nblck;
	    nc = *nr;
	} else {

/*           Restore system matrices. */

	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, m, &dwork[kwb], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
	    dlacpy_("Full", p, &nc, &dwork[kwc], &p1, &c__[c_offset], ldc, (
		    ftnlen)4);
	}
    }

    if (infcon) {

/*        Phase 2: Eliminate all infinite and all finite nonzero */
/*                 uncontrolable eigenvalues. */

	if (lspace) {

/*           Save system matrices. */

	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, m, &b[b_offset], ldb, &dwork[kwb], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", p, &nc, &c__[c_offset], ldc, &dwork[kwc], &p1, (
		    ftnlen)4);
	}

/*        Perform infinite controllability form reduction. */
/*        Workspace: need   MAX(N,2*M). */

	tg01hx_(jobq, jobz, &nc, &nc, m, p, &nc, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], &ldp, dum, 
		&ldq, dum, &ldz, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (*nr < nc || ! lspace) {
	    if (nblck > 1) {
		lbe = iwork[1] + iwork[2] - 1;
	    } else if (nblck == 1) {
		lbe = iwork[1] - 1;
	    } else {
		lbe = 0;
	    }
	    lba = 0;
	    infred[2] = nc - *nr;
	    infred[7] = nblck;
	    nc = *nr;
	} else {

/*           Restore system matrices. */

	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, m, &dwork[kwb], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
	    dlacpy_("Full", p, &nc, &dwork[kwc], &p1, &c__[c_offset], ldc, (
		    ftnlen)4);
	}
    }

    if (finobs || infobs) {

/*        Compute the pertransposed dual system exploiting matrix shapes. */

/* Computing MAX */
	i__2 = 0, i__3 = nc - 1;
	i__1 = max(i__2,i__3);
	tb01xd_("Z", &nc, m, p, &lba, &i__1, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, dum, &c__1, info, (ftnlen)1);
/* Computing MAX */
	i__2 = 0, i__3 = nc - 1;
	i__1 = max(i__2,i__3);
	ma02cd_(&nc, &lbe, &i__1, &e[e_offset], lde);
    }

    if (finobs) {

/*        Phase 3: Eliminate all finite unobservable eigenvalues. */

	if (lspace) {

/*           Save system matrices. */

	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, p, &b[b_offset], ldb, &dwork[kwc], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", m, &nc, &c__[c_offset], ldc, &dwork[kwb], &m1, (
		    ftnlen)4);
	}

/*        Perform finite observability form reduction. */
/*        Workspace: need   MAX(N,2*P). */

	tg01hx_(jobz, jobq, &nc, &nc, p, m, &nc, &lbe, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], &ldm, dum, 
		&ldz, dum, &ldq, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (*nr < nc || ! lspace) {
	    if (nblck > 1) {
		lba = iwork[1] + iwork[2] - 1;
	    } else if (nblck == 1) {
		lba = iwork[1] - 1;
	    } else {
		lba = 0;
	    }
	    lbe = 0;
	    infred[3] = nc - *nr;
	    infred[7] = nblck;
	    nc = *nr;
	} else {

/*           Restore system matrices. */

	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, p, &dwork[kwc], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
	    dlacpy_("Full", m, &nc, &dwork[kwb], &m1, &c__[c_offset], ldc, (
		    ftnlen)4);
	}
    }

    if (infobs) {

/*        Phase 4: Eliminate all infinite and all finite nonzero */
/*                 unobservable eigenvalues. */

	if (lspace) {

/*           Save system matrices. */

	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, p, &b[b_offset], ldb, &dwork[kwc], &n1, (
		    ftnlen)4);
	    dlacpy_("Full", m, &nc, &c__[c_offset], ldc, &dwork[kwb], &m1, (
		    ftnlen)4);
	}

/*        Perform infinite observability form reduction. */
/*        Workspace: need   MAX(N,2*P). */

	tg01hx_(jobz, jobq, &nc, &nc, p, m, &nc, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], &ldm, dum, 
		&ldz, dum, &ldq, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (*nr < nc || ! lspace) {
	    if (nblck > 1) {
		lbe = iwork[1] + iwork[2] - 1;
	    } else if (nblck == 1) {
		lbe = iwork[1] - 1;
	    } else {
		lbe = 0;
	    }
	    lba = 0;
	    infred[4] = nc - *nr;
	    infred[7] = nblck;
	    nc = *nr;
	} else {

/*           Restore system matrices. */

	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
	    dlacpy_("Full", &nc, p, &dwork[kwc], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
	    dlacpy_("Full", m, &nc, &dwork[kwb], &m1, &c__[c_offset], ldc, (
		    ftnlen)4);
	}
    }

    if (finobs || infobs) {

/*        Compute the pertransposed dual system exploiting matrix shapes. */

/* Computing MAX */
	i__2 = 0, i__3 = nc - 1;
	i__1 = max(i__2,i__3);
	tb01xd_("Z", &nc, p, m, &lba, &i__1, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, dum, &c__1, info, (ftnlen)1);
/* Computing MAX */
	i__2 = 0, i__3 = nc - 1;
	i__1 = max(i__2,i__3);
	ma02cd_(&nc, &lbe, &i__1, &e[e_offset], lde);
    }

/*     Set structural information on A and E. */

    infred[5] = lba;
    infred[6] = lbe;

    return 0;
/* *** Last line of TG01JD *** */
} /* tg01jd_ */

