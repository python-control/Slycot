/* AB09CX.f -- translated by f2c (version 20100827).
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

static doublereal c_b16 = 1.;
static integer c__1 = 1;
static doublereal c_b49 = -1.;
static doublereal c_b54 = 0.;

/* Subroutine */ int ab09cx_(char *dico, char *ordsel, integer *n, integer *m,
	 integer *p, integer *nr, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *hsv, doublereal *tol1, doublereal *tol2, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *iwarn, integer *
	info, ftnlen dico_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, na, kr, kt, ku, kw, nu, kb1, kb2, kc1, nr1, 
	    kw1, kw2, kti;
    static doublereal skp;
    static integer ldb1, ldb2, ldc1, kc2t, nkr1;
    static doublereal skp2;
    static integer ndim;
    static doublereal atol;
    static integer ierr;
    static doublereal rtol;
    static integer ldc2t;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab04md_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), ab09ax_(char *, char 
	    *, char *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), tb01kd_(char *, char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    mb01sd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, ftnlen), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01wd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static logical discr;
    static integer irank;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer nminr;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer khsvp, khsvp2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dgelsy_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static logical fixord;
    static doublereal srrtol;
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable */
/*     original state-space representation (A,B,C,D) by using the optimal */
/*     Hankel-norm approximation method in conjunction with square-root */
/*     balancing. The state dynamics matrix A of the original system is */
/*     an upper quasi-triangular matrix in real Schur canonical form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of */
/*             the resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(MAX(0,NR-KR+1),NMIN), */
/*             where KR is the multiplicity of the Hankel singular value */
/*             HSV(NR+1), NR is the desired order on entry, and NMIN is */
/*             the order of a minimal realization of the given system; */
/*             NMIN is determined as the number of Hankel singular values */
/*             greater than N*EPS*HNORM(A,B,C), where EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH) and */
/*             HNORM(A,B,C) is the Hankel norm of the system (computed */
/*             in HSV(1)); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A in a real Schur */
/*             canonical form. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, if INFO = 0, the leading NR-by-M part of this */
/*             array contains the input/state matrix Br of the reduced */
/*             order system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, if INFO = 0, the leading P-by-NR part of this */
/*             array contains the state/output matrix Cr of the reduced */
/*             order system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the original input/output matrix D. */
/*             On exit, if INFO = 0, the leading P-by-M part of this */
/*             array contains the input/output matrix Dr of the reduced */
/*             order system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the Hankel singular values of */
/*             the original system ordered decreasingly. HSV(1) is the */
/*             Hankel norm of the system. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(A,B,C), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(A,B,C) is the */
/*             Hankel-norm of the given system (computed in HSV(1)). */
/*             For computing a minimal realization, the recommended */
/*             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH). */
/*             This value is used by default if TOL1 <= 0 on entry. */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the given system. The recommended value is */
/*             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default */
/*             if TOL2 <= 0 on entry. */
/*             If TOL2 > 0, then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = MAX(1,M),   if DICO = 'C'; */
/*             LIWORK = MAX(1,N,M), if DICO = 'D'. */
/*             On exit, if INFO = 0, IWORK(1) contains NMIN, the order of */
/*             the computed minimal realization. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LDW1,LDW2 ), where */
/*             LDW1 = N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2, */
/*             LDW2 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                    MAX( 3*M+1, MIN(N,M)+P ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than the order of a minimal realization of the */
/*                   given system. In this case, the resulting NR is set */
/*                   automatically to a value corresponding to the order */
/*                   of a minimal realization of the system. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 2:  the computation of Hankel singular values failed; */
/*             = 3:  the computation of stable projection failed; */
/*             = 4:  the order of computed stable projection differs */
/*                   from the order of Hankel-norm approximation. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                           (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09CX determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                       (2) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     The optimal Hankel-norm approximation method of [1], based on the */
/*     square-root balancing projection formulas of [2], is employed. */

/*     REFERENCES */

/*     [1] Glover, K. */
/*         All optimal Hankel norm approximation of linear */
/*         multivariable systems and their L-infinity error bounds. */
/*         Int. J. Control, Vol. 36, pp. 1145-1193, 1984. */

/*     [2] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on an accuracy enhancing square-root */
/*     technique. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, April 1998. */
/*     Based on the RASP routine OHNAP1. */

/*     REVISIONS */

/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     April 24, 2000, A. Varga, DLR Oberpfaffenhofen. */
/*     April  8, 2001, A. Varga, DLR Oberpfaffenhofen. */
/*     March 26, 2005, V. Sima, Research Institute for Informatics. */

/*     KEYWORDS */

/*     Balancing, Hankel-norm approximation, model reduction, */
/*     multivariable system, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars */
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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --hsv;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    *iwarn = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Check the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldc < max(1,*p)) {
	*info = -12;
    } else if (*ldd < max(1,*p)) {
	*info = -14;
    } else if (*tol2 > 0. && *tol2 > *tol1) {
	*info = -17;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = max(*n,*m);
/* Computing MAX */
	i__4 = *m * 3 + 1, i__5 = min(*n,*m) + *p;
	i__1 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) / 2, i__2 =
		 *n * (*m + *p + 2) + (*m << 1) * *p + min(*n,*m) + max(i__4,
		i__5);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -20;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09CX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nr = 0;
	iwork[1] = 0;
	dwork[1] = 1.;
	return 0;
    }

    rtol = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
    srrtol = sqrt(rtol);

/*     Allocate working storage. */

    kt = 1;
    kti = kt + *n * *n;
    kw = kti + *n * *n;

/*     Compute a minimal order balanced realization of the given system. */
/*     Workspace: need   N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2; */
/*                prefer larger. */

    i__1 = *ldwork - kw + 1;
    ab09ax_(dico, "Balanced", "Automatic", n, m, p, &nminr, &a[a_offset], lda,
	     &b[b_offset], ldb, &c__[c_offset], ldc, &hsv[1], &dwork[kt], n, &
	    dwork[kti], n, tol2, &iwork[1], &dwork[kw], &i__1, iwarn, info, (
	    ftnlen)1, (ftnlen)8, (ftnlen)9);

    if (*info != 0) {
	return 0;
    }
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Compute the order of reduced system. */

    atol = rtol * hsv[1];
    if (fixord) {
	if (*nr > 0) {
	    if (*nr > nminr) {
		*nr = nminr;
		*iwarn = 1;
	    }
	}
    } else {
	atol = max(*tol1,atol);
	*nr = 0;
	i__1 = nminr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (hsv[i__] <= atol) {
		goto L20;
	    }
	    ++(*nr);
/* L10: */
	}
L20:
	;
    }

    if (*nr == nminr) {
	iwork[1] = nminr;
	dwork[1] = (doublereal) wrkopt;
	kw = *n * (*n + 2) + 1;

/*        Reduce Ar to a real Schur form. */

	i__1 = *ldwork - kw + 1;
	tb01wd_(&nminr, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &dwork[(*n << 1) + 1], n, &dwork[1], &dwork[*
		n + 1], &dwork[kw], &i__1, &ierr);
	if (ierr != 0) {
	    *info = 3;
	    return 0;
	}
	return 0;
    }
    skp = hsv[*nr + 1];

/*     If necessary, reduce the order such that HSV(NR) > HSV(NR+1). */

L30:
    if (*nr > 0) {
	if ((d__1 = hsv[*nr] - skp, abs(d__1)) <= srrtol * skp) {
	    --(*nr);
	    goto L30;
	}
    }

/*     Determine KR, the multiplicity of HSV(NR+1). */

    kr = 1;
    i__1 = nminr;
    for (i__ = *nr + 2; i__ <= i__1; ++i__) {
	if ((d__1 = hsv[i__] - skp, abs(d__1)) > srrtol * skp) {
	    goto L50;
	}
	++kr;
/* L40: */
    }
L50:

/*     For discrete-time case, apply the discrete-to-continuous bilinear */
/*     transformation. */

    if (discr) {

/*        Workspace: need   N; */
/*                   prefer larger. */

	ab04md_("Discrete", &nminr, m, p, &c_b16, &c_b16, &a[a_offset], lda, &
		b[b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &
		iwork[1], &dwork[1], ldwork, info, (ftnlen)8);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);
    }

/*     Define leading dimensions and offsets for temporary data. */

    nu = nminr - *nr - kr;
    na = *nr + nu;
    ldb1 = na;
    ldc1 = *p;
    ldb2 = kr;
    ldc2t = max(kr,*m);
    nr1 = *nr + 1;
/* Computing MIN */
    i__1 = nminr, i__2 = nr1 + kr;
    nkr1 = min(i__1,i__2);

    khsvp = 1;
    khsvp2 = khsvp + na;
    ku = khsvp2 + na;
    kb1 = ku + *p * *m;
    kb2 = kb1 + ldb1 * *m;
    kc1 = kb2 + ldb2 * *m;
    kc2t = kc1 + ldc1 * na;
    kw = kc2t + ldc2t * *p;

/*     Save B2 and C2'. */

    dlacpy_("Full", &kr, m, &b[nr1 + b_dim1], ldb, &dwork[kb2], &ldb2, (
	    ftnlen)4);
    ma02ad_("Full", p, &kr, &c__[nr1 * c_dim1 + 1], ldc, &dwork[kc2t], &ldc2t,
	     (ftnlen)4);
    if (*nr > 0) {

/*        Permute the elements of HSV and of matrices A, B, C. */

	dcopy_(nr, &hsv[1], &c__1, &dwork[khsvp], &c__1);
	dcopy_(&nu, &hsv[nkr1], &c__1, &dwork[khsvp + *nr], &c__1);
	dlacpy_("Full", &nminr, &nu, &a[nkr1 * a_dim1 + 1], lda, &a[nr1 * 
		a_dim1 + 1], lda, (ftnlen)4);
	dlacpy_("Full", &nu, &na, &a[nkr1 + a_dim1], lda, &a[nr1 + a_dim1], 
		lda, (ftnlen)4);
	dlacpy_("Full", &nu, m, &b[nkr1 + b_dim1], ldb, &b[nr1 + b_dim1], ldb,
		 (ftnlen)4);
	dlacpy_("Full", p, &nu, &c__[nkr1 * c_dim1 + 1], ldc, &c__[nr1 * 
		c_dim1 + 1], ldc, (ftnlen)4);

/*        Save B1 and C1. */

	dlacpy_("Full", &na, m, &b[b_offset], ldb, &dwork[kb1], &ldb1, (
		ftnlen)4);
	dlacpy_("Full", p, &na, &c__[c_offset], ldc, &dwork[kc1], &ldc1, (
		ftnlen)4);
    }

/*     Compute U = C2*pinv(B2'). */
/*     Workspace: need   N*(M+P+2) + 2*M*P + */
/*                       max(min(KR,M)+3*M+1,2*min(KR,M)+P); */
/*                prefer N*(M+P+2) + 2*M*P + */
/*                       max(min(KR,M)+2*M+(M+1)*NB,2*min(KR,M)+P*NB), */
/*                where  NB  is the maximum of the block sizes for */
/*                DGEQP3, DTZRZF, DTZRQF, DORMQR, and DORMRZ. */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	iwork[j] = 0;
/* L55: */
    }
    i__1 = *ldwork - kw + 1;
    dgelsy_(&kr, m, p, &dwork[kb2], &ldb2, &dwork[kc2t], &ldc2t, &iwork[1], &
	    rtol, &irank, &dwork[kw], &i__1, &ierr);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__1,i__2);
    ma02ad_("Full", m, p, &dwork[kc2t], &ldc2t, &dwork[ku], p, (ftnlen)4);

/*     Compute D <- D + HSV(NR+1)*U. */

    i__ = ku;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	daxpy_(p, &skp, &dwork[i__], &c__1, &d__[j * d_dim1 + 1], &c__1);
	i__ += *p;
/* L60: */
    }

    if (*nr > 0) {
	skp2 = skp * skp;

/*        Compute G = inv(S1*S1-skp*skp*I), where S1 is the diagonal */
/*        matrix of relevant singular values (of order NMINR - KR). */

	i1 = khsvp2;
	i__1 = khsvp + na - 1;
	for (i__ = khsvp; i__ <= i__1; ++i__) {
	    dwork[i1] = 1. / (dwork[i__] * dwork[i__] - skp2);
	    ++i1;
/* L70: */
	}

/*        Compute C <- C1*S1-skp*U*B1'. */

	mb01sd_("Column", p, &na, &c__[c_offset], ldc, &dwork[1], &dwork[
		khsvp], (ftnlen)6);
	d__1 = -skp;
	dgemm_("NoTranspose", "Transpose", p, &na, m, &d__1, &dwork[ku], p, &
		dwork[kb1], &ldb1, &c_b16, &c__[c_offset], ldc, (ftnlen)11, (
		ftnlen)9);

/*        Compute B <- G*(S1*B1-skp*C1'*U). */

	mb01sd_("Row", &na, m, &b[b_offset], ldb, &dwork[khsvp], &dwork[1], (
		ftnlen)3);
	d__1 = -skp;
	dgemm_("Transpose", "NoTranspose", &na, m, p, &d__1, &dwork[kc1], &
		ldc1, &dwork[ku], p, &c_b16, &b[b_offset], ldb, (ftnlen)9, (
		ftnlen)11);
	mb01sd_("Row", &na, m, &b[b_offset], ldb, &dwork[khsvp2], &dwork[1], (
		ftnlen)3);

/*        Compute A <- -A1' - B*B1'. */

	i__1 = na;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    dswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[j + a_dim1], lda);
/* L80: */
	}
	dgemm_("NoTranspose", "Transpose", &na, &na, m, &c_b49, &b[b_offset], 
		ldb, &dwork[kb1], &ldb1, &c_b49, &a[a_offset], lda, (ftnlen)
		11, (ftnlen)9);

/*        Extract stable part. */
/*        Workspace:  need   N*N+5*N; */
/*                    prefer larger. */

	kw1 = na * na + 1;
	kw2 = kw1 + na;
	kw = kw2 + na;
	i__1 = *ldwork - kw + 1;
	tb01kd_("Continuous", "Stability", "General", &na, m, p, &c_b54, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, &ndim,
		 &dwork[1], &na, &dwork[kw1], &dwork[kw2], &dwork[kw], &i__1, 
		&ierr, (ftnlen)10, (ftnlen)9, (ftnlen)7);
	if (ierr != 0) {
	    *info = 3;
	    return 0;
	}

	if (ndim != *nr) {
	    *info = 4;
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);

/*        For discrete-time case, apply the continuous-to-discrete */
/*        bilinear transformation. */

	if (discr) {
	    ab04md_("Continuous", nr, m, p, &c_b16, &c_b16, &a[a_offset], lda,
		     &b[b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], 
		    ldd, &iwork[1], &dwork[1], ldwork, info, (ftnlen)10);
	}
    }
    iwork[1] = nminr;
    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of AB09CX *** */
} /* ab09cx_ */

