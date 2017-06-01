/* AG08BD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static doublereal c_b25 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ag08bd_(char *equil, integer *l, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *e, integer *lde, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *nfz, integer *nrank, integer *
	niz, integer *dinfz, integer *nkror, integer *ninfe, integer *nkrol, 
	integer *infz, integer *kronr, integer *infe, integer *kronl, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, i0, i1, n2, nb, ii, mm, nn, pp, mu, nu, ipd;
    static doublereal dum[1];
    static integer ldw, itau, numu, kabcd;
    extern /* Subroutine */ int ma02bd_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), ma02cd_(integer *, integer *, 
	    integer *, doublereal *, integer *), tg01ad_(char *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, ftnlen), tg01fd_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), ag08by_(logical *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer labcd2;
    static doublereal toler;
    static integer jwork, ldabcd;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer nsinfe;
    static logical lequil;
    static doublereal svlmax;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
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

/*     To extract from the system pencil */

/*                       ( A-lambda*E B ) */
/*           S(lambda) = (              ) */
/*                       (      C     D ) */

/*     a regular pencil Af-lambda*Ef which has the finite Smith zeros of */
/*     S(lambda) as generalized eigenvalues. The routine also computes */
/*     the orders of the infinite Smith zeros and determines the singular */
/*     and infinite Kronecker structure of system pencil, i.e., the right */
/*     and left Kronecker indices, and the multiplicities of infinite */
/*     eigenvalues. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the system */
/*             matrix as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Af of the reduced pencil. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Ef of the reduced pencil. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0; */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             direct transmission matrix D of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NFZ     (output) INTEGER */
/*             The number of finite zeros. */

/*     NRANK   (output) INTEGER */
/*             The normal rank of the system pencil. */

/*     NIZ     (output) INTEGER */
/*             The number of infinite zeros. */

/*     DINFZ   (output) INTEGER */
/*             The maximal multiplicity of infinite Smith zeros. */

/*     NKROR   (output) INTEGER */
/*             The number of right Kronecker indices. */

/*     NINFE   (output) INTEGER */
/*             The number of elementary infinite blocks. */

/*     NKROL   (output) INTEGER */
/*             The number of left Kronecker indices. */

/*     INFZ    (output) INTEGER array, dimension (N+1) */
/*             The leading DINFZ elements of INFZ contain information */
/*             on the infinite elementary divisors as follows: */
/*             the system has INFZ(i) infinite elementary divisors of */
/*             degree i in the Smith form, where i = 1,2,...,DINFZ. */

/*     KRONR   (output) INTEGER array, dimension (N+M+1) */
/*             The leading NKROR elements of this array contain the */
/*             right Kronecker (column) indices. */

/*     INFE    (output) INTEGER array, dimension (1+MIN(L+P,N+M)) */
/*             The leading NINFE elements of INFE contain the */
/*             multiplicities of infinite eigenvalues. */

/*     KRONL   (output) INTEGER array, dimension (L+P+1) */
/*             The leading NKROL elements of this array contain the */
/*             left Kronecker (row) indices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             If the user sets TOL <= 0, then default tolerances are */
/*             used instead, as follows: TOLDEF = L*N*EPS in TG01FD */
/*             (to determine the rank of E) and TOLDEF = (L+P)*(N+M)*EPS */
/*             in the rest, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N+max(1,M) */
/*             On output, IWORK(1) contains the normal rank of the */
/*             transfer function matrix. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 4*(L+N), LDW ), if EQUIL = 'S', */
/*             LDWORK >= LDW,                 if EQUIL = 'N', where */
/*             LDW = max(L+P,M+N)*(M+N) + max(1,5*max(L+P,M+N)). */
/*             For optimum performance LDWORK should be larger. */

/*             If LDWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             DWORK array, returns this value as the first entry of */
/*             the DWORK array, and no error message related to LDWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine extracts from the system matrix of a descriptor */
/*     system (A-lambda*E,B,C,D) a regular pencil Af-lambda*Ef which */
/*     has the finite zeros of the system as generalized eigenvalues. */
/*     The procedure has the following main computational steps: */

/*        (a) construct the (L+P)-by-(N+M) system pencil */

/*             S(lambda) = ( B  A )-lambda*( 0  E ); */
/*                         ( D  C )        ( 0  0 ) */

/*        (b) reduce S(lambda) to S1(lambda) with the same finite */
/*            zeros and right Kronecker structure but with E */
/*            upper triangular and nonsingular; */

/*        (c) reduce S1(lambda) to S2(lambda) with the same finite */
/*            zeros and right Kronecker structure but with D of */
/*            full row rank; */

/*        (d) reduce S2(lambda) to S3(lambda) with the same finite zeros */
/*            and with D square invertible; */

/*        (e) perform a unitary transformation on the columns of */

/*            S3(lambda) = (A-lambda*E   B) in order to reduce it to */
/*                         (     C       D) */

/*            (Af-lambda*Ef   X), with Y and Ef square invertible; */
/*            (     0         Y) */

/*        (f) compute the right and left Kronecker indices of the system */
/*            matrix, which together with the multiplicities of the */
/*            finite and infinite eigenvalues constitute the */
/*            complete set of structural invariants under strict */
/*            equivalence transformations of a linear system. */

/*     REFERENCES */

/*     [1] P. Misra, P. Van Dooren and A. Varga. */
/*         Computation of structural invariants of generalized */
/*         state-space systems. */
/*         Automatica, 30, pp. 1921-1936, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable (see [1]). */

/*     FURTHER COMMENTS */

/*     In order to compute the finite Smith zeros of the system */
/*     explicitly, a call to this routine may be followed by a */
/*     call to the LAPACK Library routines DGEGV or DGGEV. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     May 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Sep. 1999, */
/*     Jan. 2009, Mar. 2009, Apr. 2009. */
/*     A. Varga, DLR Oberpfaffenhofen, Nov. 1999, Feb. 2002, Mar. 2002. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, orthogonal transformation, structural invariant. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --infz;
    --kronr;
    --infe;
    --kronl;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
/* Computing MAX */
    i__1 = *l + *p, i__2 = *n + *m;
    ldabcd = max(i__1,i__2);
    labcd2 = ldabcd * (*n + *m);
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
    lquery = *ldwork == -1;

/*     Test the input scalar arguments. */

    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*l < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*lda < max(1,*l)) {
	*info = -7;
    } else if (*lde < max(1,*l)) {
	*info = -9;
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
	*info = -11;
    } else if (*ldc < max(1,*p)) {
	*info = -13;
    } else if (*ldd < max(1,*p)) {
	*info = -15;
    } else if (*tol >= 1.) {
	*info = -27;
    } else {
/* Computing MIN */
	i__1 = *l + *p, i__2 = *m + *n;
	i0 = min(i__1,i__2);
	i1 = min(*l,*n);
	ii = min(*m,*p);
/* Computing MAX */
	i__1 = 1, i__2 = ldabcd * 5;
	ldw = labcd2 + max(i__1,i__2);
	if (lequil) {
/* Computing MAX */
	    i__1 = *l + *n << 2;
	    ldw = max(i__1,ldw);
	}
	if (lquery) {
	    tg01fd_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset]
		    , lde, &b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1,
		     dum, &c__1, &nn, &n2, tol, &iwork[1], &dwork[1], &c_n1, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
	    i__1 = ldw, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    svlmax = 0.;
	    i__1 = *m + *n;
	    i__2 = *p + *l;
	    i__3 = ldabcd + i1;
	    ag08by_(&c_true, &i1, &i__1, &i__2, &svlmax, &dwork[1], &i__3, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = labcd2 + (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    i__1 = *m + *n;
	    i__2 = ldabcd + i1;
	    ag08by_(&c_false, &i1, &ii, &i__1, &svlmax, &dwork[1], &i__2, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = labcd2 + (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    i__1 = i1 + ii;
	    nb = ilaenv_(&c__1, "ZGERQF", " ", &ii, &i__1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = labcd2 + ii + ii * nb;
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
	    i__3 = i1 + ii;
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RT", &i1, &i__3, &ii, 
		    &c_n1, (ftnlen)6, (ftnlen)2);
	    nb = min(i__1,i__2);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = labcd2 + ii + i1 * nb;
	    wrkopt = max(i__1,i__2);
	} else if (*ldwork < ldw) {
	    *info = -30;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AG08BD", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

    *niz = 0;
    *nkrol = 0;
    *nkror = 0;

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*l,*n), i__1 = max(i__1,*m);
    if (max(i__1,*p) == 0) {
	*nfz = 0;
	*dinfz = 0;
	*ninfe = 0;
	*nrank = 0;
	iwork[1] = 0;
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

    wrkopt = 1;
    kabcd = 1;
    jwork = kabcd + labcd2;

/*     If required, balance the system pencil. */
/*     Workspace: need   4*(L+N). */

    if (lequil) {
	tg01ad_("A", l, n, m, p, &c_b25, &a[a_offset], lda, &e[e_offset], lde,
		 &b[b_offset], ldb, &c__[c_offset], ldc, &dwork[1], &dwork[*l 
		+ 1], &dwork[*l + *n + 1], info, (ftnlen)1);
	wrkopt = *l + *n << 2;
    }
    svlmax = dlange_("Frobenius", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)
	    9);

/*     Reduce the system matrix to QR form, */

/*          ( A11-lambda*E11 A12 B1 ) */
/*          (     A21        A22 B2 ) , */
/*          (     C1         C2  D  ) */

/*     with E11 invertible and upper triangular. */
/*     Real workspace: need   max( 1, N+P, min(L,N)+max(3*N-1,M,L) ); */
/*                     prefer larger. */
/*     Integer workspace: N. */

    tg01fd_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset], lde, &
	    b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1, dum, &c__1, &
	    nn, &n2, tol, &iwork[1], &dwork[1], ldwork, info, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[1];
    wrkopt = max(i__1,i__2);

/*     Construct the system pencil */

/*                          MM         NN */
/*                      ( B1 A12 A11-lambda*E11 ) NN */
/*        S1(lambda) =  ( B2 A22      A21       ) L-NN */
/*                      ( D  C2       C1        ) P */

/*     of dimension (L+P)-by-(M+N). */
/*     Workspace: need  LABCD2 = max( L+P, N+M )*( N+M ). */

    n2 = *n - nn;
    mm = *m + n2;
    pp = *p + (*l - nn);
    dlacpy_("Full", l, m, &b[b_offset], ldb, &dwork[kabcd], &ldabcd, (ftnlen)
	    4);
    dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[kabcd + *l], &ldabcd, (
	    ftnlen)4);
    dlacpy_("Full", l, &n2, &a[(nn + 1) * a_dim1 + 1], lda, &dwork[kabcd + 
	    ldabcd * *m], &ldabcd, (ftnlen)4);
    dlacpy_("Full", p, &n2, &c__[(nn + 1) * c_dim1 + 1], ldc, &dwork[kabcd + 
	    ldabcd * *m + *l], &ldabcd, (ftnlen)4);
    dlacpy_("Full", l, &nn, &a[a_offset], lda, &dwork[kabcd + ldabcd * mm], &
	    ldabcd, (ftnlen)4);
    dlacpy_("Full", p, &nn, &c__[c_offset], ldc, &dwork[kabcd + ldabcd * mm + 
	    *l], &ldabcd, (ftnlen)4);

/*     If required, set tolerance. */

    toler = *tol;
    if (toler <= 0.) {
	toler = (doublereal) ((*l + *p) * (*m + *n)) * dlamch_("Precision", (
		ftnlen)9);
    }
/* Computing MAX */
    i__1 = nn + pp;
    i__2 = nn + mm;
    d__1 = svlmax, d__2 = dlange_("Frobenius", &i__1, &i__2, &dwork[kabcd], &
	    ldabcd, &dwork[jwork], (ftnlen)9);
    svlmax = max(d__1,d__2);

/*     Extract the reduced pencil S2(lambda) */

/*             ( Bc  Ac-lambda*Ec ) */
/*             ( Dc      Cc       ) */

/*     having the same finite Smith zeros as the system pencil */
/*     S(lambda) but with Dc, a MU-by-MM full row rank */
/*     left upper trapezoidal matrix, and Ec, an NU-by-NU */
/*     upper triangular nonsingular matrix. */

/*     Real workspace: need   max( min(P+L,M+N)+max(min(L,N),3*(M+N)-1), */
/*                                  5*(P+L), 1 ) + LABCD2; */
/*                     prefer larger. */
/*     Integer workspace: MM, MM <= M+N; PP <= P+L. */

    i__1 = *ldwork - jwork + 1;
    ag08by_(&c_true, &nn, &mm, &pp, &svlmax, &dwork[kabcd], &ldabcd, &e[
	    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &kronl[1], 
	    &toler, &iwork[1], &dwork[jwork], &i__1, info);

/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

/*     Set the number of simple (nondynamic) infinite eigenvalues */
/*     and the normal rank of the system pencil. */

    nsinfe = mu;
    *nrank = nn + mu;

/*     Pertranspose the system. */

/* Computing MAX */
    i__2 = 0, i__3 = nu - 1;
    i__1 = max(i__2,i__3);
/* Computing MAX */
    i__5 = 0, i__6 = nu - 1;
    i__4 = max(i__5,i__6);
    tb01xd_("D", &nu, &mm, &mm, &i__1, &i__4, &dwork[kabcd + ldabcd * mm], &
	    ldabcd, &dwork[kabcd], &ldabcd, &dwork[kabcd + ldabcd * mm + nu], 
	    &ldabcd, &dwork[kabcd + nu], &ldabcd, info, (ftnlen)1);
    i__1 = nu + mm;
    ma02bd_("Right", &i__1, &mm, &dwork[kabcd], &ldabcd, (ftnlen)5);
    i__1 = nu + mm;
    ma02bd_("Left", &mm, &i__1, &dwork[kabcd + nu], &ldabcd, (ftnlen)4);
/* Computing MAX */
    i__2 = 0, i__3 = nu - 1;
    i__1 = max(i__2,i__3);
    ma02cd_(&nu, &c__0, &i__1, &e[e_offset], lde);

    if (mu != mm) {
	nn = nu;
	pp = mm;
	mm = mu;
	kabcd += (pp - mm) * ldabcd;

/*        Extract the reduced pencil S3(lambda), */

/*             ( Br  Ar-lambda*Er ) , */
/*             ( Dr      Cr       ) */

/*        having the same finite Smith zeros as the pencil S(lambda), */
/*        but with Dr, an MU-by-MU invertible upper triangular matrix, */
/*        and Er, an NU-by-NU upper triangular nonsingular matrix. */

/*        Workspace: need   max( 1, 5*(M+N) ) + LABCD2. */
/*                   prefer larger. */
/*        No integer workspace necessary. */

	i__1 = *ldwork - jwork + 1;
	ag08by_(&c_false, &nn, &mm, &pp, &svlmax, &dwork[kabcd], &ldabcd, &e[
		e_offset], lde, &nu, &mu, &i0, &i1, nkror, &iwork[1], &kronr[
		1], &toler, &iwork[1], &dwork[jwork], &i__1, info);

/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
    }

    if (nu != 0) {

/*        Perform a unitary transformation on the columns of */
/*                     ( Br Ar-lambda*Er ) */
/*                     ( Dr     Cr       ) */
/*        in order to reduce it to */
/*                     ( *  Af-lambda*Ef ) */
/*                     ( Y       0       ) */
/*        with Y and Ef square invertible. */

/*        Compute Af by reducing  ( Br Ar ) to  ( *  Af ) . */
/*                                ( Dr Cr )     ( Y   0 ) */

	numu = nu + mu;
	ipd = kabcd + nu;
	itau = jwork;
	jwork = itau + mu;

/*        Workspace: need   LABCD2 + 2*min(M,P); */
/*                   prefer LABCD2 + min(M,P) + min(M,P)*NB. */

	i__1 = *ldwork - jwork + 1;
	dtzrzf_(&mu, &numu, &dwork[ipd], &ldabcd, &dwork[itau], &dwork[jwork],
		 &i__1, info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Workspace: need   LABCD2 + min(M,P) + min(L,N); */
/*                   prefer LABCD2 + min(M,P) + min(L,N)*NB. */

	i__1 = *ldwork - jwork + 1;
	dormrz_("Right", "Transpose", &nu, &numu, &mu, &nu, &dwork[ipd], &
		ldabcd, &dwork[itau], &dwork[kabcd], &ldabcd, &dwork[jwork], &
		i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Save Af. */

	dlacpy_("Full", &nu, &nu, &dwork[kabcd + ldabcd * mu], &ldabcd, &a[
		a_offset], lda, (ftnlen)4);

/*        Compute Ef by applying the saved transformations from previous */
/*        reduction to ( 0  Er ) . */

	dlaset_("Full", &nu, &mu, &c_b25, &c_b25, &dwork[kabcd], &ldabcd, (
		ftnlen)4);
	dlacpy_("Full", &nu, &nu, &e[e_offset], lde, &dwork[kabcd + ldabcd * 
		mu], &ldabcd, (ftnlen)4);

	i__1 = *ldwork - jwork + 1;
	dormrz_("Right", "Transpose", &nu, &numu, &mu, &nu, &dwork[ipd], &
		ldabcd, &dwork[itau], &dwork[kabcd], &ldabcd, &dwork[jwork], &
		i__1, info, (ftnlen)5, (ftnlen)9);

/*        Save Ef. */

	dlacpy_("Full", &nu, &nu, &dwork[kabcd + ldabcd * mu], &ldabcd, &e[
		e_offset], lde, (ftnlen)4);
    }

    *nfz = nu;

/*     Set right Kronecker indices (column indices). */

    i__1 = *nkror;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = kronr[i__];
/* L10: */
    }

    j = 0;
    i__1 = *nkror;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = j + iwork[i__];
	for (ii = j + 1; ii <= i__2; ++ii) {
	    kronr[ii] = i__ - 1;
/* L20: */
	}
	j += iwork[i__];
/* L30: */
    }

    *nkror = j;

/*     Set left Kronecker indices (row indices). */

    i__1 = *nkrol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = kronl[i__];
/* L40: */
    }

    j = 0;
    i__1 = *nkrol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = j + iwork[i__];
	for (ii = j + 1; ii <= i__2; ++ii) {
	    kronl[ii] = i__ - 1;
/* L50: */
	}
	j += iwork[i__];
/* L60: */
    }

    *nkrol = j;

/*     Determine the number of simple infinite blocks */
/*     as the difference between the number of infinite blocks */
/*     of order greater than one and the order of Dr. */

    *ninfe = 0;
    i__1 = *dinfz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*ninfe += infz[i__];
/* L70: */
    }
    *ninfe = nsinfe - *ninfe;
    i__1 = *ninfe;
    for (i__ = 1; i__ <= i__1; ++i__) {
	infe[i__] = 1;
/* L80: */
    }

/*     Set the structure of infinite eigenvalues. */

    i__1 = *dinfz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ninfe + infz[i__];
	for (ii = *ninfe + 1; ii <= i__2; ++ii) {
	    infe[ii] = i__ + 1;
/* L90: */
	}
	*ninfe += infz[i__];
/* L100: */
    }

    iwork[1] = nsinfe;
    dwork[1] = (doublereal) wrkopt;
    return 0;
/* *** Last line of AG08BD *** */
} /* ag08bd_ */

