/* AB09HX.f -- translated by f2c (version 20100827).
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

static doublereal c_b14 = 1.;
static integer c__1 = 1;
static doublereal c_b49 = 0.;

/* Subroutine */ int ab09hx_(char *dico, char *job, char *ordsel, integer *n, 
	integer *m, integer *p, integer *nr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *hsv, doublereal *t, 
	integer *ldt, doublereal *ti, integer *ldti, doublereal *tol1, 
	doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork, 
	logical *bwork, integer *iwarn, integer *info, ftnlen dico_len, 
	ftnlen job_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, t_dim1, t_offset, ti_dim1, ti_offset, i__1, i__2, i__3, 
	    i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, ij, ku, kv, kw, lw, ns, nr1;
    static logical bal, bta, spa;
    static integer ldw;
    static doublereal atol;
    static integer ierr, ktau;
    static doublereal temp;
    extern /* Subroutine */ int ab09dd_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), ma02ad_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), ab04md_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dscal_(integer *, doublereal *, doublereal *, integer *), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), mb03ud_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen),
	     ab09hy_(integer *, integer *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    logical *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    static doublereal rcond;
    static integer nminr;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrmv_(
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef, ricond;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgetrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static logical fixord;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     stable state-space representation (A,B,C,D) by using the */
/*     stochastic balancing approach in conjunction with the square-root */
/*     or the balancing-free square-root Balance & Truncate (B&T) or */
/*     Singular Perturbation Approximation (SPA) model reduction methods. */
/*     The state dynamics matrix A of the original system is an upper */
/*     quasi-triangular matrix in real Schur canonical form and D must be */
/*     full row rank. */

/*     For the B&T approach, the matrices of the reduced order system */
/*     are computed using the truncation formulas: */

/*          Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T .     (1) */

/*     For the SPA approach, the matrices of a minimal realization */
/*     (Am,Bm,Cm) are computed using the truncation formulas: */

/*          Am = TI * A * T ,  Bm = TI * B ,  Cm = C * T .     (2) */

/*     Am, Bm, Cm and D serve further for computing the SPA of the given */
/*     system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'F':  use the balancing-free square-root */
/*                     Balance & Truncate method; */
/*             = 'S':  use the square-root Singular Perturbation */
/*                     Approximation method; */
/*             = 'P':  use the balancing-free square-root */
/*                     Singular Perturbation Approximation method. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  M >= P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of */
/*             the resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR */
/*             is the desired order on entry and NMIN is the order of a */
/*             minimal realization of the given system; NMIN is */
/*             determined as the number of Hankel singular values greater */
/*             than N*EPS, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,N*EPS). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A in a real Schur */
/*             canonical form. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system. */

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
/*             If INFO = 0, it contains the Hankel singular values, */
/*             ordered decreasingly, of the phase system. All singular */
/*             values are less than or equal to 1. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             If INFO = 0 and NR > 0, the leading N-by-NR part of this */
/*             array contains the right truncation matrix T in (1), for */
/*             the B&T approach, or in (2), for the SPA approach. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N) */
/*             If INFO = 0 and NR > 0, the leading NR-by-N part of this */
/*             array contains the left truncation matrix TI in (1), for */
/*             the B&T approach, or in (2), for the SPA approach. */

/*     LDTI    INTEGER */
/*             The leading dimension of array TI.  LDTI >= MAX(1,N). */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value lies in the */
/*             interval [0.00001,0.001]. */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = N*EPS, where EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the phase system (see METHOD) corresponding */
/*             to the given system. */
/*             The recommended value is TOL2 = N*EPS. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(1,2*N) */
/*             On exit with INFO = 0, IWORK(1) contains the order of the */
/*             minimal realization of the system. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) contains RCOND, the reciprocal */
/*             condition number of the U11 matrix from the expression */
/*             used to compute the solution X = U21*inv(U11) of the */
/*             Riccati equation for spectral factorization. */
/*             A small value RCOND indicates possible ill-conditioning */
/*             of the respective Riccati equation. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 2, N*(MAX(N,M,P)+5), */
/*                            2*N*P+MAX(P*(M+2),10*N*(N+1) ) ). */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension 2*N */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than the order of a minimal realization of the */
/*                   given system. In this case, the resulting NR is */
/*                   set automatically to a value corresponding to the */
/*                   order of a minimal realization of the system. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'), or it is not in */
/*                   a real Schur form; */
/*             = 2:  the reduction of Hamiltonian matrix to real */
/*                   Schur form failed; */
/*             = 3:  the reordering of the real Schur form of the */
/*                   Hamiltonian matrix failed; */
/*             = 4:  the Hamiltonian matrix has less than N stable */
/*                   eigenvalues; */
/*             = 5:  the coefficient matrix U11 in the linear system */
/*                   X*U11 = U21, used to determine X, is singular to */
/*                   working precision; */
/*             = 6:  the feedthrough matrix D has not a full row rank P; */
/*             = 7:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                             (3) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09HX determines for */
/*     the given system (3), the matrices of a reduced NR-rder system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                         (4) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     If JOB = 'B', the square-root stochastic Balance & Truncate */
/*     method of [1] is used and the resulting model is balanced. */

/*     If JOB = 'F', the balancing-free square-root version of the */
/*     stochastic Balance & Truncate method [1] is used. */

/*     If JOB = 'S', the stochastic balancing method, in conjunction */
/*     with the square-root version of the Singular Perturbation */
/*     Approximation method [2,3] is used. */

/*     If JOB = 'P', the stochastic balancing method, in conjunction */
/*     with the balancing-free square-root version of the Singular */
/*     Perturbation Approximation method [2,3] is used. */

/*     By setting TOL1 = TOL2, the routine can be also used to compute */
/*     Balance & Truncate approximations. */

/*     REFERENCES */

/*     [1] Varga A. and Fasol K.H. */
/*         A new square-root balancing-free stochastic truncation */
/*         model reduction algorithm. */
/*         Proc. of 12th IFAC World Congress, Sydney, 1993. */

/*     [2] Liu Y. and Anderson B.D.O. */
/*         Singular Perturbation Approximation of balanced systems. */
/*         Int. J. Control, Vol. 50, pp. 1379-1405, 1989. */

/*     [3] Varga A. */
/*         Balancing-free square-root algorithm for computing singular */
/*         perturbation approximations. */
/*         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991, */
/*         Vol. 2, pp. 1062-1065. */

/*     NUMERICAL ASPECTS */

/*     The implemented method relies on accuracy enhancing square-root */
/*     or balancing-free square-root methods. The effectiveness of the */
/*     accuracy enhancing technique depends on the accuracy of the */
/*     solution of a Riccati equation. Ill-conditioned Riccati solution */
/*     typically results when D is nearly rank deficient. */
/*                                      3 */
/*     The algorithm requires about 100N  floating point operations. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000. */
/*     D. Sima, University of Bucharest, May 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2000. */
/*     Partly based on the RASP routine SRBFS1, by A. Varga, 1992. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001. */

/*     KEYWORDS */

/*     Balance and truncate, minimal state-space representation, */
/*     model reduction, multivariable system, */
/*     singular perturbation approximation, state-space model, */
/*     stochastic balancing. */

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
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --hsv;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    ti_dim1 = *ldti;
    ti_offset = 1 + ti_dim1;
    ti -= ti_offset;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    *info = 0;
    *iwarn = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    bta = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "F", (ftnlen)
	    1, (ftnlen)1);
    spa = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "P", (ftnlen)
	    1, (ftnlen)1);
    bal = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "S", (ftnlen)
	    1, (ftnlen)1);
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
/* Computing MAX */
/* Computing MAX */
    i__3 = max(*n,*m);
/* Computing MAX */
    i__4 = *p * (*m + 2), i__5 = *n * 10 * (*n + 1);
    i__1 = 2, i__2 = *n * (max(i__3,*p) + 5), i__1 = max(i__1,i__2), i__2 = (*
	    n << 1) * *p + max(i__4,i__5);
    lw = max(i__1,i__2);

/*     Test the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (bta || spa)) {
	*info = -2;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*p < 0 || *p > *m) {
	*info = -6;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -7;
    } else if (*lda < max(1,*n)) {
	*info = -9;
    } else if (*ldb < max(1,*n)) {
	*info = -11;
    } else if (*ldc < max(1,*p)) {
	*info = -13;
    } else if (*ldd < max(1,*p)) {
	*info = -15;
    } else if (*ldt < max(1,*n)) {
	*info = -18;
    } else if (*ldti < max(1,*n)) {
	*info = -20;
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
	*info = -22;
    } else if (*ldwork < lw) {
	*info = -25;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09HX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nr = 0;
	iwork[1] = 0;
	dwork[1] = 2.;
	dwork[2] = 1.;
	return 0;
    }

/*     For discrete-time case, apply the discrete-to-continuous bilinear */
/*     transformation. */

    if (discr) {

/*        Real workspace:    need  N, prefer larger; */
/*        Integer workspace: need  N. */

	ab04md_("Discrete", n, m, p, &c_b14, &c_b14, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &
		iwork[1], &dwork[1], ldwork, &ierr, (ftnlen)8);
	if (ierr != 0) {
	    *info = 1;
	    return 0;
	}
/* Computing MAX */
	i__1 = *n, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);
    } else {
	wrkopt = 0;
    }

/*     Compute in TI and T the Cholesky factors Su and Ru of the */
/*     controllability and observability Grammians, respectively. */
/*     Real workspace:    need  MAX( 2, N*(MAX(N,M,P)+5), */
/*                                   2*N*P+MAX(P*(M+2),10*N*(N+1) ) ); */
/*                        prefer larger. */
/*     Integer workspace: need  2*N. */

    ab09hy_(n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], 
	    ldc, &d__[d_offset], ldd, &scalec, &scaleo, &ti[ti_offset], ldti, 
	    &t[t_offset], ldt, &iwork[1], &dwork[1], ldwork, &bwork[1], info);
    if (*info != 0) {
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[1];
    wrkopt = max(i__1,i__2);
    ricond = dwork[2];

/*     Save Su in V. */

    ku = 1;
    kv = ku + *n * *n;
    kw = kv + *n * *n;
    dlacpy_("Upper", n, n, &ti[ti_offset], ldti, &dwork[kv], n, (ftnlen)5);
/*                               | x x | */
/*     Compute Ru*Su in the form | 0 x | in TI. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dtrmv_("Upper", "NoTranspose", "NonUnit", &j, &t[t_offset], ldt, &ti[
		j * ti_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
/* L10: */
    }

/*     Compute the singular value decomposition Ru*Su = V*S*UT */
/*     of the upper triangular matrix Ru*Su, with UT in TI and V in U. */

/*     Workspace:  need   2*N*N + 5*N; */
/*                 prefer larger. */

    i__1 = *ldwork - kw + 1;
    mb03ud_("Vectors", "Vectors", n, &ti[ti_offset], ldti, &dwork[ku], n, &
	    hsv[1], &dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)7);
    if (ierr != 0) {
	*info = 7;
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__1,i__2);

/*     Scale the singular values. */

    d__1 = 1. / scalec / scaleo;
    dscal_(n, &d__1, &hsv[1], &c__1);

/*     Partition S, U and V conformally as: */

/*     S = diag(S1,S2,S3),  U = [U1,U2,U3] (U' in TI) and V = [V1,V2,V3] */
/*     (in U). */

/*     Compute the order NR of reduced system, as the order of S1. */

    toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
    atol = toldef;
    if (fixord) {
	if (*nr > 0) {
	    if (hsv[*nr] <= atol) {
		*nr = 0;
		*iwarn = 1;
		fixord = FALSE_;
	    }
	}
    } else {
	atol = max(*tol1,atol);
	*nr = 0;
    }
    if (! fixord) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (hsv[j] <= atol) {
		goto L30;
	    }
	    ++(*nr);
/* L20: */
	}
L30:
	;
    }

/*     Compute the order of minimal realization as the order of [S1 S2]. */

    nr1 = *nr + 1;
    nminr = *nr;
    if (*nr < *n) {
	if (spa) {
	    atol = max(*tol2,toldef);
	}
	i__1 = *n;
	for (j = nr1; j <= i__1; ++j) {
	    if (hsv[j] <= atol) {
		goto L50;
	    }
	    ++nminr;
/* L40: */
	}
L50:
	;
    }

/*     Finish if the order is zero. */

    if (*nr == 0) {
	if (spa) {
	    ab09dd_("Continuous", n, m, p, nr, &a[a_offset], lda, &b[b_offset]
		    , ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &rcond, &
		    iwork[1], &dwork[1], &ierr, (ftnlen)10);
	    iwork[1] = nminr;
	} else {
	    iwork[1] = 0;
	}
	dwork[1] = (doublereal) wrkopt;
	dwork[2] = ricond;
	return 0;
    }

/*     Compute NS, the order of S2. */
/*     Note: For BTA, NS is always zero, because NMINR = NR. */

    ns = nminr - *nr;

/*     Compute the truncation matrices. */

/*     Compute TI' = | TI1' TI2' | = Ru'*| V1 V2 | in U. */

    dtrmm_("Left", "Upper", "Transpose", "NonUnit", n, &nminr, &c_b14, &t[
	    t_offset], ldt, &dwork[ku], n, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute  T = | T1 T2 | = Su*| U1 U2 | . */

    ma02ad_("Full", &nminr, n, &ti[ti_offset], ldti, &t[t_offset], ldt, (
	    ftnlen)4);
    dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, &nminr, &c_b14, &
	    dwork[kv], n, &t[t_offset], ldt, (ftnlen)4, (ftnlen)5, (ftnlen)11,
	     (ftnlen)7);
    ktau = kv;

    if (bal) {
	ij = ku;

/*        Square-Root B&T/SPA method. */

/*        Compute the truncation matrices for balancing */
/*                    -1/2            -1/2 */
/*               T1*S1     and TI1'*S1    . */

	i__1 = *nr;
	for (j = 1; j <= i__1; ++j) {
	    temp = 1. / sqrt(hsv[j]);
	    dscal_(n, &temp, &t[j * t_dim1 + 1], &c__1);
	    dscal_(n, &temp, &dwork[ij], &c__1);
	    ij += *n;
/* L70: */
	}
    } else {

/*        Balancing-Free B&T/SPA method. */

/*        Compute orthogonal bases for the images of matrices T1 and */
/*        TI1'. */

/*        Workspace:  need   N*MAX(N,M,P) + 2*NR; */
/*                    prefer N*MAX(N,M,P) + NR*(NB+1) */
/*                           (NB determined by ILAENV for DGEQRF). */

	kw = ktau + *nr;
	ldw = *ldwork - kw + 1;
	dgeqrf_(n, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
	dorgqr_(n, nr, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, 
		&ierr);
	dgeqrf_(n, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
	dorgqr_(n, nr, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
    }
    if (ns > 0) {

/*        Compute orthogonal bases for the images of matrices T2 and */
/*        TI2'. */

/*        Workspace:  need   N*MAX(N,M,P) + 2*NS; */
/*                    prefer N*MAX(N,M,P) + NS*(NB+1) */
/*                           (NB determined by ILAENV for DGEQRF). */
	kw = ktau + ns;
	ldw = *ldwork - kw + 1;
	dgeqrf_(n, &ns, &t[nr1 * t_dim1 + 1], ldt, &dwork[ktau], &dwork[kw], &
		ldw, &ierr);
	dorgqr_(n, &ns, &ns, &t[nr1 * t_dim1 + 1], ldt, &dwork[ktau], &dwork[
		kw], &ldw, &ierr);
	dgeqrf_(n, &ns, &dwork[ku + *n * *nr], n, &dwork[ktau], &dwork[kw], &
		ldw, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
	dorgqr_(n, &ns, &ns, &dwork[ku + *n * *nr], n, &dwork[ktau], &dwork[
		kw], &ldw, &ierr);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
    }

/*     Transpose TI' in TI. */

    ma02ad_("Full", n, &nminr, &dwork[ku], n, &ti[ti_offset], ldti, (ftnlen)4)
	    ;

    if (! bal) {
/*                        -1 */
/*        Compute (TI1*T1)  *TI1 in TI. */

	dgemm_("NoTranspose", "NoTranspose", nr, nr, n, &c_b14, &ti[ti_offset]
		, ldti, &t[t_offset], ldt, &c_b49, &dwork[ku], n, (ftnlen)11, 
		(ftnlen)11);
	dgetrf_(nr, nr, &dwork[ku], n, &iwork[1], &ierr);
	dgetrs_("NoTranspose", nr, n, &dwork[ku], n, &iwork[1], &ti[ti_offset]
		, ldti, &ierr, (ftnlen)11);

	if (ns > 0) {
/*                           -1 */
/*           Compute (TI2*T2)  *TI2 in TI2. */

	    dgemm_("NoTranspose", "NoTranspose", &ns, &ns, n, &c_b14, &ti[nr1 
		    + ti_dim1], ldti, &t[nr1 * t_dim1 + 1], ldt, &c_b49, &
		    dwork[ku], n, (ftnlen)11, (ftnlen)11);
	    dgetrf_(&ns, &ns, &dwork[ku], n, &iwork[1], &ierr);
	    dgetrs_("NoTranspose", &ns, n, &dwork[ku], n, &iwork[1], &ti[nr1 
		    + ti_dim1], ldti, &ierr, (ftnlen)11);
	}
    }

/*     Compute TI*A*T (A is in RSF). */

    ij = ku;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j + 1;
	k = min(i__2,*n);
	dgemv_("NoTranspose", &nminr, &k, &c_b14, &ti[ti_offset], ldti, &a[j *
		 a_dim1 + 1], &c__1, &c_b49, &dwork[ij], &c__1, (ftnlen)11);
	ij += *n;
/* L80: */
    }
    dgemm_("NoTranspose", "NoTranspose", &nminr, &nminr, n, &c_b14, &dwork[ku]
	    , n, &t[t_offset], ldt, &c_b49, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);

/*     Compute TI*B and C*T. */

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
    dgemm_("NoTranspose", "NoTranspose", &nminr, m, n, &c_b14, &ti[ti_offset],
	     ldti, &dwork[ku], n, &c_b49, &b[b_offset], ldb, (ftnlen)11, (
	    ftnlen)11);

    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);
    dgemm_("NoTranspose", "NoTranspose", p, &nminr, n, &c_b14, &dwork[ku], p, 
	    &t[t_offset], ldt, &c_b49, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);

/*     Compute the singular perturbation approximation if possible. */
/*     Note that IERR = 1 on exit from AB09DD cannot appear here. */

/*     Workspace:  need real    4*(NMINR-NR); */
/*                 need integer 2*(NMINR-NR). */

    ab09dd_("Continuous", &nminr, m, p, nr, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &rcond, &iwork[1], 
	    &dwork[1], &ierr, (ftnlen)10);

/*     For discrete-time case, apply the continuous-to-discrete */
/*     bilinear transformation. */

    if (discr) {
	ab04md_("Continuous", nr, m, p, &c_b14, &c_b14, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &
		iwork[1], &dwork[1], ldwork, &ierr, (ftnlen)10);

/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[1];
	wrkopt = max(i__1,i__2);
    }
    iwork[1] = nminr;
    dwork[1] = (doublereal) wrkopt;
    dwork[2] = ricond;

    return 0;
/* *** Last line of AB09HX *** */
} /* ab09hx_ */

