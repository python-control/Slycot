/* AB09IX.f -- translated by f2c (version 20100827).
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
static doublereal c_b33 = 1.;
static doublereal c_b47 = 0.;

/* Subroutine */ int ab09ix_(char *dico, char *job, char *fact, char *ordsel, 
	integer *n, integer *m, integer *p, integer *nr, doublereal *scalec, 
	doublereal *scaleo, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ti, integer *ldti, doublereal *t, integer *ldt, 
	integer *nminr, doublereal *hsv, doublereal *tol1, doublereal *tol2, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, 
	integer *info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, 
	ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, t_dim1, t_offset, ti_dim1, ti_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, ij, ku, kv, kw, lw, ns, nr1;
    static logical bal, bta, spa;
    static integer ldw;
    static logical rsf;
    static doublereal skp;
    static integer nred;
    static doublereal atol;
    static integer ierr, ktau;
    static doublereal temp;
    extern /* Subroutine */ int ab09dd_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), ma02ad_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), mb03ud_(
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    static doublereal rcond;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrmv_(
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
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
/*     state-space representation (A,B,C,D) by using the square-root or */
/*     balancing-free square-root Balance & Truncate (B&T) or */
/*     Singular Perturbation Approximation (SPA) model reduction methods. */
/*     The computation of truncation matrices TI and T is based on */
/*     the Cholesky factor S of a controllability Grammian P = S*S' */
/*     and the Cholesky factor R of an observability Grammian Q = R'*R, */
/*     where S and R are given upper triangular matrices. */

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
/*             = 'B':  use the square-root B&T method; */
/*             = 'F':  use the balancing-free square-root B&T method; */
/*             = 'S':  use the square-root SPA method; */
/*             = 'P':  use the balancing-free square-root SPA method. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not, on entry, the matrix A is in a */
/*             real Schur form, as follows: */
/*             = 'S':  A is in a real Schur form; */
/*             = 'N':  A is a general dense square matrix. */

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
/*             The number of system outputs.  P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of */
/*             the resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(NR,NMINR), where NR */
/*             is the desired order on entry and NMINR is the number of */
/*             the Hankel singular values greater than N*EPS*S1, where */
/*             EPS is the machine precision (see LAPACK Library Routine */
/*             DLAMCH) and S1 is the largest Hankel singular value */
/*             (computed in HSV(1)); */
/*             NR can be further reduced to ensure HSV(NR) > HSV(NR+1); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,N*EPS*S1). */

/*     SCALEC  (input) DOUBLE PRECISION */
/*             Scaling factor for the Cholesky factor S of the */
/*             controllability Grammian, i.e., S/SCALEC is used to */
/*             compute the Hankel singular values.  SCALEC > 0. */

/*     SCALEO  (input) DOUBLE PRECISION */
/*             Scaling factor for the Cholesky factor R of the */
/*             observability Grammian, i.e., R/SCALEO is used to */
/*             compute the Hankel singular values.  SCALEO > 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. If FACT = 'S', */
/*             A is in a real Schur form. */
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
/*             On entry, if JOB = 'S' or JOB = 'P', the leading P-by-M */
/*             part of this array must contain the original input/output */
/*             matrix D. */
/*             On exit, if INFO = 0 and JOB = 'S' or JOB = 'P', the */
/*             leading P-by-M part of this array contains the */
/*             input/output matrix Dr of the reduced order system. */
/*             If JOB = 'B' or JOB = 'F', this array is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= 1,        if JOB = 'B' or JOB = 'F'; */
/*             LDD >= MAX(1,P), if JOB = 'S' or JOB = 'P'. */

/*     TI      (input/output) DOUBLE PRECISION array, dimension (LDTI,N) */
/*             On entry, the leading N-by-N upper triangular part of */
/*             this array must contain the Cholesky factor S of a */
/*             controllability Grammian P = S*S'. */
/*             On exit, if INFO = 0, and NR > 0, the leading NMINR-by-N */
/*             part of this array contains the left truncation matrix */
/*             TI in (1), for the B&T approach, or in (2), for the */
/*             SPA approach. */

/*     LDTI    INTEGER */
/*             The leading dimension of array TI.  LDTI >= MAX(1,N). */

/*     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             On entry, the leading N-by-N upper triangular part of */
/*             this array must contain the Cholesky factor R of an */
/*             observability Grammian Q = R'*R. */
/*             On exit, if INFO = 0, and NR > 0, the leading N-by-NMINR */
/*             part of this array contains the right truncation matrix */
/*             T in (1), for the B&T approach, or in (2), for the */
/*             SPA approach. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     NMINR   (output) INTEGER */
/*             The number of Hankel singular values greater than */
/*             MAX(TOL2,N*EPS*S1). */
/*             Note: If S and R are the Cholesky factors of the */
/*             controllability and observability Grammians of the */
/*             original system (A,B,C,D), respectively, then NMINR is */
/*             the order of a minimal realization of the original system. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the Hankel singular values, */
/*             ordered decreasingly. The Hankel singular values are */
/*             singular values of the product R*S. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of the reduced system. */
/*             For model reduction, the recommended value lies in the */
/*             interval [0.00001,0.001]. */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = N*EPS*S1, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH) and S1 is the largest */
/*             Hankel singular value (computed in HSV(1)). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the system. */
/*             The recommended value is TOL2 = N*EPS*S1. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0, and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension LIWORK, where */
/*             LIWORK = 0,   if JOB = 'B'; */
/*             LIWORK = N,   if JOB = 'F'; */
/*             LIWORK = 2*N, if JOB = 'S' or 'P'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, 2*N*N + 5*N, N*MAX(M,P) ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than NMINR, the order of a minimal realization of */
/*                   the given system; in this case, the resulting NR is */
/*                   set automatically to NMINR; */
/*             = 2:  with ORDSEL = 'F', the selected order NR corresponds */
/*                   to repeated singular values, which are neither all */
/*                   included nor all excluded from the reduced model; */
/*                   in this case, the resulting NR is set automatically */
/*                   to the largest value such that HSV(NR) > HSV(NR+1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                             (3) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09IX determines for */
/*     the given system (3), the matrices of a reduced NR order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                         (4) */

/*     by using the square-root or balancing-free square-root */
/*     Balance & Truncate (B&T) or Singular Perturbation Approximation */
/*     (SPA) model reduction methods. */

/*     The projection matrices TI and T are determined using the */
/*     Cholesky factors S and R of a controllability Grammian P and an */
/*     observability Grammian Q. */
/*     The Hankel singular values HSV(1), ...., HSV(N) are computed as */
/*     singular values of the product R*S. */

/*     If JOB = 'B', the square-root Balance & Truncate technique */
/*     of [1] is used. */

/*     If JOB = 'F', the balancing-free square-root version of the */
/*     Balance & Truncate technique [2] is used. */

/*     If JOB = 'S', the square-root version of the Singular Perturbation */
/*     Approximation method [3,4] is used. */

/*     If JOB = 'P', the balancing-free square-root version of the */
/*     Singular Perturbation Approximation method [3,4] is used. */

/*     REFERENCES */

/*     [1] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [2] Varga A. */
/*         Efficient minimal realization procedure based on balancing. */
/*         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991, */
/*         A. El Moudni, P. Borne, S. G. Tzafestas (Eds.), */
/*         Vol. 2, pp. 42-46. */

/*     [3] Liu Y. and Anderson B.D.O. */
/*         Singular Perturbation Approximation of balanced systems. */
/*         Int. J. Control, Vol. 50, pp. 1379-1405, 1989. */

/*     [4] Varga A. */
/*         Balancing-free square-root algorithm for computing singular */
/*         perturbation approximations. */
/*         Proc. 30-th CDC, Brighton, Dec. 11-13, 1991, */
/*         Vol. 2, pp. 1062-1065. */

/*     NUMERICAL ASPECTS */

/*     The implemented method relies on accuracy enhancing square-root */
/*     or balancing-free square-root methods. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*              Sep. 2001. */

/*     KEYWORDS */

/*     Balance and truncate, minimal state-space representation, */
/*     model reduction, multivariable system, */
/*     singular perturbation approximation, state-space model. */

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
    ti_dim1 = *ldti;
    ti_offset = 1 + ti_dim1;
    ti -= ti_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --hsv;
    --iwork;
    --dwork;

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
    rsf = lsame_(fact, "S", (ftnlen)1, (ftnlen)1);
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/* Computing MAX */
    i__1 = 1, i__2 = (*n << 1) * *n + *n * 5, i__1 = max(i__1,i__2), i__2 = *
	    n * max(*m,*p);
    lw = max(i__1,i__2);

/*     Test the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (bta || spa)) {
	*info = -2;
    } else if (! (rsf || lsame_(fact, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*p < 0) {
	*info = -7;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -8;
    } else if (*scalec <= 0.) {
	*info = -9;
    } else if (*scaleo <= 0.) {
	*info = -10;
    } else if (*lda < max(1,*n)) {
	*info = -12;
    } else if (*ldb < max(1,*n)) {
	*info = -14;
    } else if (*ldc < max(1,*p)) {
	*info = -16;
    } else if (*ldd < 1 || spa && *ldd < *p) {
	*info = -18;
    } else if (*ldti < max(1,*n)) {
	*info = -20;
    } else if (*ldt < max(1,*n)) {
	*info = -22;
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
	*info = -26;
    } else if (*ldwork < lw) {
	*info = -29;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09IX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nr = 0;
	*nminr = 0;
	dwork[1] = 1.;
	return 0;
    }

/*     Save S in DWORK(KV). */

    kv = 1;
    ku = kv + *n * *n;
    kw = ku + *n * *n;
    dlacpy_("Upper", n, n, &ti[ti_offset], ldti, &dwork[kv], n, (ftnlen)5);
/*                             | x x | */
/*     Compute R*S in the form | 0 x | in TI. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dtrmv_("Upper", "NoTranspose", "NonUnit", &j, &t[t_offset], ldt, &ti[
		j * ti_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
/* L10: */
    }

/*     Compute the singular value decomposition R*S = V*Sigma*UT of the */
/*     upper triangular matrix R*S, with UT in TI and V in DWORK(KU). */

/*     Workspace:  need   2*N*N + 5*N; */
/*                 prefer larger. */

    i__1 = *ldwork - kw + 1;
    mb03ud_("Vectors", "Vectors", n, &ti[ti_offset], ldti, &dwork[ku], n, &
	    hsv[1], &dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)7);
    if (ierr != 0) {
	*info = 1;
	return 0;
    }
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Scale the singular values. */

    d__1 = 1. / *scalec / *scaleo;
    dscal_(n, &d__1, &hsv[1], &c__1);

/*     Partition Sigma, U and V conformally as: */

/*     Sigma = diag(Sigma1,Sigma2,Sigma3),  U = [U1,U2,U3] (U' in TI) and */
/*     V = [V1,V2,V3] (in DWORK(KU)). */

/*     Compute NMINR, the order of a minimal realization, as the order */
/*     of [Sigma1 Sigma2]. */

    toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
/* Computing MAX */
    d__1 = *tol2, d__2 = toldef * hsv[1];
    atol = max(d__1,d__2);
    *nminr = *n;
L20:
    if (*nminr > 0) {
	if (hsv[*nminr] <= atol) {
	    --(*nminr);
	    goto L20;
	}
    }

/*     Compute the order NR of reduced system, as the order of Sigma1. */

    if (fixord) {

/*        Check if the desired order is less than the order of a minimal */
/*        realization. */

	if (*nr > *nminr) {

/*           Reduce the order to NMINR. */

	    *nr = *nminr;
	    *iwarn = 1;
	}

/*        Check for singular value multiplicity at cut-off point. */

	if (*nr > 0 && *nr < *nminr) {
	    skp = hsv[*nr];
	    if (skp - hsv[*nr + 1] <= toldef * skp) {
		*iwarn = 2;

/*              Reduce the order such that HSV(NR) > HSV(NR+1). */

L30:
		--(*nr);
		if (*nr > 0) {
		    if (hsv[*nr] - skp <= toldef * skp) {
			goto L30;
		    }
		}
	    }
	}
    } else {

/*        The order is given as the number of singular values */
/*        exceeding MAX( TOL1, N*EPS*HSV(1) ). */

	atol = max(*tol1,atol);
	*nr = 0;
	i__1 = *nminr;
	for (j = 1; j <= i__1; ++j) {
	    if (hsv[j] <= atol) {
		goto L50;
	    }
	    ++(*nr);
/* L40: */
	}
L50:
	;
    }

/*     Finish if the order is zero. */

    if (*nr == 0) {
	if (spa) {
	    ab09dd_(dico, n, m, p, nr, &a[a_offset], lda, &b[b_offset], ldb, &
		    c__[c_offset], ldc, &d__[d_offset], ldd, &rcond, &iwork[1]
		    , &dwork[1], &ierr, (ftnlen)1);
	}
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

/*     Compute NS, the order of Sigma2. For BTA, NS = 0. */

    if (spa) {
	nred = *nminr;
    } else {
	nred = *nr;
    }
    ns = nred - *nr;

/*     Compute the truncation matrices. */

/*     Compute TI' = | TI1' TI2' | = R'*| V1 V2 | in DWORK(KU). */

    dtrmm_("Left", "Upper", "Transpose", "NonUnit", n, &nred, &c_b33, &t[
	    t_offset], ldt, &dwork[ku], n, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute  T = | T1 T2 | = S*| U1 U2 | . */

    ma02ad_("Full", &nred, n, &ti[ti_offset], ldti, &t[t_offset], ldt, (
	    ftnlen)4);
    dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, &nred, &c_b33, &
	    dwork[kv], n, &t[t_offset], ldt, (ftnlen)4, (ftnlen)5, (ftnlen)11,
	     (ftnlen)7);

    ktau = kw;
    if (bal) {
	ij = ku;

/*        Square-Root B&T/SPA method. */

/*        Compute the truncation matrices for balancing */
/*                        -1/2                -1/2 */
/*               T1*Sigma1     and TI1'*Sigma1    . */

	i__1 = *nr;
	for (j = 1; j <= i__1; ++j) {
	    temp = 1. / sqrt(hsv[j]);
	    dscal_(n, &temp, &t[j * t_dim1 + 1], &c__1);
	    dscal_(n, &temp, &dwork[ij], &c__1);
	    ij += *n;
/* L60: */
	}

    } else {

/*        Balancing-Free B&T/SPA method. */

/*        Compute orthogonal bases for the images of matrices T1 and */
/*        TI1'. */

/*        Workspace:  need   2*N*N + 2*N; */
/*                    prefer larger. */

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

/*        Workspace:  need   2*N*N + 2*N; */
/*                    prefer larger. */

	nr1 = *nr + 1;
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

    ma02ad_("Full", n, &nred, &dwork[ku], n, &ti[ti_offset], ldti, (ftnlen)4);

    if (! bal) {
/*                        -1 */
/*        Compute (TI1*T1)  *TI1 in TI. */

	dgemm_("NoTranspose", "NoTranspose", nr, nr, n, &c_b33, &ti[ti_offset]
		, ldti, &t[t_offset], ldt, &c_b47, &dwork[ku], n, (ftnlen)11, 
		(ftnlen)11);
	dgetrf_(nr, nr, &dwork[ku], n, &iwork[1], &ierr);
	dgetrs_("NoTranspose", nr, n, &dwork[ku], n, &iwork[1], &ti[ti_offset]
		, ldti, &ierr, (ftnlen)11);

	if (ns > 0) {
/*                           -1 */
/*           Compute (TI2*T2)  *TI2 in TI2. */

	    dgemm_("NoTranspose", "NoTranspose", &ns, &ns, n, &c_b33, &ti[nr1 
		    + ti_dim1], ldti, &t[nr1 * t_dim1 + 1], ldt, &c_b47, &
		    dwork[ku], n, (ftnlen)11, (ftnlen)11);
	    dgetrf_(&ns, &ns, &dwork[ku], n, &iwork[1], &ierr);
	    dgetrs_("NoTranspose", &ns, n, &dwork[ku], n, &iwork[1], &ti[nr1 
		    + ti_dim1], ldti, &ierr, (ftnlen)11);
	}
    }

/*     Compute TI*A*T. Exploit RSF of A if possible. */
/*     Workspace:  need   N*N. */

    if (rsf) {
	ij = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__2 = j + 1;
	    k = min(i__2,*n);
	    dgemv_("NoTranspose", &nred, &k, &c_b33, &ti[ti_offset], ldti, &a[
		    j * a_dim1 + 1], &c__1, &c_b47, &dwork[ij], &c__1, (
		    ftnlen)11);
	    ij += *n;
/* L80: */
	}
    } else {
	dgemm_("NoTranspose", "NoTranspose", &nred, n, n, &c_b33, &ti[
		ti_offset], ldti, &a[a_offset], lda, &c_b47, &dwork[1], n, (
		ftnlen)11, (ftnlen)11);
    }
    dgemm_("NoTranspose", "NoTranspose", &nred, &nred, n, &c_b33, &dwork[1], 
	    n, &t[t_offset], ldt, &c_b47, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);

/*     Compute TI*B and C*T. */
/*     Workspace:  need   N*MAX(M,P). */

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
    dgemm_("NoTranspose", "NoTranspose", &nred, m, n, &c_b33, &ti[ti_offset], 
	    ldti, &dwork[1], n, &c_b47, &b[b_offset], ldb, (ftnlen)11, (
	    ftnlen)11);

    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], p, (ftnlen)4);
    dgemm_("NoTranspose", "NoTranspose", p, &nred, n, &c_b33, &dwork[1], p, &
	    t[t_offset], ldt, &c_b47, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);

/*     Compute the singular perturbation approximation if possible. */
/*     Note that IERR = 1 on exit from AB09DD cannot appear here. */

/*     Workspace:  need real    4*(NMINR-NR); */
/*                 need integer 2*(NMINR-NR). */

    if (spa) {
	ab09dd_(dico, &nred, m, p, nr, &a[a_offset], lda, &b[b_offset], ldb, &
		c__[c_offset], ldc, &d__[d_offset], ldd, &rcond, &iwork[1], &
		dwork[1], &ierr, (ftnlen)1);
    } else {
	*nminr = *nr;
    }
    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of AB09IX *** */
} /* ab09ix_ */

