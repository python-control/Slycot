/* AB09GD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ab09gd_(char *dico, char *jobcf, char *fact, char *jobmr,
	 char *equil, char *ordsel, integer *n, integer *m, integer *p, 
	integer *nr, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *nq, doublereal *hsv, 
	doublereal *tol1, doublereal *tol2, doublereal *tol3, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *iwarn, integer *info, 
	ftnlen dico_len, ftnlen jobcf_len, ftnlen fact_len, ftnlen jobmr_len, 
	ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer kb, kc, kd, mp, pm, kt, kw, lw1, lw2, lw3, lw4, kbr, kcr, 
	    kbt, kdr, kdt, ndr, kti, lwr;
    static logical left;
    static integer ierr;
    extern /* Subroutine */ int sb08cd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), ab09bx_(
	    char *, char *, char *, integer *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen), sb08dd_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), sb08ed_(
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, ftnlen), sb08fd_(char *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, ftnlen), sb08gd_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *), sb08hd_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), tb01id_(char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static logical stabd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static integer maxmp, nminr;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnk, wrkopt;


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
/*     state-space representation (A,B,C,D) by using either the */
/*     square-root or the balancing-free square-root Singular */
/*     Perturbation Approximation (SPA) model reduction method in */
/*     conjunction with stable coprime factorization techniques. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBCF   CHARACTER*1 */
/*             Specifies whether left or right coprime factorization is */
/*             to be used as follows: */
/*             = 'L':  use left coprime factorization; */
/*             = 'R':  use right coprime factorization. */

/*     FACT    CHARACTER*1 */
/*             Specifies the type of coprime factorization to be computed */
/*             as follows: */
/*             = 'S':  compute a coprime factorization with prescribed */
/*                     stability degree ALPHA; */
/*             = 'I':  compute a coprime factorization with inner */
/*                     denominator. */

/*     JOBMR   CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'N':  use the balancing-free square-root */
/*                     Balance & Truncate method. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

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
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(NR,NQ,NMIN), where NR */
/*             is the desired order on entry, NQ is the order of the */
/*             computed coprime factorization of the given system, and */
/*             NMIN is the order of a minimal realization of the extended */
/*             system (see METHOD); NMIN is determined as the number of */
/*             Hankel singular values greater than NQ*EPS*HNORM(Ge), */
/*             where EPS is the machine precision (see LAPACK Library */
/*             Routine DLAMCH) and HNORM(Ge) is the Hankel norm of the */
/*             extended system (computed in HSV(1)); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,NQ*EPS*HNORM(Ge)). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             If FACT = 'S', the desired stability degree for the */
/*             factors of the coprime factorization (see SLICOT Library */
/*             routines SB08ED/SB08FD). */
/*             ALPHA < 0 for a continuous-time system (DICO = 'C'), and */
/*             0 <= ALPHA < 1 for a discrete-time system (DICO = 'D'). */
/*             If FACT = 'I', ALPHA is not used. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the reduced */
/*             order system. */

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

/*     NQ      (output) INTEGER */
/*             The order of the computed extended system Ge (see METHOD). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the NQ Hankel singular values of */
/*             the extended system Ge ordered decreasingly (see METHOD). */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced extended system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(Ge), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(Ge) is the */
/*             Hankel-norm of the extended system (computed in HSV(1)). */
/*             The value TOL1 = NQ*EPS*HNORM(Ge) is used by default if */
/*             TOL1 <= 0 on entry, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the extended system Ge (see METHOD). */
/*             The recommended value is TOL2 = NQ*EPS*HNORM(Ge). */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0, then TOL2 <= TOL1. */

/*     TOL3    DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of */
/*             B or C are considered zero (used for controllability or */
/*             observability tests by the coprime factorization method). */
/*             If the user sets TOL3 <= 0, then an implicitly computed, */
/*             default tolerance TOLDEF is used: */
/*             TOLDEF = N*EPS*NORM(C'), if JOBCF = 'L', or */
/*             TOLDEF = N*EPS*NORM(B),  if JOBCF = 'R', */
/*             where EPS is the machine precision, and NORM(.) denotes */
/*             the 1-norm of a matrix. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX(1,2*N,PM)) */
/*             where  PM = P, if JOBCF = 'L', */
/*                    PM = M, if JOBCF = 'R'. */
/*             On exit with INFO = 0, IWORK(1) contains the order of the */
/*             minimal realization of the system. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,LW1) if JOBCF = 'L' and FACT = 'S', */
/*             LDWORK >= MAX(1,LW2) if JOBCF = 'L' and FACT = 'I', */
/*             LDWORK >= MAX(1,LW3) if JOBCF = 'R' and FACT = 'S', */
/*             LDWORK >= MAX(1,LW4) if JOBCF = 'R' and FACT = 'I', where */
/*             LW1 = N*(2*MAX(M,P) + P) + MAX(M,P)*(MAX(M,P) + P) + */
/*                   MAX( N*P+MAX(N*(N+5), 5*P, 4*M), LWR ), */
/*             LW2 = N*(2*MAX(M,P) + P) + MAX(M,P)*(MAX(M,P) + P) + */
/*                   MAX( N*P+MAX(N*(N+5), P*(P+2), 4*P, 4*M), LWR ), */
/*             LW3 = (N+M)*(M+P) + MAX( 5*M, 4*P, LWR ), */
/*             LW4 = (N+M)*(M+P) + MAX( M*(M+2), 4*M, 4*P, LWR ), and */
/*             LWR = 2*N*N + N*(MAX(N,M+P)+5) + N*(N+1)/2. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 10*K+I: */
/*               I = 1:  with ORDSEL = 'F', the selected order NR is */
/*                       greater than the order of the computed coprime */
/*                       factorization of the given system. In this case, */
/*                       the resulting NR is set automatically to a value */
/*                       corresponding to the order of a minimal */
/*                       realization of the system; */
/*               K > 0:  K violations of the numerical stability */
/*                       condition occured when computing the coprime */
/*                       factorization using pole assignment (see SLICOT */
/*                       Library routines SB08CD/SB08ED, SB08DD/SB08FD). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to a real Schur form failed; */
/*             = 2:  a failure was detected during the ordering of the */
/*                   real Schur form of A, or in the iterative process */
/*                   for reordering the eigenvalues of Z'*(A + H*C)*Z */
/*                   (or Z'*(A + B*F)*Z) along the diagonal; see SLICOT */
/*                   Library routines SB08CD/SB08ED (or SB08DD/SB08FD); */
/*             = 3:  the matrix A has an observable or controllable */
/*                   eigenvalue on the imaginary axis if DICO = 'C' or */
/*                   on the unit circle if DICO = 'D'; */
/*             = 4:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                       (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system, and let G be the corresponding */
/*     transfer-function matrix. The subroutine AB09GD determines */
/*     the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                   (2) */

/*     with the transfer-function matrix Gr, by using the */
/*     singular perturbation approximation (SPA) method in conjunction */
/*     with a left coprime factorization (LCF) or a right coprime */
/*     factorization (RCF) technique: */

/*     1. Compute the appropriate stable coprime factorization of G: */
/*                     -1                   -1 */
/*                G = R  *Q (LCF) or G = Q*R   (RCF). */

/*     2. Perform the model reduction algorithm on the extended system */
/*                                           ( Q ) */
/*                Ge = ( Q R ) (LCF) or Ge = ( R )  (RCF) */

/*        to obtain a reduced extended system with reduced factors */
/*                                               ( Qr ) */
/*                Ger = ( Qr Rr ) (LCF) or Ger = ( Rr )  (RCF). */

/*     3. Recover the reduced system from the reduced factors as */
/*                       -1                       -1 */
/*                Gr = Rr  *Qr (LCF) or Gr = Qr*Rr   (RCF). */

/*     The approximation error for the extended system satisfies */

/*        HSV(NR) <= INFNORM(Ge-Ger) <= 2*[HSV(NR+1) + ... + HSV(NQ)], */

/*     where INFNORM(G) is the infinity-norm of G. */

/*     If JOBMR = 'B', the balancing-based square-root SPA method of [1] */
/*     is used for model reduction. */
/*     If JOBMR = 'N', the balancing-free square-root SPA method of [2] */
/*     is used for model reduction. */
/*     By setting TOL1 = TOL2, the routine can be used to compute */
/*     Balance & Truncate approximations. */

/*     If FACT = 'S', the stable coprime factorization with prescribed */
/*     stability degree ALPHA is computed by using the algorithm of [3]. */
/*     If FACT = 'I', the stable coprime factorization with inner */
/*     denominator is computed by using the algorithm of [4]. */

/*     REFERENCES */

/*     [1] Liu Y. and Anderson B.D.O. */
/*         Singular Perturbation Approximation of Balanced Systems. */
/*         Int. J. Control, Vol. 50, pp. 1379-1405, 1989. */

/*     [2] Varga A. */
/*         Balancing-free square-root algorithm for computing singular */
/*         perturbation approximations. */
/*         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991, Vol. 2, */
/*         pp. 1062-1065. */

/*     [3] Varga A. */
/*         Coprime factors model reduction method based on square-root */
/*         balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, Vol. 11, */
/*         pp. 303-311, 1993. */

/*     [4] Varga A. */
/*         A Schur method for computing coprime factorizations with */
/*         inner denominators and applications in model reduction. */
/*         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, August 1998. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */

/*     KEYWORDS */

/*     Balancing, coprime factorization, minimal realization, */
/*     model reduction, multivariable system, state-space model. */

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
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    *iwarn = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
    left = lsame_(jobcf, "L", (ftnlen)1, (ftnlen)1);
    stabd = lsame_(fact, "S", (ftnlen)1, (ftnlen)1);
    maxmp = max(*m,*p);

/* Computing MAX */
    i__1 = *n, i__2 = *m + *p;
    lwr = (*n << 1) * *n + *n * (max(i__1,i__2) + 5) + *n * (*n + 1) / 2;
    lw1 = *n * ((maxmp << 1) + *p) + maxmp * (maxmp + *p);
/* Computing MAX */
/* Computing MAX */
    i__2 = *n * (*n + 5), i__3 = *p * (*p + 2), i__2 = max(i__2,i__3), i__3 = 
	    *p << 2, i__2 = max(i__2,i__3), i__3 = *m << 2;
    i__1 = *n * *p + max(i__2,i__3);
    lw2 = lw1 + max(i__1,lwr);
/* Computing MAX */
/* Computing MAX */
    i__2 = *n * (*n + 5), i__3 = *p * 5, i__2 = max(i__2,i__3), i__3 = *m << 
	    2;
    i__1 = *n * *p + max(i__2,i__3);
    lw1 += max(i__1,lwr);
/* Computing MAX */
    i__1 = *m * 5, i__2 = *p << 2, i__1 = max(i__1,i__2);
    lw3 = (*n + *m) * (*m + *p) + max(i__1,lwr);
/* Computing MAX */
    i__1 = *m * (*m + 2), i__2 = *m << 2, i__1 = max(i__1,i__2), i__2 = *p << 
	    2, i__1 = max(i__1,i__2);
    lw4 = (*n + *m) * (*m + *p) + max(i__1,lwr);

/*     Test the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (left || lsame_(jobcf, "R", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (stabd || lsame_(fact, "I", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -5;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -6;
    } else if (stabd && (! discr && *alpha >= 0. || discr && (*alpha < 0. || *
	    alpha >= 1.))) {
	*info = -7;
    } else if (*n < 0) {
	*info = -8;
    } else if (*m < 0) {
	*info = -9;
    } else if (*p < 0) {
	*info = -10;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -11;
    } else if (*lda < max(1,*n)) {
	*info = -13;
    } else if (*ldb < max(1,*n)) {
	*info = -15;
    } else if (*ldc < max(1,*p)) {
	*info = -17;
    } else if (*ldd < max(1,*p)) {
	*info = -19;
    } else if (*tol2 > 0. && *tol2 > *tol1) {
	*info = -23;
    } else if (*ldwork < 1 || stabd && left && *ldwork < lw1 || ! stabd && 
	    left && *ldwork < lw2 || stabd && ! left && *ldwork < lw3 || ! 
	    stabd && ! left && *ldwork < lw4) {
	*info = -27;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09GD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nr = 0;
	*nq = 0;
	iwork[1] = 0;
	dwork[1] = 1.;
	return 0;
    }

    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */

	maxred = 100.;
	tb01id_("A", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb, &
		c__[c_offset], ldc, &dwork[1], info, (ftnlen)1);
    }

/*     Perform the coprime factor model reduction procedure. */

    kd = 1;
    if (left) {
/*                           -1 */
/*        Compute a LCF G = R  *Q. */

	mp = *m + *p;
	kdr = kd + maxmp * maxmp;
	kc = kdr + maxmp * *p;
	kb = kc + maxmp * *n;
	kbr = kb + *n * maxmp;
	kw = kbr + *n * *p;
	lwr = *ldwork - kw + 1;
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[kb], n, (ftnlen)4);
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[kc], &maxmp, (
		ftnlen)4);
	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[kd], &maxmp, (
		ftnlen)4);

	if (stabd) {

/*           Compute a LCF with prescribed stability degree. */

/*           Workspace needed:      N*(2*MAX(M,P)+P) + */
/*                                  MAX(M,P)*(MAX(M,P)+P); */
/*           Additional workspace:  need   N*P+MAX(N*(N+5),5*P,4*M); */
/*                                  prefer larger. */

	    sb08ed_(dico, n, m, p, alpha, &a[a_offset], lda, &dwork[kb], n, &
		    dwork[kc], &maxmp, &dwork[kd], &maxmp, nq, &ndr, &dwork[
		    kbr], n, &dwork[kdr], &maxmp, tol3, &dwork[kw], &lwr, 
		    iwarn, info, (ftnlen)1);
	} else {

/*           Compute a LCF with inner denominator. */

/*           Workspace needed:      N*(2*MAX(M,P)+P) + */
/*                                  MAX(M,P)*(MAX(M,P)+P); */
/*           Additional workspace:  need   N*P + */
/*                                         MAX(N*(N+5),P*(P+2),4*P,4*M); */
/*                                  prefer larger. */

	    sb08cd_(dico, n, m, p, &a[a_offset], lda, &dwork[kb], n, &dwork[
		    kc], &maxmp, &dwork[kd], &maxmp, nq, &ndr, &dwork[kbr], n,
		     &dwork[kdr], &maxmp, tol3, &dwork[kw], &lwr, iwarn, info,
		     (ftnlen)1);
	}

	*iwarn *= 10;
	if (*info != 0) {
	    return 0;
	}

	wrkopt = (integer) dwork[kw] + kw - 1;

	if (*nq == 0) {
	    *nr = 0;
	    iwork[1] = 0;
	    dwork[1] = (doublereal) wrkopt;
	    return 0;
	}

	if (maxmp > *m) {

/*           Form the matrices ( BQ, BR ) and ( DQ, DR ) in consecutive */
/*           columns (see SLICOT Library routines SB08CD/SB08ED). */

	    kbt = kbr;
	    kbr = kb + *n * *m;
	    kdt = kdr;
	    kdr = kd + maxmp * *m;
	    dlacpy_("Full", nq, p, &dwork[kbt], n, &dwork[kbr], n, (ftnlen)4);
	    dlacpy_("Full", p, p, &dwork[kdt], &maxmp, &dwork[kdr], &maxmp, (
		    ftnlen)4);
	}

/*        Perform model reduction on ( Q, R ) to determine ( Qr, Rr ). */

/*        Workspace needed:      N*(2*MAX(M,P)+P) + */
/*                               MAX(M,P)*(MAX(M,P)+P) + 2*N*N; */
/*        Additional workspace:  need   N*(MAX(N,M+P)+5) + N*(N+1)/2; */
/*                               prefer larger. */

	kt = kw;
	kti = kt + *nq * *nq;
	kw = kti + *nq * *nq;
	i__1 = *ldwork - kw + 1;
	ab09bx_(dico, jobmr, ordsel, nq, &mp, p, nr, &a[a_offset], lda, &
		dwork[kb], n, &dwork[kc], &maxmp, &dwork[kd], &maxmp, &hsv[1],
		 &dwork[kt], n, &dwork[kti], n, tol1, tol2, &iwork[1], &dwork[
		kw], &i__1, &iwarnk, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

	*iwarn += iwarnk;
	if (ierr != 0) {
	    *info = 4;
	    return 0;
	}

	nminr = iwork[1];
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
/*                                                -1 */
/*        Compute the reduced order system Gr = Rr  *Qr. */

/*        Workspace needed:      N*(2*MAX(M,P)+P) + */
/*                               MAX(M,P)*(MAX(M,P)+P); */
/*        Additional workspace:  need   4*P. */

	kw = kt;
	sb08gd_(nr, m, p, &a[a_offset], lda, &dwork[kb], n, &dwork[kc], &
		maxmp, &dwork[kd], &maxmp, &dwork[kbr], n, &dwork[kdr], &
		maxmp, &iwork[1], &dwork[kw], info);

/*        Copy the reduced system matrices Br, Cr, and Dr to B, C, and D, */
/*        respectively. */

	dlacpy_("Full", nr, m, &dwork[kb], n, &b[b_offset], ldb, (ftnlen)4);
	dlacpy_("Full", p, nr, &dwork[kc], &maxmp, &c__[c_offset], ldc, (
		ftnlen)4);
	dlacpy_("Full", p, m, &dwork[kd], &maxmp, &d__[d_offset], ldd, (
		ftnlen)4);
    } else {
/*                             -1 */
/*        Compute a RCF G = Q*R  . */

	pm = *p + *m;
	kdr = kd + *p;
	kc = kd + pm * *m;
	kcr = kc + *p;
	kw = kc + pm * *n;
	lwr = *ldwork - kw + 1;
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[kc], &pm, (ftnlen)4)
		;
	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[kd], &pm, (ftnlen)4)
		;

	if (stabd) {

/*           Compute a RCF with prescribed stability degree. */

/*           Workspace needed:      (N+M)*(M+P); */
/*           Additional workspace:  need  MAX( N*(N+5), 5*M, 4*P ); */
/*                                  prefer larger. */

	    sb08fd_(dico, n, m, p, alpha, &a[a_offset], lda, &b[b_offset], 
		    ldb, &dwork[kc], &pm, &dwork[kd], &pm, nq, &ndr, &dwork[
		    kcr], &pm, &dwork[kdr], &pm, tol3, &dwork[kw], &lwr, 
		    iwarn, info, (ftnlen)1);
	} else {

/*           Compute a RCF with inner denominator. */

/*           Workspace needed:      (N+M)*(M+P); */
/*           Additional workspace:  need  MAX(N*(N+5),M*(M+2),4*M,4*P); */
/*                                  prefer larger. */

	    sb08dd_(dico, n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &
		    dwork[kc], &pm, &dwork[kd], &pm, nq, &ndr, &dwork[kcr], &
		    pm, &dwork[kdr], &pm, tol3, &dwork[kw], &lwr, iwarn, info,
		     (ftnlen)1);
	}

	*iwarn *= 10;
	if (*info != 0) {
	    return 0;
	}

	wrkopt = (integer) dwork[kw] + kw - 1;

	if (*nq == 0) {
	    *nr = 0;
	    iwork[1] = 0;
	    dwork[1] = (doublereal) wrkopt;
	    return 0;
	}
/*                                   ( Q )              ( Qr ) */
/*        Perform model reduction on ( R ) to determine ( Rr ). */

/*        Workspace needed:      (N+M)*(M+P) + 2*N*N; */
/*        Additional workspace:  need   N*(MAX(N,M+P)+5) + N*(N+1)/2; */
/*                               prefer larger. */

	kt = kw;
	kti = kt + *nq * *nq;
	kw = kti + *nq * *nq;
	i__1 = *ldwork - kw + 1;
	ab09bx_(dico, jobmr, ordsel, nq, m, &pm, nr, &a[a_offset], lda, &b[
		b_offset], ldb, &dwork[kc], &pm, &dwork[kd], &pm, &hsv[1], &
		dwork[kt], n, &dwork[kti], n, tol1, tol2, &iwork[1], &dwork[
		kw], &i__1, &iwarnk, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

	*iwarn += iwarnk;
	if (ierr != 0) {
	    *info = 4;
	    return 0;
	}

	nminr = iwork[1];
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
/*                                                   -1 */
/*        Compute the reduced order system Gr = Qr*Rr  . */

/*        Workspace needed:      (N+M)*(M+P); */
/*        Additional workspace:  need 4*M. */

	kw = kt;
	sb08hd_(nr, m, p, &a[a_offset], lda, &b[b_offset], ldb, &dwork[kc], &
		pm, &dwork[kd], &pm, &dwork[kcr], &pm, &dwork[kdr], &pm, &
		iwork[1], &dwork[kw], info);

/*        Copy the reduced system matrices Cr and Dr to C and D. */

	dlacpy_("Full", p, nr, &dwork[kc], &pm, &c__[c_offset], ldc, (ftnlen)
		4);
	dlacpy_("Full", p, m, &dwork[kd], &pm, &d__[d_offset], ldd, (ftnlen)4)
		;
    }

    iwork[1] = nminr;
    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of AB09GD *** */
} /* ab09gd_ */

