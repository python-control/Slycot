/* SB16AY.f -- translated by f2c (version 20100827).
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

static logical c_true = TRUE_;
static logical c_false = FALSE_;
static doublereal c_b21 = 0.;
static doublereal c_b24 = 1.;
static doublereal c_b28 = -1.;
static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int sb16ay_(char *dico, char *jobc, char *jobo, char *weight,
	 integer *n, integer *m, integer *p, integer *nc, integer *ncs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *ac, 
	integer *ldac, doublereal *bc, integer *ldbc, doublereal *cc, integer 
	*ldcc, doublereal *dc, integer *lddc, doublereal *scalec, doublereal *
	scaleo, doublereal *s, integer *lds, doublereal *r__, integer *ldr, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen jobc_len, ftnlen jobo_len, ftnlen weight_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, ac_dim1, ac_offset, b_dim1, b_offset, bc_dim1, 
	    bc_offset, c_dim1, c_offset, cc_dim1, cc_offset, d_dim1, d_offset,
	     dc_dim1, dc_offset, r_dim1, r_offset, s_dim1, s_offset, i__1, 
	    i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer me, ne, jj, ki, pe, kl, kq, kr, mp, ku, kw, lw, nnc, kwa, 
	    kwb, kwc, kwd, ldu, ncu;
    static doublereal dum[1], tol;
    static integer ncu1;
    static logical perf;
    static integer ierr, nncu, ktau;
    extern /* Subroutine */ int ab05pd_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), ab05qd_(char *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen), ab07nd_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    static integer mbbar;
    extern /* Subroutine */ int mb04od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    static integer pcbar;
    extern /* Subroutine */ int mb01wd_(char *, char *, char *, char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen), sb03od_(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dcopy_(integer *, doublereal *, integer *,
	     doublereal *, integer *);
    static logical leftw;
    extern /* Subroutine */ int dsyev_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static char jobfac[1];
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical frwght, rightw;
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

/*     To compute for given state-space representations (A,B,C,D) and */
/*     (Ac,Bc,Cc,Dc) of the transfer-function matrices of the */
/*     open-loop system G and feedback controller K, respectively, */
/*     the Cholesky factors of the frequency-weighted */
/*     controllability and observability Grammians corresponding */
/*     to a frequency-weighted model reduction problem. */
/*     The controller must stabilize the closed-loop system. */
/*     The state matrix Ac must be in a block-diagonal real Schur form */
/*     Ac = diag(Ac1,Ac2), where Ac1 contains the unstable eigenvalues */
/*     of Ac and Ac2 contains the stable eigenvalues of Ac. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G and K are continuous-time systems; */
/*             = 'D':  G and K are discrete-time systems. */

/*     JOBC    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted controllability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to standard Enns' method [1]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified Enns' method of [2]. */

/*     JOBO    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted observability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to standard Enns' method [1]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [2]. */

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency-weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'O':  stability enforcing left (output) weighting */
/*                               -1 */
/*                     V = (I-G*K) *G is used (W = I); */
/*             = 'I':  stability enforcing right (input) weighting */
/*                               -1 */
/*                     W = (I-G*K) *G is used (V = I); */
/*             = 'P':  stability and performance enforcing weightings */
/*                               -1                -1 */
/*                     V = (I-G*K) *G ,  W = (I-G*K)  are used. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the open-loop system state-space */
/*             representation, i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NC      (input) INTEGER */
/*             The order of the controller state-space representation, */
/*             i.e., the order of the matrix AC.  NC >= 0. */

/*     NCS     (input) INTEGER */
/*             The dimension of the stable part of the controller, i.e., */
/*             the order of matrix Ac2.  NC >= NCS >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system with the transfer-function */
/*             matrix G. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             input/output matrix D of the open-loop system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     AC      (input) DOUBLE PRECISION array, dimension (LDAC,NC) */
/*             The leading NC-by-NC part of this array must contain */
/*             the state dynamics matrix Ac of the controller in a */
/*             block diagonal real Schur form Ac = diag(Ac1,Ac2), where */
/*             Ac1 is (NC-NCS)-by-(NC-NCS) and contains the unstable */
/*             eigenvalues of Ac, and Ac2 is NCS-by-NCS and contains */
/*             the stable eigenvalues of Ac. */

/*     LDAC    INTEGER */
/*             The leading dimension of array AC.  LDAC >= MAX(1,NC). */

/*     BC      (input) DOUBLE PRECISION array, dimension (LDBC,P) */
/*             The leading NC-by-P part of this array must contain */
/*             the input/state matrix Bc of the controller. */

/*     LDBC    INTEGER */
/*             The leading dimension of array BC.  LDBC >= MAX(1,NC). */

/*     CC      (input) DOUBLE PRECISION array, dimension (LDCC,NC) */
/*             The leading M-by-NC part of this array must contain */
/*             the state/output matrix Cc of the controller. */

/*     LDCC    INTEGER */
/*             The leading dimension of array CC.  LDCC >= MAX(1,M). */

/*     DC      (input) DOUBLE PRECISION array, dimension (LDDC,P) */
/*             The leading M-by-P part of this array must contain */
/*             the input/output matrix Dc of the controller. */

/*     LDDC    INTEGER */
/*             The leading dimension of array DC.  LDDC >= MAX(1,M). */

/*     SCALEC  (output) DOUBLE PRECISION */
/*             Scaling factor for the controllability Grammian. */
/*             See METHOD. */

/*     SCALEO  (output) DOUBLE PRECISION */
/*             Scaling factor for the observability Grammian. See METHOD. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,NCS) */
/*             The leading NCS-by-NCS upper triangular part of this array */
/*             contains the Cholesky factor S of the frequency-weighted */
/*             controllability Grammian P = S*S'. See METHOD. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,NCS). */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,NCS) */
/*             The leading NCS-by-NCS upper triangular part of this array */
/*             contains the Cholesky factor R of the frequency-weighted */
/*             observability Grammian Q = R'*R. See METHOD. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,NCS). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(LIWRK) */
/*             LIWRK = 0,       if WEIGHT = 'N'; */
/*             LIWRK = 2(M+P),  if WEIGHT = 'O', 'I', or 'P'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, LFREQ ), */
/*             where */
/*             LFREQ = (N+NC)*(N+NC+2*M+2*P)+ */
/*                     MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7), (M+P)*(M+P+4)) */
/*                                      if WEIGHT = 'I' or 'O' or 'P'; */
/*             LFREQ  = NCS*(MAX(M,P)+5) if WEIGHT = 'N'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the closed-loop system is not well-posed; */
/*                   its feedthrough matrix is (numerically) singular; */
/*             = 2:  the computation of the real Schur form of the */
/*                   closed-loop state matrix failed; */
/*             = 3:  the closed-loop state matrix is not stable; */
/*             = 4:  the solution of a symmetric eigenproblem failed; */
/*             = 5:  the NCS-by-NCS trailing part Ac2 of the state */
/*                   matrix Ac is not stable or not in a real Schur form. */

/*     METHOD */

/*     If JOBC = 'S', the controllability Grammian P is determined as */
/*     follows: */

/*     - if WEIGHT = 'O' or 'N', P satisfies for a continuous-time */
/*       controller the Lyapunov equation */

/*            Ac2*P + P*Ac2' +  scalec^2*Bc*Bc' = 0 */

/*       and for a discrete-time controller */

/*            Ac2*P*Ac2' - P +  scalec^2*Bc*Bc' = 0; */

/*     - if WEIGHT = 'I' or 'P', let Pi be the solution of the */
/*       continuous-time Lyapunov equation */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0 */

/*       or of the discrete-time Lyapunov equation */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0, */

/*       where Ai and Bi are the state and input matrices of a special */
/*       state-space realization of the input frequency weight (see [2]); */
/*       P results as the trailing NCS-by-NCS part of Pi partitioned as */

/*           Pi = ( *  * ). */
/*                ( *  P ) */

/*     If JOBC = 'E', a modified controllability Grammian P1 >= P is */
/*     determined to guarantee stability for a modified Enns' method [2]. */

/*     If JOBO = 'S', the observability Grammian Q is determined as */
/*     follows: */

/*     - if WEIGHT = 'I' or 'N', Q satisfies for a continuous-time */
/*       controller the Lyapunov equation */

/*            Ac2'*Q + Q*Ac2 +  scaleo^2*Cc'*Cc = 0 */

/*       and for a discrete-time controller */

/*            Ac2'*Q*Ac2 - Q +  scaleo^2*Cc'*Cc = 0; */

/*     - if WEIGHT = 'O' or 'P', let Qo be the solution of the */
/*       continuous-time Lyapunov equation */

/*            Ao'*Qo + Qo*Ao +  scaleo^2*Co'*Co = 0 */

/*       or of the discrete-time Lyapunov equation */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0, */

/*       where Ao and Co are the state and output matrices of a */
/*       special state-space realization of the output frequency weight */
/*       (see [2]); if WEIGHT = 'O', Q results as the leading NCS-by-NCS */
/*       part of Qo partitioned as */

/*           Qo = ( Q  * ) */
/*                ( *  * ) */

/*       while if WEIGHT = 'P', Q results as the trailing NCS-by-NCS */
/*       part of Qo partitioned as */

/*           Qo = ( *  * ). */
/*                ( *  Q ) */

/*     If JOBO = 'E', a modified observability Grammian Q1 >= Q is */
/*     determined to guarantee stability for a modified Enns' method [2]. */

/*     The routine computes directly the Cholesky factors S and R */
/*     such that P = S*S' and Q = R'*R according to formulas */
/*     developed in [2]. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. CDC, Las Vegas, pp. 127-132, 1984. */

/*     [2] Varga, A. and Anderson, B.D.O. */
/*         Frequency-weighted balancing related controller reduction. */
/*         Proceedings of the 15th IFAC World Congress, July 21-26, 2002, */
/*         Barcelona, Spain, Vol.15, Part 1, 2002-07-21. */

/*     CONTRIBUTORS */

/*     A. Varga, Australian National University, Canberra, November 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*     May 2009. */
/*     A. Varga, DLR Oberpfafenhofen, June 2001. */


/*     KEYWORDS */

/*     Controller reduction, frequency weighting, multivariable system, */
/*     state-space model, state-space representation. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    ac_dim1 = *ldac;
    ac_offset = 1 + ac_dim1;
    ac -= ac_offset;
    bc_dim1 = *ldbc;
    bc_offset = 1 + bc_dim1;
    bc -= bc_offset;
    cc_dim1 = *ldcc;
    cc_offset = 1 + cc_dim1;
    cc -= cc_offset;
    dc_dim1 = *lddc;
    dc_offset = 1 + dc_dim1;
    dc -= dc_offset;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --iwork;
    --dwork;

    /* Function Body */
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    leftw = lsame_(weight, "O", (ftnlen)1, (ftnlen)1);
    rightw = lsame_(weight, "I", (ftnlen)1, (ftnlen)1);
    perf = lsame_(weight, "P", (ftnlen)1, (ftnlen)1);
    frwght = leftw || rightw || perf;

    *info = 0;
    nnc = *n + *nc;
    mp = *m + *p;
    if (frwght) {
/* Computing MAX */
/* Computing MAX */
	i__3 = max(nnc,*m);
	i__1 = nnc * (nnc + max(i__3,*p) + 7), i__2 = mp * (mp + 4);
	lw = nnc * (nnc + (mp << 1)) + max(i__1,i__2);
    } else {
	lw = *ncs * (max(*m,*p) + 5);
    }
    lw = max(1,lw);

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (lsame_(jobc, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobc, 
	    "E", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lsame_(jobo, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobo, 
	    "E", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*p < 0) {
	*info = -7;
    } else if (*nc < 0) {
	*info = -8;
    } else if (*ncs < 0 || *ncs > *nc) {
	*info = -9;
    } else if (*lda < max(1,*n)) {
	*info = -11;
    } else if (*ldb < max(1,*n)) {
	*info = -13;
    } else if (*ldc < max(1,*p)) {
	*info = -15;
    } else if (*ldd < max(1,*p)) {
	*info = -17;
    } else if (*ldac < max(1,*nc)) {
	*info = -19;
    } else if (*ldbc < max(1,*nc)) {
	*info = -21;
    } else if (*ldcc < max(1,*m)) {
	*info = -23;
    } else if (*lddc < max(1,*m)) {
	*info = -25;
    } else if (*lds < max(1,*ncs)) {
	*info = -29;
    } else if (*ldr < max(1,*ncs)) {
	*info = -31;
    } else if (*ldwork < lw) {
	*info = -34;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB16AY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *scalec = 1.;
    *scaleo = 1.;
/* Computing MIN */
    i__1 = min(*ncs,*m);
    if (min(i__1,*p) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    wrkopt = 1;
    ncu = *nc - *ncs;
    ncu1 = ncu + 1;

    if (! perf) {

/*        Compute the Grammians in the case of no weighting or */
/*        one-sided weighting. */

	if (leftw || lsame_(weight, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute the standard controllability Grammian. */

/*           Solve for the Cholesky factor S of P, P = S*S', */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */

/*               Ac2*P + P*Ac2' +  scalec^2*Bc2*Bc2' = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */

/*               Ac2*P*Ac2' - P +  scalec^2*Bc2*Bc2' = 0, */

/*           where Bc2 is the matrix formed from the last NCS rows of Bc. */

/*           Workspace:  need   NCS*(P+5); */
/*                              prefer larger. */
	    ku = 1;
	    ktau = ku + *ncs * *p;
	    kw = ktau + *ncs;

	    dlacpy_("Full", ncs, p, &bc[ncu1 + bc_dim1], ldbc, &dwork[ku], 
		    ncs, (ftnlen)4);
	    i__1 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_true, ncs, p, &ac[ncu1 + ncu1 * ac_dim1], ldac,
		     &dwork[ku], ncs, &dwork[ktau], &s[s_offset], lds, scalec,
		     &dwork[kw], &i__1, &ierr);
	    if (ierr != 0) {
		*info = 5;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	}

	if (rightw || lsame_(weight, "N", (ftnlen)1, (ftnlen)1)) {

/*           Compute the standard observability Grammian. */

/*           Solve for the Cholesky factor R of Q, Q = R'*R, */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */

/*               Ac2'*Q + Q*Ac2  +  scaleo^2*Cc2'*Cc2 = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */

/*               Ac2'*Q*Ac2 - Q +  scaleo^2*Cc2'*Cc2 = 0, */

/*           where Cc2 is the matrix formed from the last NCS columns */
/*           of Cc. */

/*           Workspace:  need   NCS*(M + 5); */
/*                              prefer larger. */
	    ku = 1;
	    ktau = ku + *m * *ncs;
	    kw = ktau + *ncs;

	    dlacpy_("Full", m, ncs, &cc[ncu1 * cc_dim1 + 1], ldcc, &dwork[ku],
		     m, (ftnlen)4);
	    i__1 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_false, ncs, m, &ac[ncu1 + ncu1 * ac_dim1], 
		    ldac, &dwork[ku], m, &dwork[ktau], &r__[r_offset], ldr, 
		    scaleo, &dwork[kw], &i__1, &ierr);
	    if (ierr != 0) {
		*info = 5;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	}

/*        Finish if there are no weights. */

	if (lsame_(weight, "N", (ftnlen)1, (ftnlen)1)) {
	    dwork[1] = (doublereal) wrkopt;
	    return 0;
	}
    }

    if (frwght) {

/*        Allocate working storage for computing the weights. */

/*        Real workspace:    need MAX(1,NNC*NNC+2*NNC*MP+MP*(MP+4)); */
/*        Integer workspace: need 2*MP. */

	kwa = 1;
	kwb = kwa + nnc * nnc;
	kwc = kwb + nnc * mp;
	kwd = kwc + nnc * mp;
	kw = kwd + mp * mp;
	kl = kwd;

	if (leftw) {

/*           Build the extended matrices */

/*           Ao = ( Ac+Bc*inv(R)*D*Cc   Bc*inv(R)*C   ), */
/*                (     B*inv(Rt)*Cc  A+B*Dc*inv(R)*C ) */

/*           Co = ( -inv(R)*D*Cc  -inv(R)*C ) , */

/*           where  R = I-D*Dc and Rt = I-Dc*D. */
/*                             -1 */
/*           Method: Compute Ge  = ( Ge11 Ge12 ), where Ge = ( K   -Im ). */
/*                                 ( Ge21 Ge22 )             ( -Ip  G  ) */

/*                               -1 */
/*           Then  Ge11 = -(I-G*K) *G . */

/*           Construct first Ge = (  K  -Im ) such that the stable part */
/*                                ( -Ip  G  ) */
/*           of K is in the leading position (to avoid updating of */
/*           QR factorization). */

	    dlaset_("Full", m, p, &c_b21, &c_b21, &dwork[kwd], &mp, (ftnlen)4)
		    ;
	    ab05pd_("N", ncs, p, m, &ncu, &c_b24, &ac[ncu1 + ncu1 * ac_dim1], 
		    ldac, &bc[ncu1 + bc_dim1], ldbc, &cc[ncu1 * cc_dim1 + 1], 
		    ldcc, &dwork[kwd], &mp, &ac[ac_offset], ldac, &bc[
		    bc_offset], ldbc, &cc[cc_offset], ldcc, &dc[dc_offset], 
		    lddc, &ne, &dwork[kwa], &nnc, &dwork[kwb], &nnc, &dwork[
		    kwc], &mp, &dwork[kwd], &mp, &ierr, (ftnlen)1);
	    ab05qd_("Over", nc, p, m, n, m, p, &dwork[kwa], &nnc, &dwork[kwb],
		     &nnc, &dwork[kwc], &mp, &dwork[kwd], &mp, &a[a_offset], 
		    lda, &b[b_offset], ldb, &c__[c_offset], ldc, &d__[
		    d_offset], ldd, &ne, &me, &pe, &dwork[kwa], &nnc, &dwork[
		    kwb], &nnc, &dwork[kwc], &mp, &dwork[kwd], &mp, &ierr, (
		    ftnlen)4);
	    dlaset_("Full", m, m, &c_b21, &c_b28, &dwork[kwd + mp * *p], &mp, 
		    (ftnlen)4);
	    dlaset_("Full", p, p, &c_b21, &c_b28, &dwork[kwd + *m], &mp, (
		    ftnlen)4);

	} else {

/*           Build the extended matrices */

/*           Ai = ( A+B*Dc*inv(R)*C   B*inv(Rt)*Cc   ) , */
/*                (   Bc*inv(R)*C  Ac+Bc*inv(R)*D*Cc ) */

/*           Bi = ( B*Dc*inv(R)    B*inv(Rt)  ) , */
/*                ( Bc*inv(R)    Bc*D*inv(Rt) ) */

/*           Ci = (  -inv(R)*C   -inv(R)*D*Cc ) , where */

/*           R = I-D*Dc and Rt = I-Dc*D. */

/*                             -1 */
/*           Method: Compute Ge  = ( Ge11 Ge12 ), where Ge = ( G   -Ip ). */
/*                                 ( Ge21 Ge22 )             ( -Im  K  ) */

/*                              -1                     -1 */
/*           Then Ge22 = -(I-G*K) *G and Ge21 = -(I-G*K) . */

/*           Construct first Ge = (  G  -Ip ). */
/*                                ( -Im  K  ) */

	    ab05qd_("N", n, m, p, nc, p, m, &a[a_offset], lda, &b[b_offset], 
		    ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &ac[
		    ac_offset], ldac, &bc[bc_offset], ldbc, &cc[cc_offset], 
		    ldcc, &dc[dc_offset], lddc, &ne, &me, &pe, &dwork[kwa], &
		    nnc, &dwork[kwb], &nnc, &dwork[kwc], &mp, &dwork[kwd], &
		    mp, &ierr, (ftnlen)1);
	    dlaset_("Full", p, p, &c_b21, &c_b28, &dwork[kwd + mp * *m], &mp, 
		    (ftnlen)4);
	    dlaset_("Full", m, m, &c_b21, &c_b28, &dwork[kwd + *p], &mp, (
		    ftnlen)4);
	}
/*                  -1 */
/*        Compute Ge   = ( Ge11 Ge12 ). */
/*                       ( Ge21 Ge22 ) */

/*        Additional real workspace: need 4*MP; */
/*        Integer workspace:         need 2*MP. */

	i__1 = *ldwork - kw + 1;
	ab07nd_(&nnc, &mp, &dwork[kwa], &nnc, &dwork[kwb], &nnc, &dwork[kwc], 
		&mp, &dwork[kwd], &mp, &rcond, &iwork[1], &dwork[kw], &i__1, &
		ierr);
	if (ierr != 0) {
	    *info = 1;
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);

/*                     -1   ( A1 | B1  B2  ) */
/*        Partition  Ge   = (--------------) and select appropriate */
/*                          ( C1 | D11 D12 ) */
/*                          ( C2 | D21 D22 ) */

/*        pointers to matrices and column dimensions to define weights. */

	if (rightw) {

/*           Define B2 for Ge22. */

	    me = *m;
	    kwb += nnc * *p;
	} else if (perf) {

/*           Define B1 and C2 for Ge21. */

	    me = *p;
	    kwc += *m;
	}
    }

    if (leftw || perf) {

/*        Compute the frequency-weighted observability Grammian. */

/*        Solve for the Cholesky factor Ro of Qo, Qo = Ro'*Ro, */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ao'*Qo + Qo*Ao  +  scaleo^2*Co'*Co = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0. */

/*        Additional workspace:  need   NNC*(NNC+MAX(NNC,P)+7); */
/*                               prefer larger. */

	ldu = max(nnc,*p);
	ku = kl;
	kq = ku + nnc * ldu;
	kr = kq + nnc * nnc;
	ki = kr + nnc;
	kw = ki + nnc;

	*(unsigned char *)jobfac = 'N';
	dlacpy_("Full", p, &nnc, &dwork[kwc], &mp, &dwork[ku], &ldu, (ftnlen)
		4);
	i__1 = *ldwork - kw + 1;
	sb03od_(dico, jobfac, "No-transpose", &nnc, p, &dwork[kwa], &nnc, &
		dwork[kq], &nnc, &dwork[ku], &ldu, scaleo, &dwork[kr], &dwork[
		ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
		12);
	if (ierr != 0) {
	    if (ierr == 6) {
		*info = 2;
	    } else {
		*info = 3;
	    }
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);

/*        Partition Ro as Ro = ( R11 R12 ). */
/*                             (  0  R22 ) */

	if (leftw) {

/*           R = R11 (NCS-by-NCS). */

	    dlacpy_("Upper", ncs, ncs, &dwork[ku], &ldu, &r__[r_offset], ldr, 
		    (ftnlen)5);
	} else {

/*           Compute R such that R'*R = R22'*R22 + R12'*R12, where */
/*           R22 is NCS-by-NCS and R12 is (N+NCU)-by-NCS. */
/*           R22 corresponds to the stable part of the controller. */

	    nncu = *n + ncu;
	    dlacpy_("Upper", ncs, ncs, &dwork[ku + (ldu + 1) * nncu], &ldu, &
		    r__[r_offset], ldr, (ftnlen)5);
	    ktau = ku;
	    mb04od_("Full", ncs, &c__0, &nncu, &r__[r_offset], ldr, &dwork[ku 
		    + ldu * nncu], &ldu, dum, &c__1, dum, &c__1, &dwork[ktau],
		     &dwork[kw], (ftnlen)4);

	    i__1 = *ncs;
	    for (j = 1; j <= i__1; ++j) {
		if (r__[j + j * r_dim1] < 0.) {
		    i__2 = *ncs - j + 1;
		    dscal_(&i__2, &c_b28, &r__[j + j * r_dim1], ldr);
		}
/* L10: */
	    }
	}
    }

    if (rightw || perf) {

/*        Compute the frequency-weighted controllability Grammian. */

/*        Solve for the Cholesky factor Si of Pi, Pi = Si*Si', */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0. */

/*        Additional workspace:  need   NNC*(NNC+MAX(NNC,P,M)+7); */
/*                               prefer larger. */

	ku = kl;
	kq = ku + nnc * max(nnc,me);
	kr = kq + nnc * nnc;
	ki = kr + nnc;
	kw = ki + nnc;

	dlacpy_("Full", &nnc, &me, &dwork[kwb], &nnc, &dwork[ku], &nnc, (
		ftnlen)4);
	*(unsigned char *)jobfac = 'F';
	if (rightw) {
	    *(unsigned char *)jobfac = 'N';
	}
	i__1 = *ldwork - kw + 1;
	sb03od_(dico, jobfac, "Transpose", &nnc, &me, &dwork[kwa], &nnc, &
		dwork[kq], &nnc, &dwork[ku], &nnc, scalec, &dwork[kr], &dwork[
		ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
		9);
	if (ierr != 0) {
	    if (ierr == 6) {
		*info = 2;
	    } else {
		*info = 3;
	    }
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);

/*        Partition Si as Si = ( S11 S12 ) with S22 NCS-by-NCS and */
/*                             (  0  S22 ) */
/*        set S = S22. */

	nncu = *n + ncu;
	dlacpy_("Upper", ncs, ncs, &dwork[ku + (nnc + 1) * nncu], &nnc, &s[
		s_offset], lds, (ftnlen)5);
    }

    ku = 1;
    if (leftw || perf) {
	if (lsame_(jobo, "E", (ftnlen)1, (ftnlen)1)) {

/*           Form Y = -Ac2'*(R'*R)-(R'*R)*Ac2 if DICO = 'C', or */
/*                Y = -Ac2'*(R'*R)*Ac2+(R'*R) if DICO = 'D'. */

/*           Workspace:  need   2*NCS*NCS. */

	    dlacpy_("Upper", ncs, ncs, &r__[r_offset], ldr, &dwork[ku], ncs, (
		    ftnlen)5);
	    dlacpy_("Full", ncs, ncs, &ac[ncu1 + ncu1 * ac_dim1], ldac, &
		    dwork[ku + *ncs * *ncs], ncs, (ftnlen)4);
	    mb01wd_(dico, "Upper", "No-transpose", "Hessenberg", ncs, &c_b28, 
		    &c_b21, &r__[r_offset], ldr, &dwork[ku + *ncs * *ncs], 
		    ncs, &dwork[ku], ncs, &ierr, (ftnlen)1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)10);

/*           Compute the eigendecomposition of Y as Y = Z*Sigma*Z'. */

	    kw = ku + *ncs;
	    i__1 = *ldwork - kw + 1;
	    dsyev_("Vectors", "Upper", ncs, &r__[r_offset], ldr, &dwork[ku], &
		    dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)5);
	    if (ierr > 0) {
		*info = 4;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 <= 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
	    d__3 = (d__1 = dwork[ku], abs(d__1)), d__4 = (d__2 = dwork[ku + *
		    ncs - 1], abs(d__2));
	    tol = max(d__3,d__4) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form Cc = [ sqrt(Sigma2)*Z2' ] */

	    pcbar = 0;
	    jj = ku;
	    i__1 = *ncs;
	    for (j = 1; j <= i__1; ++j) {
		if (dwork[jj] > tol) {
		    d__1 = sqrt(dwork[jj]);
		    dscal_(ncs, &d__1, &r__[j * r_dim1 + 1], &c__1);
		    dcopy_(ncs, &r__[j * r_dim1 + 1], &c__1, &dwork[kw + 
			    pcbar], ncs);
		    ++pcbar;
		}
		++jj;
/* L20: */
	    }

/*           Solve for the Cholesky factor R of Q, Q = R'*R, */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                               _   _ */
/*                   Ac2'*Q + Q*Ac2  +  scaleo^2*Cc'*Cc = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                              _   _ */
/*                   Ac2'*Q*Ac2 - Q +  scaleo^2*Cc'*Cc = 0. */

/*           Workspace:  need   NCS*(NCS + 6); */
/*                              prefer larger. */

	    ku = kw;
	    ktau = ku + *ncs * *ncs;
	    kw = ktau + *ncs;

	    i__1 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_false, ncs, &pcbar, &ac[ncu1 + ncu1 * ac_dim1],
		     ldac, &dwork[ku], ncs, &dwork[ktau], &r__[r_offset], ldr,
		     &t, &dwork[kw], &i__1, &ierr);
	    if (ierr != 0) {
		*info = 5;
		return 0;
	    }
	    *scaleo *= t;
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	}

    }

    if (rightw || perf) {
	if (lsame_(jobc, "E", (ftnlen)1, (ftnlen)1)) {

/*           Form X = -A2c*(S*S')-(S*S')*Ac2' if DICO = 'C', or */
/*                X = -Ac2*(S*S')*Ac2'+(S*S') if DICO = 'D'. */

/*           Workspace:  need   2*NCS*NCS. */

	    dlacpy_("Upper", ncs, ncs, &s[s_offset], lds, &dwork[ku], ncs, (
		    ftnlen)5);
	    dlacpy_("Full", ncs, ncs, &ac[ncu1 + ncu1 * ac_dim1], ldac, &
		    dwork[ku + *ncs * *ncs], ncs, (ftnlen)4);
	    mb01wd_(dico, "Upper", "Transpose", "Hessenberg", ncs, &c_b28, &
		    c_b21, &s[s_offset], lds, &dwork[ku + *ncs * *ncs], ncs, &
		    dwork[ku], ncs, &ierr, (ftnlen)1, (ftnlen)5, (ftnlen)9, (
		    ftnlen)10);

/*           Compute the eigendecomposition of X as X = Z*Sigma*Z'. */

	    kw = ku + *ncs;
	    i__1 = *ldwork - kw + 1;
	    dsyev_("Vectors", "Upper", ncs, &s[s_offset], lds, &dwork[ku], &
		    dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)5);
	    if (ierr > 0) {
		*info = 4;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 =< 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
	    d__3 = (d__1 = dwork[ku], abs(d__1)), d__4 = (d__2 = dwork[ku + *
		    ncs - 1], abs(d__2));
	    tol = max(d__3,d__4) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form Bc = [ Z2*sqrt(Sigma2) ] */

	    mbbar = 0;
	    i__ = kw;
	    jj = ku;
	    i__1 = *ncs;
	    for (j = 1; j <= i__1; ++j) {
		if (dwork[jj] > tol) {
		    ++mbbar;
		    d__1 = sqrt(dwork[jj]);
		    dscal_(ncs, &d__1, &s[j * s_dim1 + 1], &c__1);
		    dcopy_(ncs, &s[j * s_dim1 + 1], &c__1, &dwork[i__], &c__1)
			    ;
		    i__ += *ncs;
		}
		++jj;
/* L30: */
	    }

/*           Solve for the Cholesky factor S of P, P = S*S', */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                               _  _ */
/*                   Ac2*P + P*Ac2'  +  scalec^2*Bc*Bc' = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                              _  _ */
/*                   Ac2*P*Ac2' - P +  scalec^2*Bc*Bc' = 0. */

/*           Workspace:  need   maximum NCS*(NCS + 6); */
/*                       prefer larger. */

	    ku = kw;
	    ktau = ku + mbbar * *ncs;
	    kw = ktau + *ncs;

	    i__1 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_true, ncs, &mbbar, &ac[ncu1 + ncu1 * ac_dim1], 
		    ldac, &dwork[ku], ncs, &dwork[ktau], &s[s_offset], lds, &
		    t, &dwork[kw], &i__1, &ierr);
	    if (ierr != 0) {
		*info = 5;
		return 0;
	    }
	    *scalec *= t;
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	}

    }

/*     Save optimal workspace. */

    dwork[1] = (doublereal) wrkopt;

    return 0;
/* *** Last line of SB16AY *** */
} /* sb16ay_ */

