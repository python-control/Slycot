/* AB09KX.f -- translated by f2c (version 20100827).
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

static doublereal c_b24 = -1.;
static doublereal c_b25 = 0.;
static integer c_n1 = -1;
static doublereal c_b31 = 1.;
static integer c__1 = 1;

/* Subroutine */ int ab09kx_(char *job, char *dico, char *weight, integer *n, 
	integer *nv, integer *nw, integer *m, integer *p, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *av, integer *ldav, 
	doublereal *bv, integer *ldbv, doublereal *cv, integer *ldcv, 
	doublereal *dv, integer *lddv, doublereal *aw, integer *ldaw, 
	doublereal *bw, integer *ldbw, doublereal *cw, integer *ldcw, 
	doublereal *dw, integer *lddw, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen job_len, ftnlen dico_len, 
	ftnlen weight_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, av_dim1, av_offset, bv_dim1, bv_offset, cv_dim1, 
	    cv_offset, dv_dim1, dv_offset, aw_dim1, aw_offset, bw_dim1, 
	    bw_offset, cw_dim1, cw_offset, dw_dim1, dw_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, ia, ib, kw, lw, ldw, ierr, ldwn;
    static doublereal work, scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01wd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static logical discr, conjs, leftw;
    extern /* Subroutine */ int sb04py_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical frwght, rightw;
    extern /* Subroutine */ int dtrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen);


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

/*     To construct a state-space representation (A,BS,CS,DS) of the */
/*     stable projection of V*G*W or conj(V)*G*conj(W) from the */
/*     state-space representations (A,B,C,D), (AV,BV,CV,DV), and */
/*     (AW,BW,CW,DW) of the transfer-function matrices G, V and W, */
/*     respectively. G is assumed to be a stable transfer-function */
/*     matrix and the state matrix A must be in a real Schur form. */
/*     When computing the stable projection of V*G*W, V and W are assumed */
/*     to be completely unstable transfer-function matrices. */
/*     When computing the stable projection of conj(V)*G*conj(W), */
/*     V and W are assumed to be stable transfer-function matrices. */

/*     For a transfer-function matrix G, conj(G) denotes the conjugate */
/*     of G given by G'(-s) for a continuous-time system or G'(1/z) */
/*     for a discrete-time system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies which projection to be computed as follows: */
/*             = 'N':  compute the stable projection of V*G*W; */
/*             = 'C':  compute the stable projection of */
/*                     conj(V)*G*conj(W). */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G, V and W are continuous-time systems; */
/*             = 'D':  G, V and W are discrete-time systems. */

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'L':  only left weighting V is used (W = I); */
/*             = 'R':  only right weighting W is used (V = I); */
/*             = 'B':  both left and right weightings V and W are used. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. Also the number of rows of */
/*             the matrix B and the number of columns of the matrix C. */
/*             N represents the dimension of the state vector of the */
/*             system with the transfer-function matrix G.  N >= 0. */

/*     NV      (input) INTEGER */
/*             The order of the matrix AV. Also the number of rows of */
/*             the matrix BV and the number of columns of the matrix CV. */
/*             NV represents the dimension of the state vector of the */
/*             system with the transfer-function matrix V.  NV >= 0. */

/*     NW      (input) INTEGER */
/*             The order of the matrix AW. Also the number of rows of */
/*             the matrix BW and the number of columns of the matrix CW. */
/*             NW represents the dimension of the state vector of the */
/*             system with the transfer-function matrix W.  NW >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices B, D, BW and DW */
/*             and number of rows of the matrices CW and DW.  M >= 0. */
/*             M represents the dimension of input vectors of the */
/*             systems with the transfer-function matrices G and W and */
/*             also the dimension of the output vector of the system */
/*             with the transfer-function matrix W. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrices C, D, CV and DV and the */
/*             number of columns of the matrices BV and DV.  P >= 0. */
/*             P represents the dimension of output vectors of the */
/*             systems with the transfer-function matrices G and V and */
/*             also the dimension of the input vector of the system */
/*             with the transfer-function matrix V. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must */
/*             contain the state matrix A of the system with the */
/*             transfer-function matrix G in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading N-by-M part of this */
/*             array contains the input matrix BS of the stable */
/*             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W) */
/*             if JOB = 'C'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading P-by-N part of this */
/*             array contains the output matrix CS of the stable */
/*             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W) */
/*             if JOB = 'C'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the feedthrough matrix D of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading P-by-M part of this */
/*             array contains the feedthrough matrix DS of the stable */
/*             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W) */
/*             if JOB = 'C'. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,P). */

/*     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV */
/*             part of this array must contain the state matrix AV of */
/*             the system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             NV-by-NV part of this array contains a real Schur form */
/*             of AV. */
/*             AV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDAV    INTEGER */
/*             The leading dimension of the array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part */
/*             of this array must contain the input matrix BV of the */
/*             system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             NV-by-P part of this array contains the transformed input */
/*             matrix BV. */
/*             BV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDBV    INTEGER */
/*             The leading dimension of the array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading P-by-NV part */
/*             of this array must contain the output matrix CV of the */
/*             system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             P-by-NV part of this array contains the transformed output */
/*             matrix CV. */
/*             CV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDCV    INTEGER */
/*             The leading dimension of the array CV. */
/*             LDCV >= MAX(1,P), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,        if WEIGHT = 'R' or 'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If WEIGHT = 'L' or 'B', the leading P-by-P part of this */
/*             array must contain the feedthrough matrix DV of the system */
/*             with the transfer-function matrix V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDDV    INTEGER */
/*             The leading dimension of the array DV. */
/*             LDDV >= MAX(1,P), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,        if WEIGHT = 'R' or 'N'. */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW */
/*             part of this array must contain the state matrix AW of */
/*             the system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             NW-by-NW part of this array contains a real Schur form */
/*             of AW. */
/*             AW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDAW    INTEGER */
/*             The leading dimension of the array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-M part */
/*             of this array must contain the input matrix BW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             NW-by-M part of this array contains the transformed input */
/*             matrix BW. */
/*             BW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDBW    INTEGER */
/*             The leading dimension of the array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part */
/*             of this array must contain the output matrix CW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             M-by-NW part of this array contains the transformed output */
/*             matrix CW. */
/*             CW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDCW    INTEGER */
/*             The leading dimension of the array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,M) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-M part of this */
/*             array must contain the feedthrough matrix DW of the system */
/*             with the transfer-function matrix W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDDW    INTEGER */
/*             The leading dimension of the array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, LDW1, LDW2 ), where */
/*               LDW1 = 0 if WEIGHT = 'R' or 'N' and */
/*               LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                      if WEIGHT = 'L' or WEIGHT = 'B', */
/*               LDW2 = 0 if WEIGHT = 'L' or 'N' and */
/*               LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) ) */
/*                      if WEIGHT = 'R' or WEIGHT = 'B', */
/*               a = 0,    b = 0,     if DICO = 'C' or  JOB = 'N', */
/*               a = 2*NV, b = 2*NW,  if DICO = 'D' and JOB = 'C'. */
/*             For good performance, LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             =  0:  no warning; */
/*             =  1:  JOB = 'N' and AV is not completely unstable, or */
/*                    JOB = 'C' and AV is not stable; */
/*             =  2:  JOB = 'N' and AW is not completely unstable, or */
/*                    JOB = 'C' and AW is not stable; */
/*             =  3:  both above conditions appear. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  the reduction of AV to a real Schur form failed; */
/*             =  2:  the reduction of AW to a real Schur form failed; */
/*             =  3:  the solution of the Sylvester equation failed */
/*                    because the matrices A and AV have common */
/*                    eigenvalues (if JOB = 'N'), or -AV and A have */
/*                    common eigenvalues (if JOB = 'C' and DICO = 'C'), */
/*                    or AV has an eigenvalue which is the reciprocal of */
/*                    one of the eigenvalues of A (if JOB = 'C' and */
/*                    DICO = 'D'); */
/*             =  4:  the solution of the Sylvester equation failed */
/*                    because the matrices A and AW have common */
/*                    eigenvalues (if JOB = 'N'), or -AW and A have */
/*                    common eigenvalues (if JOB = 'C' and DICO = 'C'), */
/*                    or AW has an eigenvalue which is the reciprocal of */
/*                    one of the eigenvalues of A (if JOB = 'C' and */
/*                    DICO = 'D'). */

/*     METHOD */

/*     The matrices of the stable projection of V*G*W are computed as */

/*       BS = B*DW + Y*BW,  CS = CV*X + DV*C,  DS = DV*D*DW, */

/*     where X and Y satisfy the continuous-time Sylvester equations */

/*       AV*X - X*A  + BV*C = 0, */
/*       -A*Y + Y*AW + B*CW = 0. */

/*     The matrices of the stable projection of conj(V)*G*conj(W) are */
/*     computed using the explicit formulas established in [1]. */

/*     For a continuous-time system, the matrices BS, CS and DS of */
/*     the stable projection are computed as */

/*       BS = B*DW' + Y*CW',  CS = BV'*X + DV'*C,  DS = DV'*D*DW', */

/*     where X and Y satisfy the continuous-time Sylvester equations */

/*       AV'*X + X*A   + CV'*C = 0, */
/*         A*Y + Y*AW' + B*BW' = 0. */

/*     For a discrete-time system, the matrices BS, CS and DS of */
/*     the stable projection are computed as */

/*       BS = B*DW' + A*Y*CW',  CS = BV'*X*A + DV'*C, */
/*       DS = DV'*D*DW' + BV'*X*B*DW' + DV'*C*Y*CW' + BV'*X*A*Y*CW', */

/*     where X and Y satisfy the discrete-time Sylvester equations */

/*       AV'*X*A + CV'*C = X, */
/*       A*Y*AW' + B*BW' = Y. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         Explicit formulas for an efficient implementation */
/*         of the frequency-weighting model reduction approach. */
/*         Proc. 1993 European Control Conference, Groningen, NL, */
/*         pp. 693-696, 1993. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on numerically stable algorithms. */

/*     FURTHER COMMENTS */

/*     The matrix A must be stable, but its stability is not checked by */
/*     this routine. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, April 2000. */
/*     D. Sima, University of Bucharest, May 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2000. */
/*     Based on the RASP routines SFRLW, SFRLW1, SFRRW and SFRRW1, */
/*     by A. Varga, 1992. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Frequency weighting, model reduction, multivariable system, */
/*     state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
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
    av_dim1 = *ldav;
    av_offset = 1 + av_dim1;
    av -= av_offset;
    bv_dim1 = *ldbv;
    bv_offset = 1 + bv_dim1;
    bv -= bv_offset;
    cv_dim1 = *ldcv;
    cv_offset = 1 + cv_dim1;
    cv -= cv_offset;
    dv_dim1 = *lddv;
    dv_offset = 1 + dv_dim1;
    dv -= dv_offset;
    aw_dim1 = *ldaw;
    aw_offset = 1 + aw_dim1;
    aw -= aw_offset;
    bw_dim1 = *ldbw;
    bw_offset = 1 + bw_dim1;
    bw -= bw_offset;
    cw_dim1 = *ldcw;
    cw_offset = 1 + cw_dim1;
    cw -= cw_offset;
    dw_dim1 = *lddw;
    dw_offset = 1 + dw_dim1;
    dw -= dw_offset;
    --dwork;

    /* Function Body */
    conjs = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
    frwght = leftw || rightw;

    *iwarn = 0;
    *info = 0;
    if (discr && conjs) {
	ia = *nv << 1;
	ib = *nw << 1;
    } else {
	ia = 0;
	ib = 0;
    }
    lw = 1;
    if (leftw) {
/* Computing MAX */
/* Computing MAX */
	i__3 = ia, i__4 = *p * *n, i__3 = max(i__3,i__4), i__4 = *p * *m;
	i__1 = lw, i__2 = *nv * (*nv + 5), i__1 = max(i__1,i__2), i__2 = *nv *
		 *n + max(i__3,i__4);
	lw = max(i__1,i__2);
    }
    if (rightw) {
/* Computing MAX */
/* Computing MAX */
	i__3 = ib, i__4 = *m * *n, i__3 = max(i__3,i__4), i__4 = *p * *m;
	i__1 = lw, i__2 = *nw * (*nw + 5), i__1 = max(i__1,i__2), i__2 = *nw *
		 *n + max(i__3,i__4);
	lw = max(i__1,i__2);
    }

/*     Test the input scalar arguments. */

    if (! (lsame_(job, "N", (ftnlen)1, (ftnlen)1) || conjs)) {
	*info = -1;
    } else if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -2;
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*nv < 0) {
	*info = -5;
    } else if (*nw < 0) {
	*info = -6;
    } else if (*m < 0) {
	*info = -7;
    } else if (*p < 0) {
	*info = -8;
    } else if (*lda < max(1,*n)) {
	*info = -10;
    } else if (*ldb < max(1,*n)) {
	*info = -12;
    } else if (*ldc < max(1,*p)) {
	*info = -14;
    } else if (*ldd < max(1,*p)) {
	*info = -16;
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
	*info = -18;
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
	*info = -20;
    } else if (*ldcv < 1 || leftw && *ldcv < *p) {
	*info = -22;
    } else if (*lddv < 1 || leftw && *lddv < *p) {
	*info = -24;
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
	*info = -26;
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
	*info = -28;
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
	*info = -30;
    } else if (*lddw < 1 || rightw && *lddw < *m) {
	*info = -32;
    } else if (*ldwork < lw) {
	*info = -34;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09KX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (! frwght || min(*m,*p) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    work = 1.;
    if (leftw && *nv > 0) {

/*        Reduce AV to a real Schur form using an orthogonal similarity */
/*        transformation AV <- Q'*AV*Q and apply the transformation to */
/*        BV and CV: BV <- Q'*BV and CV <- CV*Q. */

/*        Workspace needed:  NV*(NV+5); */
/*                           prefer larger. */

	kw = *nv * (*nv + 2) + 1;
	i__1 = *ldwork - kw + 1;
	tb01wd_(nv, p, p, &av[av_offset], ldav, &bv[bv_offset], ldbv, &cv[
		cv_offset], ldcv, &dwork[(*nv << 1) + 1], nv, &dwork[1], &
		dwork[*nv + 1], &dwork[kw], &i__1, &ierr);
	if (ierr != 0) {
	    *info = 1;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);

	if (conjs) {

/*           Check the stability of the eigenvalues of AV. */

	    if (discr) {
		i__1 = *nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dlapy2_(&dwork[i__], &dwork[*nv + i__]) >= 1.) {
			*iwarn = 1;
			goto L50;
		    }
/* L10: */
		}
	    } else {
		i__1 = *nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dwork[i__] >= 0.) {
			*iwarn = 1;
			goto L50;
		    }
/* L20: */
		}
	    }
	} else {

/*           Check the anti-stability of the eigenvalues of AV. */

	    if (discr) {
		i__1 = *nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dlapy2_(&dwork[i__], &dwork[*nv + i__]) <= 1.) {
			*iwarn = 1;
			goto L50;
		    }
/* L30: */
		}
	    } else {
		i__1 = *nv;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dwork[i__] <= 0.) {
			*iwarn = 1;
			goto L50;
		    }
/* L40: */
		}
	    }
	}
L50:

	;
    }

    if (rightw && *nw > 0) {

/*        Reduce AW to a real Schur form using an orthogonal similarity */
/*        transformation AW <- T'*AW*T and apply the transformation to */
/*        BW and CW: BW <- T'*BW and CW <- CW*T. */

/*        Workspace needed:  NW*(NW+5); */
/*                           prefer larger. */

	kw = *nw * (*nw + 2) + 1;
	i__1 = *ldwork - kw + 1;
	tb01wd_(nw, m, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw, &cw[
		cw_offset], ldcw, &dwork[(*nw << 1) + 1], nw, &dwork[1], &
		dwork[*nw + 1], &dwork[kw], &i__1, &ierr);
	if (ierr != 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);

	if (conjs) {

/*           Check the stability of the eigenvalues of AW. */

	    if (discr) {
		i__1 = *nw;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dlapy2_(&dwork[i__], &dwork[*nw + i__]) >= 1.) {
			*iwarn += 2;
			goto L100;
		    }
/* L60: */
		}
	    } else {
		i__1 = *nw;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dwork[i__] >= 0.) {
			*iwarn += 2;
			goto L100;
		    }
/* L70: */
		}
	    }
	} else {

/*           Check the anti-stability of the eigenvalues of AW. */

	    if (discr) {
		i__1 = *nw;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dlapy2_(&dwork[i__], &dwork[*nw + i__]) <= 1.) {
			*iwarn += 2;
			goto L100;
		    }
/* L80: */
		}
	    } else {
		i__1 = *nw;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dwork[i__] <= 0.) {
			*iwarn += 2;
			goto L100;
		    }
/* L90: */
		}
	    }
	}
L100:
	;
    }

    if (leftw) {
	ldw = max(*nv,1);
	kw = *nv * *n + 1;
	if (conjs) {

/*           Compute the projection of conj(V)*G. */

/*           Total workspace needed:  NV*N + MAX( a, P*N, P*M ), where */
/*                                    a = 0,    if DICO = 'C', */
/*                                    a = 2*NV, if DICO = 'D'. */

/*           Compute -CV'*C. */
/*           Workspace needed: NV*N. */

	    dgemm_("T", "N", nv, n, p, &c_b24, &cv[cv_offset], ldcv, &c__[
		    c_offset], ldc, &c_b25, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

	    if (discr) {

/*              Compute X and SCALE satisfying */

/*              AV'*X*A - X = -SCALE*CV'*C. */

/*              Additional workspace needed: 2*NV. */

		sb04py_("T", "N", &c_n1, nv, n, &av[av_offset], ldav, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &dwork[kw], &
			ierr, (ftnlen)1, (ftnlen)1);
		if (ierr != 0) {
		    *info = 3;
		    return 0;
		}

/*              Construct C <- DV'*C + BV'*X*A/SCALE, */
/*                        D <- DV'*D + BV'*X*B/SCALE. */

/*              Additional workspace needed: MAX( P*N, P*M ). */

/*              C <- DV'*C. */

		dgemm_("T", "N", p, n, p, &c_b31, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, n, &dwork[kw], p, &c__[c_offset], ldc, (
			ftnlen)1);

/*              D <- DV'*D. */

		dgemm_("T", "N", p, m, p, &c_b31, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (
			ftnlen)1);

/*              C <- C + BV'*X*A/SCALE. */

		d__1 = 1. / scale;
		dgemm_("T", "N", p, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[1], &ldw, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dgemm_("N", "N", p, n, n, &c_b31, &dwork[kw], p, &a[a_offset],
			 lda, &c_b31, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)
			1);

/*              D <- D + BV'*X*B/SCALE. */

		dgemm_("N", "N", p, m, n, &c_b31, &dwork[kw], p, &b[b_offset],
			 ldb, &c_b31, &d__[d_offset], ldd, (ftnlen)1, (ftnlen)
			1);
	    } else {

/*              Compute X and SCALE satisfying */

/*              AV'*X + X*A + SCALE*CV'*C = 0. */

		dtrsyl_("T", "N", &c__1, nv, n, &av[av_offset], ldav, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &ierr, (
			ftnlen)1, (ftnlen)1);
		if (ierr != 0) {
		    *info = 3;
		    return 0;
		}

/*              Construct C and D. */
/*              Additional workspace needed: MAX( P*N, P*M ). */

/*              Construct C <- BV'*X/SCALE + DV'*C. */

		dgemm_("T", "N", p, n, p, &c_b31, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, n, &dwork[kw], p, &c__[c_offset], ldc, (
			ftnlen)1);
		d__1 = 1. / scale;
		dgemm_("T", "N", p, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[1], &ldw, &c_b31, &c__[c_offset], ldc, (ftnlen)
			1, (ftnlen)1);

/*              Construct D <- DV'*D. */

		dgemm_("T", "N", p, m, p, &c_b31, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (
			ftnlen)1);
	    }
	} else {

/*           Compute the projection of V*G. */

/*           Total workspace needed:  NV*N + MAX( P*N, P*M ). */

/*           Compute -BV*C. */
/*           Workspace needed: NV*N. */

	    dgemm_("N", "N", nv, n, p, &c_b24, &bv[bv_offset], ldbv, &c__[
		    c_offset], ldc, &c_b25, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

/*           Compute X and SCALE satisfying */

/*           AV*X - X*A + SCALE*BV*C = 0. */

	    dtrsyl_("N", "N", &c_n1, nv, n, &av[av_offset], ldav, &a[a_offset]
		    , lda, &dwork[1], &ldw, &scale, &ierr, (ftnlen)1, (ftnlen)
		    1);
	    if (ierr != 0) {
		*info = 3;
		return 0;
	    }

/*           Construct C <- CV*X/SCALE + DV*C. */

	    dgemm_("N", "N", p, n, p, &c_b31, &dv[dv_offset], lddv, &c__[
		    c_offset], ldc, &c_b25, &dwork[kw], p, (ftnlen)1, (ftnlen)
		    1);
	    dlacpy_("F", p, n, &dwork[kw], p, &c__[c_offset], ldc, (ftnlen)1);
	    d__1 = 1. / scale;
	    dgemm_("N", "N", p, n, nv, &d__1, &cv[cv_offset], ldcv, &dwork[1],
		     &ldw, &c_b31, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);

/*           Construct D <- DV*D. */

	    dgemm_("N", "N", p, m, p, &c_b31, &dv[dv_offset], lddv, &d__[
		    d_offset], ldd, &c_b25, &dwork[kw], p, (ftnlen)1, (ftnlen)
		    1);
	    dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (ftnlen)1);
	}
    }

    if (rightw) {
	ldwn = max(*n,1);
	kw = *n * *nw + 1;
	if (conjs) {

/*           Compute the projection of G*conj(W) or of conj(V)*G*conj(W). */

/*           Total workspace needed:  NW*N + MAX( b, M*N, P*M ), where */
/*                                    b = 0,    if DICO = 'C', */
/*                                    b = 2*NW, if DICO = 'D'. */

/*           Compute -BW*B'. */
/*           Workspace needed: N*NW. */

	    ldw = max(*nw,1);
	    dgemm_("N", "T", nw, n, m, &c_b24, &bw[bw_offset], ldbw, &b[
		    b_offset], ldb, &c_b25, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

	    if (discr) {

/*              Compute Y' and SCALE satisfying */

/*              AW*Y'*A' - Y' = -SCALE*BW*B'. */

/*              Additional workspace needed: 2*NW. */

		sb04py_("N", "T", &c_n1, nw, n, &aw[aw_offset], ldaw, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &dwork[kw], &
			ierr, (ftnlen)1, (ftnlen)1);
		if (ierr != 0) {
		    *info = 4;
		    return 0;
		}

/*              Construct B <- B*DW' + A*Y*CW'/SCALE, */
/*                        D <- D*DW' + C*Y*CW'/SCALE. */

/*              Additional workspace needed: MAX( N*M, P*M ). */

/*              B <- B*DW'. */

		dgemm_("N", "T", n, m, m, &c_b31, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b25, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
		dlacpy_("F", n, m, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)1);

/*              D <- D*DW'. */

		dgemm_("N", "T", p, m, m, &c_b31, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (
			ftnlen)1);

/*              B <- B + A*Y*CW'/SCALE. */

		d__1 = 1. / scale;
		dgemm_("T", "T", n, m, nw, &d__1, &dwork[1], &ldw, &cw[
			cw_offset], ldcw, &c_b25, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
		dgemm_("N", "N", n, m, n, &c_b31, &a[a_offset], lda, &dwork[
			kw], &ldwn, &c_b31, &b[b_offset], ldb, (ftnlen)1, (
			ftnlen)1);

/*              D <- D + C*Y*CW'/SCALE. */

		dgemm_("N", "N", p, m, n, &c_b31, &c__[c_offset], ldc, &dwork[
			kw], &ldwn, &c_b31, &d__[d_offset], ldd, (ftnlen)1, (
			ftnlen)1);
	    } else {

/*              Compute Y' and SCALE satisfying */

/*              AW*Y' + Y'*A' + SCALE*BW*B' = 0. */

		dtrsyl_("N", "T", &c__1, nw, n, &aw[aw_offset], ldaw, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &ierr, (
			ftnlen)1, (ftnlen)1);
		if (ierr != 0) {
		    *info = 4;
		    return 0;
		}

/*              Construct B and D. */
/*              Additional workspace needed: MAX( N*M, P*M ). */

/*              Construct B <- B*DW' + Y*CW'/SCALE. */

		dgemm_("N", "T", n, m, m, &c_b31, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b25, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
		dlacpy_("F", n, m, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)1);
		d__1 = 1. / scale;
		dgemm_("T", "T", n, m, nw, &d__1, &dwork[1], &ldw, &cw[
			cw_offset], ldcw, &c_b31, &b[b_offset], ldb, (ftnlen)
			1, (ftnlen)1);

/*              D <- D*DW'. */

		dgemm_("N", "T", p, m, m, &c_b31, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b25, &dwork[kw], p, (ftnlen)1, (
			ftnlen)1);
		dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (
			ftnlen)1);
	    }
	} else {

/*           Compute the projection of G*W or of V*G*W. */

/*           Total workspace needed:  NW*N + MAX( M*N, P*M ). */

/*           Compute B*CW. */
/*           Workspace needed: N*NW. */

	    dgemm_("N", "N", n, nw, m, &c_b31, &b[b_offset], ldb, &cw[
		    cw_offset], ldcw, &c_b25, &dwork[1], &ldwn, (ftnlen)1, (
		    ftnlen)1);

/*           Compute Y and SCALE satisfying */

/*           A*Y - Y*AW - SCALE*B*CW = 0. */

	    dtrsyl_("N", "N", &c_n1, n, nw, &a[a_offset], lda, &aw[aw_offset],
		     ldaw, &dwork[1], &ldwn, &scale, &ierr, (ftnlen)1, (
		    ftnlen)1);
	    if (ierr != 0) {
		*info = 4;
		return 0;
	    }

/*           Construct B and D. */
/*           Additional workspace needed: MAX( N*M, P*M ). */
/*           Construct B <- B*DW + Y*BW/SCALE. */

	    dgemm_("N", "N", n, m, m, &c_b31, &b[b_offset], ldb, &dw[
		    dw_offset], lddw, &c_b25, &dwork[kw], &ldwn, (ftnlen)1, (
		    ftnlen)1);
	    dlacpy_("F", n, m, &dwork[kw], &ldwn, &b[b_offset], ldb, (ftnlen)
		    1);
	    d__1 = 1. / scale;
	    dgemm_("N", "N", n, m, nw, &d__1, &dwork[1], &ldwn, &bw[bw_offset]
		    , ldbw, &c_b31, &b[b_offset], ldb, (ftnlen)1, (ftnlen)1);

/*           D <- D*DW. */

	    dgemm_("N", "N", p, m, m, &c_b31, &d__[d_offset], ldd, &dw[
		    dw_offset], lddw, &c_b25, &dwork[kw], p, (ftnlen)1, (
		    ftnlen)1);
	    dlacpy_("F", p, m, &dwork[kw], p, &d__[d_offset], ldd, (ftnlen)1);
	}
    }

/* Computing MAX */
    d__1 = work, d__2 = (doublereal) lw;
    dwork[1] = max(d__1,d__2);

    return 0;
/* *** Last line of AB09KX *** */
} /* ab09kx_ */

