/* AB09IY.f -- translated by f2c (version 20100827).
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

static doublereal c_b16 = 0.;
static doublereal c_b20 = 1.;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b43 = -1.;
static logical c_true = TRUE_;

/* Subroutine */ int ab09iy_(char *dico, char *jobc, char *jobo, char *weight,
	 integer *n, integer *m, integer *p, integer *nv, integer *pv, 
	integer *nw, integer *mw, doublereal *alphac, doublereal *alphao, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *av, integer *ldav, doublereal *bv, 
	integer *ldbv, doublereal *cv, integer *ldcv, doublereal *dv, integer 
	*lddv, doublereal *aw, integer *ldaw, doublereal *bw, integer *ldbw, 
	doublereal *cw, integer *ldcw, doublereal *dw, integer *lddw, 
	doublereal *scalec, doublereal *scaleo, doublereal *s, integer *lds, 
	doublereal *r__, integer *ldr, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen dico_len, ftnlen jobc_len, ftnlen jobo_len, 
	ftnlen weight_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, aw_dim1, aw_offset, b_dim1, 
	    b_offset, bv_dim1, bv_offset, bw_dim1, bw_offset, c_dim1, 
	    c_offset, cv_dim1, cv_offset, cw_dim1, cw_offset, dv_dim1, 
	    dv_offset, dw_dim1, dw_offset, r_dim1, r_offset, s_dim1, s_offset,
	     i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer ku, kw, lw, kaw, ldu;
    static doublereal dum[1], tol;
    static integer nnv, nnw, ierr, ktau;
    static doublereal work;
    static integer mbbar;
    extern /* Subroutine */ int mb04nd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), mb04od_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    static integer pcbar;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb01wd_(char *, char *, char *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dcopy_(integer *, doublereal *, integer *,
	     doublereal *, integer *);
    static logical leftw;
    extern /* Subroutine */ int dsyev_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical frwght, rightw;


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

/*     To compute for given state-space representations */
/*     (A,B,C,0), (AV,BV,CV,DV), and (AW,BW,CW,DW) of the */
/*     transfer-function matrices G, V and W, respectively, */
/*     the Cholesky factors of the frequency-weighted */
/*     controllability and observability Grammians corresponding */
/*     to a frequency-weighted model reduction problem. */
/*     G, V and W must be stable transfer-function matrices with */
/*     the state matrices A, AV, and AW in real Schur form. */
/*     It is assumed that the state space realizations (AV,BV,CV,DV) */
/*     and (AW,BW,CW,DW) are minimal. In case of possible pole-zero */
/*     cancellations in forming V*G and/or G*W, the parameters for the */
/*     choice of frequency-weighted Grammians ALPHAO and/or ALPHAC, */
/*     respectively, must be different from 1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G, V and W are continuous-time systems; */
/*             = 'D':  G, V and W are discrete-time systems. */

/*     JOBC    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted controllability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to a combination method [4] */
/*                    of the approaches of Enns [1] and Lin-Chiu [2,3]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [4]. */

/*     JOBO    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted observability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to a combination method [4] */
/*                    of the approaches of Enns [1] and Lin-Chiu [2,3]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [4]. */

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'L':  only left weighting V is used (W = I); */
/*             = 'R':  only right weighting W is used (V = I); */
/*             = 'B':  both left and right weightings V and W are used. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation of G, i.e., */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B and */
/*             the number of rows of the matrices CW and DW.  M >= 0. */
/*             M represents the dimension of the input vector of the */
/*             system with the transfer-function matrix G and */
/*             also the dimension of the output vector of the system */
/*             with the transfer-function matrix W. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrix C and the */
/*             number of columns of the matrices BV and DV.  P >= 0. */
/*             P represents the dimension of the output vector of the */
/*             system with the transfer-function matrix G and */
/*             also the dimension of the input vector of the system */
/*             with the transfer-function matrix V. */

/*     NV      (input) INTEGER */
/*             The order of the matrix AV. Also the number of rows of */
/*             the matrix BV and the number of columns of the matrix CV. */
/*             NV represents the dimension of the state vector of the */
/*             system with the transfer-function matrix V.  NV >= 0. */

/*     PV      (input) INTEGER */
/*             The number of rows of the matrices CV and DV.  PV >= 0. */
/*             PV represents the dimension of the output vector of the */
/*             system with the transfer-function matrix V. */

/*     NW      (input) INTEGER */
/*             The order of the matrix AW. Also the number of rows of */
/*             the matrix BW and the number of columns of the matrix CW. */
/*             NW represents the dimension of the state vector of the */
/*             system with the transfer-function matrix W.  NW >= 0. */

/*     MW      (input) INTEGER */
/*             The number of columns of the matrices BW and DW.  MW >= 0. */
/*             MW represents the dimension of the input vector of the */
/*             system with the transfer-function matrix W. */

/*     ALPHAC  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted controllability Grammian (see METHOD); */
/*             ABS(ALPHAC) <= 1. */

/*     ALPHAO  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted observability Grammian (see METHOD); */
/*             ABS(ALPHAO) <= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must */
/*             contain the state matrix A (of the system with the */
/*             transfer-function matrix G) in a real Schur form. */

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

/*     AV      (input) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             If WEIGHT = 'L' or 'B', the leading NV-by-NV part of this */
/*             array must contain the state matrix AV (of the system with */
/*             the transfer-function matrix V) in a real Schur form. */
/*             AV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDAV    INTEGER */
/*             The leading dimension of array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             If WEIGHT = 'L' or 'B', the leading NV-by-P part of this */
/*             array must contain the input matrix BV of the system with */
/*             the transfer-function matrix V. */
/*             BV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDBV    INTEGER */
/*             The leading dimension of array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-NV part of this */
/*             array must contain the output matrix CV of the system with */
/*             the transfer-function matrix V. */
/*             CV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDCV    INTEGER */
/*             The leading dimension of array CV. */
/*             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this */
/*             array must contain the feedthrough matrix DV of the system */
/*             with the transfer-function matrix V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDDV    INTEGER */
/*             The leading dimension of array DV. */
/*             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     AW      (input) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             If WEIGHT = 'R' or 'B', the leading NW-by-NW part of this */
/*             array must contain the state matrix AW (of the system with */
/*             the transfer-function matrix W) in a real Schur form. */
/*             AW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDAW    INTEGER */
/*             The leading dimension of array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input) DOUBLE PRECISION array, dimension (LDBW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading NW-by-MW part of this */
/*             array must contain the input matrix BW of the system with */
/*             the transfer-function matrix W. */
/*             BW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDBW    INTEGER */
/*             The leading dimension of array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-NW part of this */
/*             array must contain the output matrix CW of the system with */
/*             the transfer-function matrix W. */
/*             CW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDCW    INTEGER */
/*             The leading dimension of array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this */
/*             array must contain the feedthrough matrix DW of the system */
/*             with the transfer-function matrix W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDDW    INTEGER */
/*             The leading dimension of array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     SCALEC  (output) DOUBLE PRECISION */
/*             Scaling factor for the controllability Grammian in (1) */
/*             or (3). See METHOD. */

/*     SCALEO  (output) DOUBLE PRECISION */
/*             Scaling factor for the observability Grammian in (2) */
/*             or (4). See METHOD. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor S of the frequency-weighted */
/*             cotrollability Grammian P = S*S'. See METHOD. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor R of the frequency-weighted */
/*             observability Grammian Q = R'*R. See METHOD. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, LLEFT, LRIGHT ), */
/*             where */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the state matrices A and/or AV are not stable or */
/*                   not in a real Schur form; */
/*             = 2:  if the state matrices A and/or AW are not stable or */
/*                   not in a real Schur form; */
/*             = 3:  eigenvalues computation failure. */

/*     METHOD */

/*     Let Pi = Si*Si' and Qo = Ro'*Ro be the Cholesky factored */
/*     controllability and observability Grammians satisfying */
/*     in the continuous-time case */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0,       (1) */

/*            Ao'*Qo + Qo*Ao +  scaleo^2*Co'*Co = 0,       (2) */

/*     and in the discrete-time case */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0,       (3) */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0,       (4) */

/*     where */

/*           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) , */
/*                ( 0   Aw  )          (  Bw  ) */

/*           Ao = (  A   0  ) ,   Co = ( Dv*C  Cv ) . */
/*                ( Bv*C Av ) */

/*     Consider the partitioned Grammians */

/*           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) , */
/*                ( P12' P22 )               ( Q12' Q22 ) */

/*     where P11 and Q11 are the leading N-by-N parts of Pi and Qo, */
/*     respectively, and let P0 and Q0 be non-negative definite matrices */
/*     defined in the combination method [4] */
/*                                        -1 */
/*            P0 = P11 - ALPHAC**2*P12*P22 *P21 , */
/*                                        -1 */
/*            Q0 = Q11 - ALPHAO**2*Q12*Q22 *Q21. */

/*     The frequency-weighted controllability and observability */
/*     Grammians, P and Q, respectively, are defined as follows: */
/*     P = P0 if JOBC = 'S' (standard combination method [4]); */
/*     P = P1 >= P0 if JOBC = 'E', where P1 is the controllability */
/*     Grammian defined to enforce stability for a modified combination */
/*     method of [4]; */
/*     Q = Q0 if JOBO = 'S' (standard combination method [4]); */
/*     Q = Q1 >= Q0 if JOBO = 'E', where Q1 is the observability */
/*     Grammian defined to enforce stability for a modified combination */
/*     method of [4]. */

/*     If JOBC = JOBO = 'S' and ALPHAC = ALPHAO = 0, the choice of */
/*     Grammians corresponds to the method of Enns [1], while if */
/*     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds to the */
/*     method of Lin and Chiu [2,3]. */

/*     The routine computes directly the Cholesky factors S and R */
/*     such that P = S*S' and Q = R'*R according to formulas */
/*     developed in [4]. No matrix inversions are involved. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. CDC, Las Vegas, pp. 127-132, 1984. */

/*     [2] Lin, C.-A. and Chiu, T.-Y. */
/*         Model reduction via frequency-weighted balanced realization. */
/*         Control Theory and Advanced Technology, vol. 8, */
/*         pp. 341-351, 1992. */

/*     [3] Sreeram, V., Anderson, B.D.O and Madievski, A.G. */
/*         New results on frequency weighted balanced reduction */
/*         technique. */
/*         Proc. ACC, Seattle, Washington, pp. 4004-4009, 1995. */

/*     [4] Varga, A. and Anderson, B.D.O. */
/*         Square-root balancing-free methods for the frequency-weighted */
/*         balancing related model reduction. */
/*         (report in preparation) */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2001. */

/*     KEYWORDS */

/*     Frequency weighting, model reduction, multivariable system, */
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
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --dwork;

    /* Function Body */
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
    frwght = leftw || rightw;

    *info = 0;
    lw = 1;
    nnv = *n + *nv;
    nnw = *n + *nw;
    if (leftw && *pv > 0) {
/* Computing MAX */
	i__1 = lw, i__2 = nnv * (nnv + max(nnv,*pv) + 5);
	lw = max(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = lw, i__2 = *n * (*p + 5);
	lw = max(i__1,i__2);
    }
    if (rightw && *mw > 0) {
/* Computing MAX */
	i__1 = lw, i__2 = nnw * (nnw + max(nnw,*mw) + 5);
	lw = max(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = lw, i__2 = *n * (*m + 5);
	lw = max(i__1,i__2);
    }

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
    } else if (*nv < 0) {
	*info = -8;
    } else if (*pv < 0) {
	*info = -9;
    } else if (*nw < 0) {
	*info = -10;
    } else if (*mw < 0) {
	*info = -11;
    } else if (abs(*alphac) > 1.) {
	*info = -12;
    } else if (abs(*alphao) > 1.) {
	*info = -13;
    } else if (*lda < max(1,*n)) {
	*info = -15;
    } else if (*ldb < max(1,*n)) {
	*info = -17;
    } else if (*ldc < max(1,*p)) {
	*info = -19;
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
	*info = -21;
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
	*info = -23;
    } else if (*ldcv < 1 || leftw && *ldcv < *pv) {
	*info = -25;
    } else if (*lddv < 1 || leftw && *lddv < *pv) {
	*info = -27;
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
	*info = -29;
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
	*info = -31;
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
	*info = -33;
    } else if (*lddw < 1 || rightw && *lddw < *m) {
	*info = -35;
    } else if (*lds < max(1,*n)) {
	*info = -39;
    } else if (*ldr < max(1,*n)) {
	*info = -41;
    } else if (*ldwork < lw) {
	*info = -43;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09IY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *scalec = 1.;
    *scaleo = 1.;
/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	dwork[1] = 1.;
	return 0;
    }

    work = 1.;
    if (leftw && *pv > 0) {

/*        Build the extended permuted matrices */

/*           Ao = ( Av  Bv*C ) ,   Co = ( Cv Dv*C ) . */
/*                ( 0     A  ) */

	kaw = 1;
	ku = kaw + nnv * nnv;
	ldu = max(nnv,*pv);
	dlacpy_("Full", nv, nv, &av[av_offset], ldav, &dwork[kaw], &nnv, (
		ftnlen)4);
	dlaset_("Full", n, nv, &c_b16, &c_b16, &dwork[kaw + *nv], &nnv, (
		ftnlen)4);
	dgemm_("No-transpose", "No-transpose", nv, n, p, &c_b20, &bv[
		bv_offset], ldbv, &c__[c_offset], ldc, &c_b16, &dwork[kaw + 
		nnv * *nv], &nnv, (ftnlen)12, (ftnlen)12);
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw + nnv * *nv + *nv]
		, &nnv, (ftnlen)4);

	dlacpy_("Full", pv, nv, &cv[cv_offset], ldcv, &dwork[ku], &ldu, (
		ftnlen)4);
	dgemm_("No-transpose", "No-transpose", pv, n, p, &c_b20, &dv[
		dv_offset], lddv, &c__[c_offset], ldc, &c_b16, &dwork[ku + 
		ldu * *nv], &ldu, (ftnlen)12, (ftnlen)12);

/*        Solve for the Cholesky factor Ro of Qo, Qo = Ro'*Ro, */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ao'*Qo + Qo*Ao  +  scaleo^2*Co'*Co = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0. */

/*        Workspace:  need   (N+NV)*(N+NV+MAX(N+NV,PV)+5); */
/*                           prefer larger. */

	ktau = ku + ldu * nnv;
	kw = ktau + nnv;

	i__1 = *ldwork - kw + 1;
	sb03ou_(&discr, &c_false, &nnv, pv, &dwork[kaw], &nnv, &dwork[ku], &
		ldu, &dwork[ktau], &dwork[ku], &ldu, scaleo, &dwork[kw], &
		i__1, &ierr);

	if (ierr != 0) {
	    *info = 1;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);

/*        Partition Ro as Ro = ( R11 R12 ) and compute R such that */
/*                             (  0  R22 ) */

/*        R'*R = R22'*R22 + (1-ALPHAO**2)*R12'*R12. */

	kw = ku + ldu * *nv + *nv;
	dlacpy_("Upper", n, n, &dwork[kw], &ldu, &r__[r_offset], ldr, (ftnlen)
		5);
	if (*alphao != 0.) {
	    t = sqrt(1. - *alphao * *alphao);
	    i__1 = ku + ldu * (nnv - 1);
	    i__2 = ldu;
	    for (j = ku + ldu * *nv; i__2 < 0 ? j >= i__1 : j <= i__1; j += 
		    i__2) {
		dscal_(nv, &t, &dwork[j], &c__1);
/* L10: */
	    }
	}
	if (*alphao < 1. && *nv > 0) {
	    ktau = 1;
	    mb04od_("Full", n, &c__0, nv, &r__[r_offset], ldr, &dwork[ku + 
		    ldu * *nv], &ldu, dum, &c__1, dum, &c__1, &dwork[ktau], &
		    dwork[kw], (ftnlen)4);

	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		dwork[j] = r__[j + j * r_dim1];
		i__1 = j;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (dwork[i__] < 0.) {
			r__[i__ + j * r_dim1] = -r__[i__ + j * r_dim1];
		    }
/* L20: */
		}
/* L30: */
	    }

	}

	if (lsame_(jobo, "E", (ftnlen)1, (ftnlen)1) && *alphao < 1.) {

/*           Form Y = -A'*(R'*R)-(R'*R)*A if DICO = 'C', or */
/*                Y = -A'*(R'*R)*A+(R'*R) if DICO = 'D'. */

	    dlacpy_("Upper", n, n, &r__[r_offset], ldr, &dwork[ku], n, (
		    ftnlen)5);
	    mb01wd_(dico, "Upper", "No-transpose", "Hessenberg", n, &c_b43, &
		    c_b16, &r__[r_offset], ldr, &dwork[kaw + nnv * *nv + *nv],
		     &nnv, &dwork[ku], n, &ierr, (ftnlen)1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)10);

/*           Compute the eigendecomposition of Y as Y = Z*Sigma*Z'. */

	    ku = *n + 1;
	    i__2 = *ldwork - *n;
	    dsyev_("Vectors", "Upper", n, &r__[r_offset], ldr, &dwork[1], &
		    dwork[ku], &i__2, &ierr, (ftnlen)7, (ftnlen)5);
	    if (ierr > 0) {
		*info = 3;
		return 0;
	    }
/* Computing MAX */
	    d__1 = work, d__2 = dwork[ku] + (doublereal) (*n);
	    work = max(d__1,d__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 <= 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
	    d__2 = abs(dwork[1]), d__3 = (d__1 = dwork[*n], abs(d__1));
	    tol = max(d__2,d__3) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form C = [ sqrt(Sigma2)*Z2' ] */

	    pcbar = 0;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		if (dwork[j] > tol) {
		    d__1 = sqrt(dwork[j]);
		    dscal_(n, &d__1, &r__[j * r_dim1 + 1], &c__1);
		    dcopy_(n, &r__[j * r_dim1 + 1], &c__1, &dwork[ku + pcbar],
			     n);
		    ++pcbar;
		}
/* L40: */
	    }

/*           Solve for the Cholesky factor R of Q, Q = R'*R, */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                      _  _ */
/*                   A'*Q + Q*A  +  t^2*C'*C = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                      _  _ */
/*                   A'*Q*A - Q  +  t^2*C'*C = 0. */

/*           Workspace:  need   N*(N + 6); */
/*                              prefer larger. */

	    ktau = ku + *n * *n;
	    kw = ktau + *n;

	    i__2 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_false, n, &pcbar, &a[a_offset], lda, &dwork[ku]
		    , n, &dwork[ktau], &r__[r_offset], ldr, &t, &dwork[kw], &
		    i__2, &ierr);
	    if (ierr != 0) {
		*info = 1;
		return 0;
	    }
	    *scaleo *= t;
/* Computing MAX */
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	    work = max(d__1,d__2);
	}

    } else {

/*        Solve for the Cholesky factor R of Q, Q = R'*R, */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            A'*Q + Q*A  +  scaleo^2*C'*C = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            A'*Q*A - Q +  scaleo^2*C'*C = 0. */

/*        Workspace:  need   N*(P + 5); */
/*                           prefer larger. */

	ku = 1;
	ktau = ku + *p * *n;
	kw = ktau + *n;

	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);
	i__2 = *ldwork - kw + 1;
	sb03ou_(&discr, &c_false, n, p, &a[a_offset], lda, &dwork[ku], p, &
		dwork[ktau], &r__[r_offset], ldr, scaleo, &dwork[kw], &i__2, &
		ierr);
	if (ierr != 0) {
	    *info = 1;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);
    }

    if (rightw && *mw > 0) {

/*        Build the extended matrices */

/*           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) . */
/*                ( 0   Aw  )          (  Bw  ) */

	kaw = 1;
	ku = kaw + nnw * nnw;
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw], &nnw, (ftnlen)4)
		;
	dlaset_("Full", nw, n, &c_b16, &c_b16, &dwork[kaw + *n], &nnw, (
		ftnlen)4);
	dgemm_("No-transpose", "No-transpose", n, nw, m, &c_b20, &b[b_offset],
		 ldb, &cw[cw_offset], ldcw, &c_b16, &dwork[kaw + nnw * *n], &
		nnw, (ftnlen)12, (ftnlen)12);
	dlacpy_("Full", nw, nw, &aw[aw_offset], ldaw, &dwork[kaw + nnw * *n + 
		*n], &nnw, (ftnlen)4);

	dgemm_("No-transpose", "No-transpose", n, mw, m, &c_b20, &b[b_offset],
		 ldb, &dw[dw_offset], lddw, &c_b16, &dwork[ku], &nnw, (ftnlen)
		12, (ftnlen)12);
	dlacpy_("Full", nw, mw, &bw[bw_offset], ldbw, &dwork[ku + *n], &nnw, (
		ftnlen)4);

/*        Solve for the Cholesky factor Si of Pi, Pi = Si*Si', */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0. */

/*        Workspace:  need   (N+NW)*(N+NW+MAX(N+NW,MW)+5); */
/*                           prefer larger. */

	ktau = ku + nnw * max(nnw,*mw);
	kw = ktau + nnw;

	i__2 = *ldwork - kw + 1;
	sb03ou_(&discr, &c_true, &nnw, mw, &dwork[kaw], &nnw, &dwork[ku], &
		nnw, &dwork[ktau], &dwork[ku], &nnw, scalec, &dwork[kw], &
		i__2, &ierr);

	if (ierr != 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);

/*        Partition Si as Si = ( S11 S12 ) and compute S such that */
/*                             (  0  S22 ) */

/*        S*S' = S11*S11' + (1-ALPHAC**2)*S12*S12'. */

	dlacpy_("Upper", n, n, &dwork[ku], &nnw, &s[s_offset], lds, (ftnlen)5)
		;
	if (*alphac != 0.) {
	    t = sqrt(1. - *alphac * *alphac);
	    i__2 = ku + nnw * (nnw - 1);
	    i__1 = nnw;
	    for (j = ku + nnw * *n; i__1 < 0 ? j >= i__2 : j <= i__2; j += 
		    i__1) {
		dscal_(n, &t, &dwork[j], &c__1);
/* L50: */
	    }
	}
	if (*alphac < 1. && *nw > 0) {
	    ktau = *n * nnw + 1;
	    kw = ktau + *n;
	    mb04nd_("Full", n, &c__0, nw, &s[s_offset], lds, &dwork[ku + nnw *
		     *n], &nnw, dum, &c__1, dum, &c__1, &dwork[ktau], &dwork[
		    kw], (ftnlen)4);

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (s[j + j * s_dim1] < 0.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			s[i__ + j * s_dim1] = -s[i__ + j * s_dim1];
/* L60: */
		    }
		}
/* L70: */
	    }
	}

	if (lsame_(jobc, "E", (ftnlen)1, (ftnlen)1) && *alphac < 1.) {

/*           Form X = -A*(S*S')-(S*S')*A' if DICO = 'C', or */
/*                X = -A*(S*S')*A'+(S*S') if DICO = 'D'. */

	    dlacpy_("Upper", n, n, &s[s_offset], lds, &dwork[ku], n, (ftnlen)
		    5);
	    mb01wd_(dico, "Upper", "Transpose", "Hessenberg", n, &c_b43, &
		    c_b16, &s[s_offset], lds, &dwork[kaw], &nnw, &dwork[ku], 
		    n, &ierr, (ftnlen)1, (ftnlen)5, (ftnlen)9, (ftnlen)10);

/*           Compute the eigendecomposition of X as X = Z*Sigma*Z'. */

	    ku = *n + 1;
	    i__1 = *ldwork - *n;
	    dsyev_("Vectors", "Upper", n, &s[s_offset], lds, &dwork[1], &
		    dwork[ku], &i__1, &ierr, (ftnlen)7, (ftnlen)5);
	    if (ierr > 0) {
		*info = 3;
		return 0;
	    }
/* Computing MAX */
	    d__1 = work, d__2 = dwork[ku] + (doublereal) (*n);
	    work = max(d__1,d__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 =< 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
	    d__2 = abs(dwork[1]), d__3 = (d__1 = dwork[*n], abs(d__1));
	    tol = max(d__2,d__3) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form B = [ Z2*sqrt(Sigma2) ] */

	    mbbar = 0;
	    i__ = ku;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (dwork[j] > tol) {
		    ++mbbar;
		    d__1 = sqrt(dwork[j]);
		    dscal_(n, &d__1, &s[j * s_dim1 + 1], &c__1);
		    dcopy_(n, &s[j * s_dim1 + 1], &c__1, &dwork[i__], &c__1);
		    i__ += *n;
		}
/* L80: */
	    }

/*           Solve for the Cholesky factor S of P, P = S*S', */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                      _ _ */
/*                   A*P + P*A'  +  t^2*B*B' = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                      _ _ */
/*                   A*P*A' - P  +  t^2*B*B' = 0. */

/*           Workspace:  need   maximum N*(N + 6); */
/*                              prefer larger. */

	    ktau = ku + mbbar * *n;
	    kw = ktau + *n;

	    i__1 = *ldwork - kw + 1;
	    sb03ou_(&discr, &c_true, n, &mbbar, &a[a_offset], lda, &dwork[ku],
		     n, &dwork[ktau], &s[s_offset], lds, &t, &dwork[kw], &
		    i__1, &ierr);
	    if (ierr != 0) {
		*info = 2;
		return 0;
	    }
	    *scalec *= t;
/* Computing MAX */
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	    work = max(d__1,d__2);
	}

    } else {

/*        Solve for the Cholesky factor S of P, P = S*S', */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            A*P + P*A' +  scalec^2*B*B' = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            A*P*A' - P +  scalec^2*B*B' = 0. */

/*        Workspace:  need   N*(M+5); */
/*                           prefer larger. */

	ku = 1;
	ktau = ku + *n * *m;
	kw = ktau + *n;

	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
	i__1 = *ldwork - kw + 1;
	sb03ou_(&discr, &c_true, n, m, &a[a_offset], lda, &dwork[ku], n, &
		dwork[ktau], &s[s_offset], lds, scalec, &dwork[kw], &i__1, &
		ierr);
	if (ierr != 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
	work = max(d__1,d__2);
    }

/*     Save optimal workspace. */

    dwork[1] = work;

    return 0;
/* *** Last line of AB09IY *** */
} /* ab09iy_ */

