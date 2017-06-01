/* AB09ID.f -- translated by f2c (version 20100827).
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

static doublereal c_b31 = 0.;

/* Subroutine */ int ab09id_(char *dico, char *jobc, char *jobo, char *job, 
	char *weight, char *equil, char *ordsel, integer *n, integer *m, 
	integer *p, integer *nv, integer *pv, integer *nw, integer *mw, 
	integer *nr, doublereal *alpha, doublereal *alphac, doublereal *
	alphao, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *av, integer *ldav, doublereal *bv, integer *ldbv, 
	doublereal *cv, integer *ldcv, doublereal *dv, integer *lddv, 
	doublereal *aw, integer *ldaw, doublereal *bw, integer *ldbw, 
	doublereal *cw, integer *ldcw, doublereal *dw, integer *lddw, integer 
	*ns, doublereal *hsv, doublereal *tol1, doublereal *tol2, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *iwarn, integer *
	info, ftnlen dico_len, ftnlen jobc_len, ftnlen jobo_len, ftnlen 
	job_len, ftnlen weight_len, ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, aw_dim1, aw_offset, b_dim1, 
	    b_offset, bv_dim1, bv_offset, bw_dim1, bw_offset, c_dim1, 
	    c_offset, cv_dim1, cv_offset, cw_dim1, cw_offset, d_dim1, 
	    d_offset, dv_dim1, dv_offset, dw_dim1, dw_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ki, kl, kt, ku, kw, lw, nn, nu, nu1;
    static logical bal;
    static integer lcf;
    static logical bta;
    static integer kbr, kcr, kdr, kbv;
    static logical spa;
    static integer kbw, kcv, kcw, kdv, kti, ldw, nra, nmr, nnq, nnr, nnv, nnw,
	     nvr, nwr, ppv, ierr;
    extern /* Subroutine */ int sb08cd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), sb08dd_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, ftnlen), tb01id_(char *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static logical scale;
    extern /* Subroutine */ int tb01kd_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen), tb01pd_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), ab09ix_(char *, char *, char *, char *, integer *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), ab09iy_(char *, char *
	    , char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr, leftw;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk;
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     state-space representation (A,B,C,D) by using the frequency */
/*     weighted square-root or balancing-free square-root */
/*     Balance & Truncate (B&T) or Singular Perturbation Approximation */
/*     (SPA) model reduction methods. The algorithm tries to minimize */
/*     the norm of the frequency-weighted error */

/*           ||V*(G-Gr)*W|| */

/*     where G and Gr are the transfer-function matrices of the original */
/*     and reduced order models, respectively, and V and W are */
/*     frequency-weighting transfer-function matrices. V and W must not */
/*     have poles on the imaginary axis for a continuous-time */
/*     system or on the unit circle for a discrete-time system. */
/*     If G is unstable, only the ALPHA-stable part of G is reduced. */
/*     In case of possible pole-zero cancellations in V*G and/or G*W, */
/*     the absolute values of parameters ALPHAO and/or ALPHAC must be */
/*     different from 1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

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

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'L':  only left weighting V is used (W = I); */
/*             = 'R':  only right weighting W is used (V = I); */
/*             = 'B':  both left and right weightings V and W are used. */

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
/*             The order of the original state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

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

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. For a system with NU ALPHA-unstable */
/*             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N), */
/*             NR is set as follows: if ORDSEL = 'F', NR is equal to */
/*             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order */
/*             on entry, NMIN is the number of frequency-weighted Hankel */
/*             singular values greater than NS*EPS*S1, EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH) */
/*             and S1 is the largest Hankel singular value (computed */
/*             in HSV(1)); NR can be further reduced to ensure */
/*             HSV(NR-NU) > HSV(NR+1-NU); */
/*             if ORDSEL = 'A', NR is the sum of NU and the number of */
/*             Hankel singular values greater than MAX(TOL1,NS*EPS*S1). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix A. For a continuous-time */
/*             system (DICO = 'C'), ALPHA <= 0 is the boundary value for */
/*             the real parts of eigenvalues, while for a discrete-time */
/*             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary. */

/*     ALPHAC  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted controllability Grammian (see METHOD); */
/*             ABS(ALPHAC) <= 1. */

/*     ALPHAO  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted observability Grammian (see METHOD); */
/*             ABS(ALPHAO) <= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system. */
/*             The resulting A has a block-diagonal form with two blocks. */
/*             For a system with NU ALPHA-unstable eigenvalues and */
/*             NS ALPHA-stable eigenvalues (NU+NS = N), the leading */
/*             NU-by-NU block contains the unreduced part of A */
/*             corresponding to ALPHA-unstable eigenvalues. */
/*             The trailing (NR+NS-N)-by-(NR+NS-N) block contains */
/*             the reduced part of A corresponding to ALPHA-stable */
/*             eigenvalues. */

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

/*     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV */
/*             part of this array must contain the state matrix AV of */
/*             the system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NVR-by-NVR part of this array */
/*             contains the state matrix of a minimal realization of V */
/*             in a real Schur form. NVR is returned in IWORK(2). */
/*             AV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDAV    INTEGER */
/*             The leading dimension of array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part */
/*             of this array must contain the input matrix BV of the */
/*             system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NVR-by-P part of this array contains */
/*             the input matrix of a minimal realization of V. */
/*             BV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDBV    INTEGER */
/*             The leading dimension of array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading PV-by-NV */
/*             part of this array must contain the output matrix CV of */
/*             the system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading PV-by-NVR part of this array */
/*             contains the output matrix of a minimal realization of V. */
/*             CV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDCV    INTEGER */
/*             The leading dimension of array CV. */
/*             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this */
/*             array must contain the feedthrough matrix DV of the system */
/*             with the transfer-function matrix V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDDV    INTEGER */
/*             The leading dimension of array DV. */
/*             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW */
/*             part of this array must contain the state matrix AW of */
/*             the system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NWR-by-NWR part of this array */
/*             contains the state matrix of a minimal realization of W */
/*             in a real Schur form. NWR is returned in IWORK(3). */
/*             AW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDAW    INTEGER */
/*             The leading dimension of array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,MW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-MW */
/*             part of this array must contain the input matrix BW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NWR-by-MW part of this array */
/*             contains the input matrix of a minimal realization of W. */
/*             BW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDBW    INTEGER */
/*             The leading dimension of array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part */
/*             of this array must contain the output matrix CW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading M-by-NWR part of this array contains */
/*             the output matrix of a minimal realization of W. */
/*             CW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDCW    INTEGER */
/*             The leading dimension of array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this */
/*             array must contain the feedthrough matrix DW of the system */
/*             with the transfer-function matrix W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDDW    INTEGER */
/*             The leading dimension of array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of this array contain */
/*             the frequency-weighted Hankel singular values, ordered */
/*             decreasingly, of the ALPHA-stable part of the original */
/*             system. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*S1, where c is a constant in the */
/*             interval [0.00001,0.001], and S1 is the largest */
/*             frequency-weighted Hankel singular value of the */
/*             ALPHA-stable part of the original system (computed */
/*             in HSV(1)). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS*S1, where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS*S1. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension */
/*             ( MAX( 3, LIWRK1, LIWRK2, LIWRK3 ) ), where */
/*             LIWRK1 = 0,             if JOB = 'B'; */
/*             LIWRK1 = N,             if JOB = 'F'; */
/*             LIWRK1 = 2*N,           if JOB = 'S' or 'P'; */
/*             LIWRK2 = 0,             if WEIGHT = 'R' or 'N' or  NV = 0; */
/*             LIWRK2 = NV+MAX(P,PV),  if WEIGHT = 'L' or 'B' and NV > 0; */
/*             LIWRK3 = 0,             if WEIGHT = 'L' or 'N' or  NW = 0; */
/*             LIWRK3 = NW+MAX(M,MW),  if WEIGHT = 'R' or 'B' and NW > 0. */
/*             On exit, if INFO = 0, IWORK(1) contains the order of a */
/*             minimal realization of the stable part of the system, */
/*             IWORK(2) and IWORK(3) contain the actual orders */
/*             of the state space realizations of V and W, respectively. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LMINL, LMINR, LRCF, */
/*                            2*N*N + MAX( 1, LLEFT, LRIGHT, 2*N*N+5*N, */
/*                                         N*MAX(M,P) ) ), */
/*             where */
/*             LMINL  = 0, if WEIGHT = 'R' or 'N' or NV = 0; otherwise, */
/*             LMINL  = MAX(LLCF,NV+MAX(NV,3*P))           if P =  PV; */
/*             LMINL  = MAX(P,PV)*(2*NV+MAX(P,PV))+ */
/*                      MAX(LLCF,NV+MAX(NV,3*P,3*PV))      if P <> PV; */
/*             LRCF   = 0, and */
/*             LMINR  = 0, if WEIGHT = 'L' or 'N' or NW = 0; otherwise, */
/*             LMINR  = NW+MAX(NW,3*M)                     if M =  MW; */
/*             LMINR  = 2*NW*MAX(M,MW)+NW+MAX(NW,3*M,3*MW) if M <> MW; */
/*             LLCF   = PV*(NV+PV)+PV*NV+MAX(NV*(NV+5), PV*(PV+2), */
/*                                           4*PV, 4*P); */
/*             LRCF   = MW*(NW+MW)+MAX(NW*(NW+5),MW*(MW+2),4*MW,4*M) */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than NSMIN, the sum of the order of the */
/*                   ALPHA-unstable part and the order of a minimal */
/*                   realization of the ALPHA-stable part of the given */
/*                   system; in this case, the resulting NR is set equal */
/*                   to NSMIN; */
/*             = 2:  with ORDSEL = 'F', the selected order NR corresponds */
/*                   to repeated singular values for the ALPHA-stable */
/*                   part, which are neither all included nor all */
/*                   excluded from the reduced model; in this case, the */
/*                   resulting NR is automatically decreased to exclude */
/*                   all repeated singular values; */
/*             = 3:  with ORDSEL = 'F', the selected order NR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given system; in this case NR is set equal to the */
/*                   order of the ALPHA-unstable part. */
/*             = 10+K:  K violations of the numerical stability condition */
/*                   occured during the assignment of eigenvalues in the */
/*                   SLICOT Library routines SB08CD and/or SB08DD. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of the ordered real Schur form of A */
/*                   failed; */
/*             = 2:  the separation of the ALPHA-stable/unstable */
/*                   diagonal blocks failed because of very close */
/*                   eigenvalues; */
/*             = 3:  the reduction to a real Schur form of the state */
/*                   matrix of a minimal realization of V failed; */
/*             = 4:  a failure was detected during the ordering of the */
/*                   real Schur form of the state matrix of a minimal */
/*                   realization of V or in the iterative process to */
/*                   compute a left coprime factorization with inner */
/*                   denominator; */
/*             = 5:  if DICO = 'C' and the matrix AV has an observable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' and */
/*                   AV has an observable eigenvalue on the unit circle; */
/*             = 6:  the reduction to a real Schur form of the state */
/*                   matrix of a minimal realization of W failed; */
/*             = 7:  a failure was detected during the ordering of the */
/*                   real Schur form of the state matrix of a minimal */
/*                   realization of W or in the iterative process to */
/*                   compute a right coprime factorization with inner */
/*                   denominator; */
/*             = 8:  if DICO = 'C' and the matrix AW has a controllable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' and */
/*                   AW has a controllable eigenvalue on the unit circle; */
/*             = 9:  the computation of eigenvalues failed; */
/*             = 10: the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let G be the transfer-function matrix of the original */
/*     linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                          (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09ID determines */
/*     the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                      (2) */

/*     such that the corresponding transfer-function matrix Gr minimizes */
/*     the norm of the frequency-weighted error */

/*             V*(G-Gr)*W,                                    (3) */

/*     where V and W are transfer-function matrices without poles on the */
/*     imaginary axis in continuous-time case or on the unit circle in */
/*     discrete-time case. */

/*     The following procedure is used to reduce G: */

/*     1) Decompose additively G, of order N, as */

/*          G = G1 + G2, */

/*        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and */
/*        G2 = (A2,B2,C2,0), of order NU, has only ALPHA-unstable poles. */

/*     2) Compute for G1 a B&T or SPA frequency-weighted approximation */
/*        G1r of order NR-NU using the combination method or the */
/*        modified combination method of [4]. */

/*     3) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     For the frequency-weighted reduction of the ALPHA-stable part, */
/*     several methods described in [4] can be employed in conjunction */
/*     with the combination method and modified combination method */
/*     proposed in [4]. */

/*     If JOB = 'B', the square-root B&T method is used. */
/*     If JOB = 'F', the balancing-free square-root version of the */
/*     B&T method is used. */
/*     If JOB = 'S', the square-root version of the SPA method is used. */
/*     If JOB = 'P', the balancing-free square-root version of the */
/*     SPA method is used. */

/*     For each of these methods, left and right truncation matrices */
/*     are determined using the Cholesky factors of an input */
/*     frequency-weighted controllability Grammian P and an output */
/*     frequency-weighted observability Grammian Q. */
/*     P and Q are computed from the controllability Grammian Pi of G*W */
/*     and the observability Grammian Qo of V*G. Using special */
/*     realizations of G*W and V*G, Pi and Qo are computed in the */
/*     partitioned forms */

/*           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) , */
/*                ( P12' P22 )               ( Q12' Q22 ) */

/*     where P11 and Q11 are the leading N-by-N parts of Pi and Qo, */
/*     respectively. Let P0 and Q0 be non-negative definite matrices */
/*     defined below */
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
/*     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds */
/*     to the method of Lin and Chiu [2,3]. */

/*     If JOBC = 'S' and ALPHAC = 1, no pole-zero cancellations must */
/*     occur in G*W. If JOBO = 'S' and ALPHAO = 1, no pole-zero */
/*     cancellations must occur in V*G. The presence of pole-zero */
/*     cancellations leads to meaningless results and must be avoided. */

/*     The frequency-weighted Hankel singular values HSV(1), ...., */
/*     HSV(N) are computed as the square roots of the eigenvalues */
/*     of the product P*Q. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. 23-th CDC, Las Vegas, pp. 127-132, 1984. */

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

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root */
/*     techniques. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*              Sep. 2001. */

/*     KEYWORDS */

/*     Frequency weighting, model reduction, multivariable system, */
/*     state-space model, state-space representation. */

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
    scale = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
    frwght = leftw || rightw;

    lw = 1;
    nn = *n * *n;
    nnv = *n + *nv;
    nnw = *n + *nw;
    ppv = max(*p,*pv);

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
/* Computing MAX */
    i__1 = lw, i__2 = (nn << 1) + *n * 5, i__1 = max(i__1,i__2), i__2 = *n * 
	    max(*m,*p);
    lw = (nn << 1) + max(i__1,i__2);

    if (leftw && *nv > 0) {
/* Computing MAX */
	i__1 = *nv * (*nv + 5), i__2 = *pv * (*pv + 2), i__1 = max(i__1,i__2),
		 i__2 = ppv << 2;
	lcf = *pv * (*nv + *pv) + *pv * *nv + max(i__1,i__2);
	if (*pv == *p) {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *nv, i__4 = *p * 3;
	    i__1 = max(lw,lcf), i__2 = *nv + max(i__3,i__4);
	    lw = max(i__1,i__2);
	} else {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	    i__5 = *nv, i__6 = ppv * 3;
	    i__3 = lcf, i__4 = *nv + max(i__5,i__6);
	    i__1 = lw, i__2 = ppv * ((*nv << 1) + ppv) + max(i__3,i__4);
	    lw = max(i__1,i__2);
	}
    }

    if (rightw && *nw > 0) {
	if (*mw == *m) {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *nw, i__4 = *m * 3;
	    i__1 = lw, i__2 = *nw + max(i__3,i__4);
	    lw = max(i__1,i__2);
	} else {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *nw, i__4 = *m * 3, i__3 = max(i__3,i__4), i__4 = *mw * 3;
	    i__1 = lw, i__2 = (*nw << 1) * max(*m,*mw) + *nw + max(i__3,i__4);
	    lw = max(i__1,i__2);
	}
/* Computing MAX */
/* Computing MAX */
	i__3 = *nw * (*nw + 5), i__4 = *mw * (*mw + 2), i__3 = max(i__3,i__4),
		 i__4 = *mw << 2, i__3 = max(i__3,i__4), i__4 = *m << 2;
	i__1 = lw, i__2 = *mw * (*nw + *mw) + max(i__3,i__4);
	lw = max(i__1,i__2);
    }

/*     Check the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (lsame_(jobc, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobc, 
	    "E", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lsame_(jobo, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobo, 
	    "E", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (bta || spa)) {
	*info = -4;
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -5;
    } else if (! (scale || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -6;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -7;
    } else if (*n < 0) {
	*info = -8;
    } else if (*m < 0) {
	*info = -9;
    } else if (*p < 0) {
	*info = -10;
    } else if (*nv < 0) {
	*info = -11;
    } else if (*pv < 0) {
	*info = -12;
    } else if (*nw < 0) {
	*info = -13;
    } else if (*mw < 0) {
	*info = -14;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -15;
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
	*info = -16;
    } else if (abs(*alphac) > 1.) {
	*info = -17;
    } else if (abs(*alphao) > 1.) {
	*info = -18;
    } else if (*lda < max(1,*n)) {
	*info = -20;
    } else if (*ldb < max(1,*n)) {
	*info = -22;
    } else if (*ldc < max(1,*p)) {
	*info = -24;
    } else if (*ldd < max(1,*p)) {
	*info = -26;
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
	*info = -28;
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
	*info = -30;
    } else if (*ldcv < 1 || leftw && *ldcv < *pv) {
	*info = -32;
    } else if (*lddv < 1 || leftw && *lddv < *pv) {
	*info = -34;
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
	*info = -36;
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
	*info = -38;
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
	*info = -40;
    } else if (*lddw < 1 || rightw && *lddw < *m) {
	*info = -42;
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
	*info = -46;
    } else if (*ldwork < lw) {
	*info = -49;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09ID", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	*nr = 0;
	*ns = 0;
	iwork[1] = 0;
	iwork[2] = *nv;
	iwork[3] = *nw;
	dwork[1] = 1.;
	return 0;
    }

    if (scale) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

	maxred = 100.;
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
    }

/*     Correct the value of ALPHA to ensure stability. */

    alpwrk = *alpha;
    if (discr) {
	if (*alpha == 1.) {
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
	}
    } else {
	if (*alpha == 0.) {
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
	}
    }

/*     Allocate working storage. */

    ku = 1;
    kl = ku + nn;
    ki = kl + *n;
    kw = ki + *n;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation, A <- inv(T)*A*T, and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

    i__1 = *ldwork - kw + 1;
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kl], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

    if (ierr != 0) {
	if (ierr != 3) {
	    *info = 1;
	} else {
	    *info = 2;
	}
	return 0;
    }

    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Determine NRA, the desired order for the reduction of stable part. */

    iwarnl = 0;
    *ns = *n - nu;
    if (fixord) {
/* Computing MAX */
	i__1 = 0, i__2 = *nr - nu;
	nra = max(i__1,i__2);
	if (*nr < nu) {
	    iwarnl = 3;
	}
    } else {
	nra = 0;
    }

/*     Finish if only unstable part is present. */

    if (*ns == 0) {
	*nr = nu;
	dwork[1] = (doublereal) wrkopt;
	iwork[1] = 0;
	iwork[2] = *nv;
	iwork[3] = *nw;
	return 0;
    }

    nvr = *nv;
    if (leftw && *nv > 0) {

/*        Compute a left-coprime factorization with inner denominator */
/*        of a minimal realization of V. The resulting AV is in */
/*        real Schur form. */
/*        Workspace needed:   real  LV+MAX( 1, LCF, */
/*                                          NV + MAX( NV, 3*P, 3*PV ) ), */
/*                                  where */
/*                                  LV = 0 if P = PV and */
/*                                  LV = MAX(P,PV)*(2*NV+MAX(P,PV)) */
/*                                         otherwise; */
/*                                  LCF = PV*(NV+PV) + */
/*                                        MAX( 1, PV*NV + MAX( NV*(NV+5), */
/*                                             PV*(PV+2),4*PV,4*P ) ); */
/*                                  prefer larger; */
/*                          integer NV + MAX(P,PV). */

	if (*p == *pv) {
	    kw = 1;
	    tb01pd_("Minimal", "Scale", nv, p, pv, &av[av_offset], ldav, &bv[
		    bv_offset], ldbv, &cv[cv_offset], ldcv, &nvr, &c_b31, &
		    iwork[1], &dwork[1], ldwork, info, (ftnlen)7, (ftnlen)5);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	    wrkopt = max(i__1,i__2);
	    kbr = 1;
	    kdr = kbr + *pv * nvr;
	    kw = kdr + *pv * *pv;
	    i__1 = max(1,nvr);
	    i__2 = *ldwork - kw + 1;
	    sb08cd_(dico, &nvr, p, pv, &av[av_offset], ldav, &bv[bv_offset], 
		    ldbv, &cv[cv_offset], ldcv, &dv[dv_offset], lddv, &nnq, &
		    nnr, &dwork[kbr], &i__1, &dwork[kdr], pv, &c_b31, &dwork[
		    kw], &i__2, iwarn, &ierr, (ftnlen)1);
	} else {
	    ldw = max(*p,*pv);
	    kbv = 1;
	    kcv = kbv + *nv * ldw;
	    kw = kcv + *nv * ldw;
	    dlacpy_("Full", nv, p, &bv[bv_offset], ldbv, &dwork[kbv], nv, (
		    ftnlen)4);
	    dlacpy_("Full", pv, nv, &cv[cv_offset], ldcv, &dwork[kcv], &ldw, (
		    ftnlen)4);
	    i__1 = *ldwork - kw + 1;
	    tb01pd_("Minimal", "Scale", nv, p, pv, &av[av_offset], ldav, &
		    dwork[kbv], nv, &dwork[kcv], &ldw, &nvr, &c_b31, &iwork[1]
		    , &dwork[kw], &i__1, info, (ftnlen)7, (ftnlen)5);
	    kdv = kw;
	    kbr = kdv + ldw * ldw;
	    kdr = kbr + *pv * nvr;
	    kw = kdr + *pv * *pv;
	    dlacpy_("Full", pv, p, &dv[dv_offset], lddv, &dwork[kdv], &ldw, (
		    ftnlen)4);
	    i__1 = max(1,nvr);
	    i__2 = *ldwork - kw + 1;
	    sb08cd_(dico, &nvr, p, pv, &av[av_offset], ldav, &dwork[kbv], nv, 
		    &dwork[kcv], &ldw, &dwork[kdv], &ldw, &nnq, &nnr, &dwork[
		    kbr], &i__1, &dwork[kdr], pv, &c_b31, &dwork[kw], &i__2, 
		    iwarn, &ierr, (ftnlen)1);
	    dlacpy_("Full", &nvr, p, &dwork[kbv], nv, &bv[bv_offset], ldbv, (
		    ftnlen)4);
	    dlacpy_("Full", pv, &nvr, &dwork[kcv], &ldw, &cv[cv_offset], ldcv,
		     (ftnlen)4);
	    dlacpy_("Full", pv, p, &dwork[kdv], &ldw, &dv[dv_offset], lddv, (
		    ftnlen)4);
	}
	if (ierr != 0) {
	    *info = ierr + 2;
	    return 0;
	}
	nvr = nnq;
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
	if (*iwarn > 0) {
	    *iwarn += 10;
	}
    }

    nwr = *nw;
    if (rightw && *nw > 0) {

/*        Compute a minimal realization of W. */
/*        Workspace needed:   real  LW+MAX(1, NW + MAX(NW, 3*M, 3*MW)); */
/*                                  where */
/*                                  LW = 0,              if M = MW and */
/*                                  LW = 2*NW*MAX(M,MW), otherwise; */
/*                                  prefer larger; */
/*                          integer NW + MAX(M,MW). */

	if (*m == *mw) {
	    kw = 1;
	    tb01pd_("Minimal", "Scale", nw, mw, m, &aw[aw_offset], ldaw, &bw[
		    bw_offset], ldbw, &cw[cw_offset], ldcw, &nwr, &c_b31, &
		    iwork[1], &dwork[1], ldwork, info, (ftnlen)7, (ftnlen)5);
	} else {
	    ldw = max(*m,*mw);
	    kbw = 1;
	    kcw = kbw + *nw * ldw;
	    kw = kcw + *nw * ldw;
	    dlacpy_("Full", nw, mw, &bw[bw_offset], ldbw, &dwork[kbw], nw, (
		    ftnlen)4);
	    dlacpy_("Full", m, nw, &cw[cw_offset], ldcw, &dwork[kcw], &ldw, (
		    ftnlen)4);
	    i__1 = *ldwork - kw + 1;
	    tb01pd_("Minimal", "Scale", nw, mw, m, &aw[aw_offset], ldaw, &
		    dwork[kbw], nw, &dwork[kcw], &ldw, &nwr, &c_b31, &iwork[1]
		    , &dwork[kw], &i__1, info, (ftnlen)7, (ftnlen)5);
	    dlacpy_("Full", &nwr, mw, &dwork[kbw], nw, &bw[bw_offset], ldbw, (
		    ftnlen)4);
	    dlacpy_("Full", m, &nwr, &dwork[kcw], &ldw, &cw[cw_offset], ldcw, 
		    (ftnlen)4);
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
    }

    if (rightw && nwr > 0) {

/*        Compute a right-coprime factorization with inner denominator */
/*        of the minimal realization of W. The resulting AW is in */
/*        real Schur form. */

/*        Workspace needed:  MW*(NW+MW) + */
/*                           MAX( 1, NW*(NW+5), MW*(MW+2), 4*MW, 4*M ); */
/*                           prefer larger. */

	ldw = max(1,*mw);
	kcr = 1;
	kdr = kcr + nwr * ldw;
	kw = kdr + *mw * ldw;
	i__1 = *ldwork - kw + 1;
	sb08dd_(dico, &nwr, mw, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw,
		 &cw[cw_offset], ldcw, &dw[dw_offset], lddw, &nnq, &nnr, &
		dwork[kcr], &ldw, &dwork[kdr], &ldw, &c_b31, &dwork[kw], &
		i__1, iwarn, &ierr, (ftnlen)1);
	if (ierr != 0) {
	    *info = ierr + 5;
	    return 0;
	}
	nwr = nnq;
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
	wrkopt = max(i__1,i__2);
	if (*iwarn > 0) {
	    *iwarn += 10;
	}
    }

    nu1 = nu + 1;

/*     Allocate working storage. */

    kt = 1;
    kti = kt + nn;
    kw = kti + nn;

/*     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R */
/*     of the controllability and observability Grammians, respectively. */
/*     Real workspace:    need  2*N*N + MAX( 1, LLEFT, LRIGHT ), */
/*             where */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
/*                        prefer larger. */

    i__1 = *ldwork - kw + 1;
    ab09iy_(dico, jobc, jobo, weight, ns, m, p, &nvr, pv, &nwr, mw, alphac, 
	    alphao, &a[nu1 + nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[
	    nu1 * c_dim1 + 1], ldc, &av[av_offset], ldav, &bv[bv_offset], 
	    ldbv, &cv[cv_offset], ldcv, &dv[dv_offset], lddv, &aw[aw_offset], 
	    ldaw, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[dw_offset], 
	    lddw, &scalec, &scaleo, &dwork[kti], n, &dwork[kt], n, &dwork[kw],
	     &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (ierr != 0) {
	*info = 9;
	return 0;
    }
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
    wrkopt = max(i__1,i__2);

/*     Compute a BTA or SPA of the stable part. */
/*     Real workspace:  need  2*N*N + MAX( 1, 2*N*N+5*N, N*MAX(M,P) ). */

    i__1 = *ldwork - kw + 1;
    ab09ix_(dico, job, "Schur", ordsel, ns, m, p, &nra, &scalec, &scaleo, &a[
	    nu1 + nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
	    c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kti], n, &dwork[kt],
	     n, &nmr, &hsv[1], tol1, tol2, &iwork[1], &dwork[kw], &i__1, 
	    iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)5, (ftnlen)1);
    *iwarn = max(*iwarn,iwarnl);
    if (ierr != 0) {
	*info = 10;
	return 0;
    }
    *nr = nra + nu;

/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
    dwork[1] = (doublereal) max(i__1,i__2);
    iwork[1] = nmr;
    iwork[2] = nvr;
    iwork[3] = nwr;

    return 0;
/* *** Last line of AB09ID *** */
} /* ab09id_ */

