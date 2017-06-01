/* IB03BD.f -- translated by f2c (version 20100827).
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

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b20 = -1.;
static integer c__3 = 3;
static integer c__0 = 0;
static doublereal c_b63 = 100.;
static doublereal c_b64 = 0.;

/* Subroutine */ int ib03bd_(char *init, integer *nobr, integer *m, integer *
	l, integer *nsmp, integer *n, integer *nn, integer *itmax1, integer *
	itmax2, integer *nprint, doublereal *u, integer *ldu, doublereal *y, 
	integer *ldy, doublereal *x, integer *lx, doublereal *tol1, 
	doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen init_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9, i__10, i__11;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, z__, n2, ac, bd, ia, ib, ik, ml, iq, ir, is, iv, 
	    ix, ns, nx, iw1, iw2, iw3, ix0, ldr, bsn, mno, isv, iry, ldac, 
	    isad;
    static doublereal seed[4], rcnd[16];
    static integer ipar[7], nfev, njev, lnol, nsml, lths, nths;
    static doublereal work[4];
    static logical init1, init2;
    extern /* Subroutine */ int ib01ad_(char *, char *, char *, char *, char *
	    , char *, integer *, integer *, integer *, integer *, doublereal *
	    , integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen), ib01bd_(char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, logical *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), ib01cd_(char *, char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int md03ba_(), md03bb_();
    extern /* Subroutine */ int md03bd_(char *, char *, char *, U_fp, U_fp, 
	    U_fp, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern /* Subroutine */ int nf01be_(), nf01bf_();
    static integer idiag;
    extern /* Subroutine */ int nf01bp_(), nf01bs_();
    static integer ircnd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01vd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);
    static integer infol, lipar;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), tf01mx_(integer *, integer *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static logical bwork[1];
    extern /* Subroutine */ int tb01vy_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);
    static integer jwork, ircndb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer iwarnl, wrkopt;


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

/*     To compute a set of parameters for approximating a Wiener system */
/*     in a least-squares sense, using a neural network approach and a */
/*     MINPACK-like Levenberg-Marquardt algorithm. The Wiener system */
/*     consists of a linear part and a static nonlinearity, and it is */
/*     represented as */

/*        x(t+1) = A*x(t) + B*u(t) */
/*        z(t)   = C*x(t) + D*u(t), */

/*        y(t)   = f(z(t),wb(1:L)), */

/*     where t = 1, 2, ..., NSMP, and f is a nonlinear function, */
/*     evaluated by the SLICOT Library routine NF01AY. The parameter */
/*     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ), */
/*     where theta corresponds to the linear part, and wb(i), i = 1 : L, */
/*     correspond to the nonlinear part. See SLICOT Library routine */
/*     NF01AD for further details. */

/*     The sum of squares of the error functions, defined by */

/*        e(t) = y(t) - Y(t),  t = 1, 2, ..., NSMP, */

/*     is minimized, where Y(t) is the measured output vector. The */
/*     functions and their Jacobian matrices are evaluated by SLICOT */
/*     Library routine NF01BF (the FCN routine in the call of MD03BD). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     INIT    CHARACTER*1 */
/*             Specifies which parts have to be initialized, as follows: */
/*             = 'L' : initialize the linear part only, X already */
/*                     contains an initial approximation of the */
/*                     nonlinearity; */
/*             = 'S' : initialize the static nonlinearity only, X */
/*                     already contains an initial approximation of the */
/*                     linear part; */
/*             = 'B' : initialize both linear and nonlinear parts; */
/*             = 'N' : do not initialize anything, X already contains */
/*                     an initial approximation. */
/*             If INIT = 'S' or 'B', the error functions for the */
/*             nonlinear part, and their Jacobian matrices, are evaluated */
/*             by SLICOT Library routine NF01BE (used as a second FCN */
/*             routine in the MD03BD call for the initialization step, */
/*             see METHOD). */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             If INIT = 'L' or 'B', NOBR is the number of block rows, s, */
/*             in the input and output block Hankel matrices to be */
/*             processed for estimating the linear part.  NOBR > 0. */
/*             (In the MOESP theory,  NOBR  should be larger than  n, */
/*             the estimated dimension of state vector.) */
/*             This parameter is ignored if INIT is 'S' or 'N'. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L >= 0, and L > 0, if */
/*             INIT = 'L' or 'B'. */

/*     NSMP    (input) INTEGER */
/*             The number of input and output samples, t.  NSMP >= 0, and */
/*             NSMP >= 2*(M+L+1)*NOBR - 1, if INIT = 'L' or 'B'. */

/*     N       (input/output) INTEGER */
/*             The order of the linear part. */
/*             If INIT = 'L' or 'B', and N < 0 on entry, the order is */
/*             assumed unknown and it will be found by the routine. */
/*             Otherwise, the input value will be used. If INIT = 'S' */
/*             or 'N', N must be non-negative. The values N >= NOBR, */
/*             or N = 0, are not acceptable if INIT = 'L' or 'B'. */

/*     NN      (input) INTEGER */
/*             The number of neurons which shall be used to approximate */
/*             the nonlinear part.  NN >= 0. */

/*     ITMAX1  (input) INTEGER */
/*             The maximum number of iterations for the initialization of */
/*             the static nonlinearity. */
/*             This parameter is ignored if INIT is 'N' or 'L'. */
/*             Otherwise, ITMAX1 >= 0. */

/*     ITMAX2  (input) INTEGER */
/*             The maximum number of iterations.  ITMAX2 >= 0. */

/*     NPRINT  (input) INTEGER */
/*             This parameter enables controlled printing of iterates if */
/*             it is positive. In this case, FCN is called with IFLAG = 0 */
/*             at the beginning of the first iteration and every NPRINT */
/*             iterations thereafter and immediately prior to return, */
/*             and the current error norm is printed. Other intermediate */
/*             results could be printed by modifying the corresponding */
/*             FCN routine (NF01BE and/or NF01BF). If NPRINT <= 0, no */
/*             special calls of FCN with IFLAG = 0 are made. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU, M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             set of input samples, */
/*             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NSMP). */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             set of output samples, */
/*             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ). */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,NSMP). */

/*     X       (input/output) DOUBLE PRECISION array dimension (LX) */
/*             On entry, if INIT = 'L', the leading (NN*(L+2) + 1)*L part */
/*             of this array must contain the initial parameters for */
/*             the nonlinear part of the system. */
/*             On entry, if INIT = 'S', the elements lin1 : lin2 of this */
/*             array must contain the initial parameters for the linear */
/*             part of the system, corresponding to the output normal */
/*             form, computed by SLICOT Library routine TB01VD, where */
/*                lin1 = (NN*(L+2) + 1)*L + 1; */
/*                lin2 = (NN*(L+2) + 1)*L + N*(L+M+1) + L*M. */
/*             On entry, if INIT = 'N', the elements 1 : lin2 of this */
/*             array must contain the initial parameters for the */
/*             nonlinear part followed by the initial parameters for the */
/*             linear part of the system, as specified above. */
/*             This array need not be set on entry if INIT = 'B'. */
/*             On exit, the elements 1 : lin2 of this array contain the */
/*             optimal parameters for the nonlinear part followed by the */
/*             optimal parameters for the linear part of the system, as */
/*             specified above. */

/*     LX      (input/output) INTEGER */
/*             On entry, this parameter must contain the intended length */
/*             of X. If N >= 0, then LX >= NX := lin2 (see parameter X). */
/*             If N is unknown (N < 0 on entry), a large enough estimate */
/*             of N should be used in the formula of lin2. */
/*             On exit, if N < 0 on entry, but LX is not large enough, */
/*             then this parameter contains the actual length of X, */
/*             corresponding to the computed N. Otherwise, its value */
/*             is unchanged. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If INIT = 'S' or 'B' and TOL1 >= 0, TOL1 is the tolerance */
/*             which measures the relative error desired in the sum of */
/*             squares, as well as the relative error desired in the */
/*             approximate solution, for the initialization step of */
/*             nonlinear part. Termination occurs when either both the */
/*             actual and predicted relative reductions in the sum of */
/*             squares, or the relative error between two consecutive */
/*             iterates are at most TOL1. If the user sets  TOL1 < 0, */
/*             then  SQRT(EPS)  is used instead TOL1, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is ignored if INIT is 'N' or 'L'. */

/*     TOL2    DOUBLE PRECISION */
/*             If TOL2 >= 0, TOL2 is the tolerance which measures the */
/*             relative error desired in the sum of squares, as well as */
/*             the relative error desired in the approximate solution, */
/*             for the whole optimization process. Termination occurs */
/*             when either both the actual and predicted relative */
/*             reductions in the sum of squares, or the relative error */
/*             between two consecutive iterates are at most TOL2. If the */
/*             user sets TOL2 < 0, then  SQRT(EPS)  is used instead TOL2. */
/*             This default value could require many iterations, */
/*             especially if TOL1 is larger. If INIT = 'S' or 'B', it is */
/*             advisable that TOL2 be larger than TOL1, and spend more */
/*             time with cheaper iterations. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX( LIW1, LIW2, LIW3 )), where */
/*             LIW1 = LIW2 = 0,  if INIT = 'S' or 'N'; otherwise, */
/*             LIW1 = M+L; */
/*             LIW2 = MAX(M*NOBR+N,M*(N+L)); */
/*             LIW3 = 3+MAX(NN*(L+2)+2,NX+L), if INIT = 'S' or 'B'; */
/*             LIW3 = 3+NX+L,                 if INIT = 'L' or 'N'. */
/*             On output, if INFO = 0, IWORK(1) and IWORK(2) return the */
/*             (total) number of function and Jacobian evaluations, */
/*             respectively (including the initialization step, if it was */
/*             performed), and if INIT = 'L' or INIT = 'B', IWORK(3) */
/*             specifies how many locations of DWORK contain reciprocal */
/*             condition number estimates (see below); otherwise, */
/*             IWORK(3) = 0. If INFO = 0, the entries 4 to 3+NX of IWORK */
/*             define a permutation matrix P such that J*P = Q*R, where */
/*             J is the final calculated Jacobian, Q is an orthogonal */
/*             matrix (not stored), and R is upper triangular with */
/*             diagonal elements of nonincreasing magnitude (possibly */
/*             for each block column of J). Column j of P is column */
/*             IWORK(3+j) of the identity matrix. Moreover, the entries */
/*             4+NX:3+NX+L of this array contain the ranks of the final */
/*             submatrices S_k (see description of LMPARM in MD03BD). */

/*     DWORK   DOUBLE PRECISION array dimesion (LDWORK) */
/*             On entry, if desired, and if INIT = 'S' or 'B', the */
/*             entries DWORK(1:4) are set to initialize the random */
/*             numbers generator for the nonlinear part parameters (see */
/*             the description of the argument XINIT of SLICOT Library */
/*             routine MD03BD); this enables to obtain reproducible */
/*             results. The same seed is used for all outputs. */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the residual error norm (the */
/*             sum of squares), DWORK(3) returns the number of iterations */
/*             performed, and DWORK(4) returns the final Levenberg */
/*             factor, for optimizing the parameters of both the linear */
/*             part and the static nonlinearity part. If INIT = 'S' or */
/*             INIT = 'B' and INFO = 0, then the elements DWORK(5) to */
/*             DWORK(8) contain the corresponding four values for the */
/*             initialization step (see METHOD). (If L > 1, DWORK(8) */
/*             contains the maximum of the Levenberg factors for all */
/*             outputs.) If INIT = 'L' or INIT = 'B', and INFO = 0, */
/*             DWORK(9) to DWORK(8+IWORK(3)) contain reciprocal condition */
/*             number estimates set by SLICOT Library routines IB01AD, */
/*             IB01BD, and IB01CD. */
/*             On exit, if  INFO = -21,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             In the formulas below, N should be taken not larger than */
/*             NOBR - 1, if N < 0 on entry. */
/*             LDWORK = MAX( LW1, LW2, LW3, LW4 ), where */
/*             LW1 = 0, if INIT = 'S' or 'N'; otherwise, */
/*             LW1 = MAX( 2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR, */
/*                        4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) + */
/*                        MAX( LDW1, LDW2 ), */
/*                        (N+L)*(N+M) + N + N*N + 2 + N*(N+M+L) + */
/*                        MAX( 5*N, 2, MIN( LDW3, LDW4 ), LDW5, LDW6 ), */
/*                 where, */
/*                 LDW1 >= MAX( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N, */
/*                              L*NOBR*N + */
/*                              MAX( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L, */
/*                                   2*(L*NOBR-L)*N+N*N+8*N, */
/*                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) ) */
/*                 LDW2 >= 0,                                  if M = 0; */
/*                 LDW2 >= L*NOBR*N + M*NOBR*(N+L)*(M*(N+L)+1) + */
/*                         MAX( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0; */
/*                 LDW3 = NSMP*L*(N+1) + 2*N + MAX( 2*N*N, 4*N ), */
/*                 LDW4 = N*(N+1) + 2*N + */
/*                        MAX( N*L*(N+1) + 2*N*N + L*N, 4*N ); */
/*                 LDW5 = NSMP*L + (N+L)*(N+M) + 3*N+M+L; */
/*                 LDW6 = NSMP*L + (N+L)*(N+M) + N + */
/*                        MAX(1, N*N*L + N*L + N, N*N + */
/*                            MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), */
/*                                N*M)); */
/*             LW2 = LW3 = 0, if INIT = 'L' or 'N'; otherwise, */
/*             LW2 = NSMP*L + BSN + */
/*                   MAX( 4, NSMP + */
/*                           MAX( NSMP*BSN + MAX( 2*NN, 5*BSN + 1 ), */
/*                                BSN**2 + BSN + */
/*                                MAX( NSMP + 2*NN, 5*BSN ) ) ); */
/*             LW3 = MAX( LDW7, NSMP*L + (N+L)*(2*N+M) + 2*N ); */
/*                 LDW7 = NSMP*L + (N+L)*(N+M) + 3*N+M+L,  if M > 0; */
/*                 LDW7 = NSMP*L + (N+L)*N + 2*N+L,        if M = 0; */
/*             LW4 = NSMP*L + NX + */
/*                   MAX( 4, NSMP*L + */
/*                           MAX( NSMP*L*( BSN + LTHS ) + */
/*                                MAX( NSMP*L + L1, L2 + NX ), */
/*                                     NX*( BSN + LTHS ) + NX + */
/*                                     MAX( NSMP*L + L1, NX + L3 ) ) ), */
/*                  L0 = MAX( N*(N+L), N+M+L ),    if M > 0; */
/*                  L0 = MAX( N*(N+L), L ),        if M = 0; */
/*                  L1 = NSMP*L + MAX( 2*NN, (N+L)*(N+M) + 2*N + L0); */
/*                  L2 = 4*NX + 1,  if L <= 1 or BSN = 0; otherwise, */
/*                  L2 = BSN + MAX(3*BSN+1,LTHS); */
/*                  L2 = MAX(L2,4*LTHS+1),         if NSMP > BSN; */
/*                  L2 = MAX(L2,(NSMP-BSN)*(L-1)), if BSN < NSMP < 2*BSN; */
/*                  L3 = 4*NX,                     if L <= 1 or BSN = 0; */
/*                  L3 = LTHS*BSN + 2*NX + 2*MAX(BSN,LTHS), */
/*                                                 if L > 1 and BSN > 0, */
/*                  with BSN  = NN*( L + 2 ) + 1, */
/*                       LTHS = N*( L + M + 1 ) + L*M. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             < 0:  the user set IFLAG = IWARN in (one of) the */
/*                   subroutine(s) FCN, i.e., NF01BE, if INIT = 'S' */
/*                   or 'B', and/or NF01BF; this value cannot be returned */
/*                   without changing the FCN routine(s); */
/*                   otherwise, IWARN has the value k*100 + j*10 + i, */
/*                   where k is defined below, i refers to the whole */
/*                   optimization process, and j refers to the */
/*                   initialization step (j = 0, if INIT = 'L' or 'N'), */
/*                   and the possible values for i and j have the */
/*                   following meaning (where TOL* denotes TOL1 or TOL2, */
/*                   and similarly for ITMAX*): */
/*             = 1:  both actual and predicted relative reductions in */
/*                   the sum of squares are at most TOL*; */
/*             = 2:  relative error between two consecutive iterates is */
/*                   at most TOL*; */
/*             = 3:  conditions for i or j = 1 and i or j = 2 both hold; */
/*             = 4:  the cosine of the angle between the vector of error */
/*                   function values and any column of the Jacobian is at */
/*                   most EPS in absolute value; */
/*             = 5:  the number of iterations has reached ITMAX* without */
/*                   satisfying any convergence condition; */
/*             = 6:  TOL* is too small: no further reduction in the sum */
/*                   of squares is possible; */
/*             = 7:  TOL* is too small: no further improvement in the */
/*                   approximate solution X is possible; */
/*             = 8:  the vector of function values e is orthogonal to the */
/*                   columns of the Jacobian to machine precision. */
/*             The digit k is normally 0, but if INIT = 'L' or 'B', it */
/*             can have a value in the range 1 to 6 (see IB01AD, IB01BD */
/*             and IB01CD). In all these cases, the entries DWORK(1:4), */
/*             DWORK(5:8) (if INIT = 'S' or 'B'), and DWORK(9:8+IWORK(3)) */
/*             (if INIT = 'L' or 'B'), are set as described above. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*                   otherwise, INFO has the value k*100 + j*10 + i, */
/*                   where k is defined below, i refers to the whole */
/*                   optimization process, and j refers to the */
/*                   initialization step (j = 0, if INIT = 'L' or 'N'), */
/*                   and the possible values for i and j have the */
/*                   following meaning: */
/*             = 1:  the routine FCN returned with INFO <> 0 for */
/*                   IFLAG = 1; */
/*             = 2:  the routine FCN returned with INFO <> 0 for */
/*                   IFLAG = 2; */
/*             = 3:  the routine QRFACT returned with INFO <> 0; */
/*             = 4:  the routine LMPARM returned with INFO <> 0. */
/*             In addition, if INIT = 'L' or 'B', i could also be */
/*             = 5:  if a Lyapunov equation could not be solved; */
/*             = 6:  if the identified linear system is unstable; */
/*             = 7:  if the QR algorithm failed on the state matrix */
/*                   of the identified linear system. */
/*             QRFACT and LMPARM are generic names for SLICOT Library */
/*             routines NF01BS and NF01BP, respectively, for the whole */
/*             optimization process, and MD03BA and MD03BB, respectively, */
/*             for the initialization step (if INIT = 'S' or 'B'). */
/*             The digit k is normally 0, but if INIT = 'L' or 'B', it */
/*             can have a value in the range 1 to 10 (see IB01AD/IB01BD). */

/*     METHOD */

/*     If INIT = 'L' or 'B', the linear part of the system is */
/*     approximated using the combined MOESP and N4SID algorithm. If */
/*     necessary, this algorithm can also choose the order, but it is */
/*     advantageous if the order is already known. */

/*     If INIT = 'S' or 'B', the output of the approximated linear part */
/*     is computed and used to calculate an approximation of the static */
/*     nonlinearity using the Levenberg-Marquardt algorithm [1,3]. */
/*     This step is referred to as the (nonlinear) initialization step. */

/*     As last step, the Levenberg-Marquardt algorithm is used again to */
/*     optimize the parameters of the linear part and the static */
/*     nonlinearity as a whole. Therefore, it is necessary to parametrise */
/*     the matrices of the linear part. The output normal form [2] */
/*     parameterisation is used. */

/*     The Jacobian is computed analytically, for the nonlinear part, and */
/*     numerically, for the linear part. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     [2] Peeters, R.L.M., Hanzon, B., and Olivi, M. */
/*         Balanced realizations of discrete-time stable all-pass */
/*         systems and the tangential Schur algorithm. */
/*         Proceedings of the European Control Conference, */
/*         31 August - 3 September 1999, Karlsruhe, Germany. */
/*         Session CP-6, Discrete-time Systems, 1999. */

/*     [3] More, J.J. */
/*         The Levenberg-Marquardt algorithm: implementation and theory. */
/*         In Watson, G.A. (Ed.), Numerical Analysis, Lecture Notes in */
/*         Mathematics, vol. 630, Springer-Verlag, Berlin, Heidelberg */
/*         and New York, pp. 105-116, 1978. */

/*     NUMERICAL ASPECTS */

/*     The Levenberg-Marquardt algorithm described in [3] is scaling */
/*     invariant and globally convergent to (maybe local) minima. */
/*     The convergence rate near a local minimum is quadratic, if the */
/*     Jacobian is computed analytically, and linear, if the Jacobian */
/*     is computed numerically. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, March, 2002, Apr. 2002, Feb. 2004, March 2005. */

/*     KEYWORDS */

/*     Least-squares approximation, Levenberg-Marquardt algorithm, */
/*     matrix operations, optimization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     FACTOR is a scaling factor for variables (see MD03BD). */
/*     Condition estimation and internal scaling of variables are used */
/*     (see MD03BD). */
/*     Default tolerances are used in MD03BD for measuring the */
/*     orthogonality between the vector of function values and columns */
/*     of the Jacobian (GTOL), and for the rank estimations (TOL). */
/*     For INIT = 'L' or 'B', additional parameters are set: */
/*     The following six parameters are used in the call of IB01AD; */
/*     The following three parameters are used in the call of IB01BD; */
/*     The following two parameters are used in the call of IB01CD; */
/*     TOLN controls the estimated order in IB01AD (default value); */
/*     RCOND controls the rank decisions in IB01AD, IB01BD, and IB01CD */
/*     (default); */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --x;
    --iwork;
    --dwork;

    /* Function Body */
    init1 = lsame_(init, "B", (ftnlen)1, (ftnlen)1) || lsame_(init, "L", (
	    ftnlen)1, (ftnlen)1);
    init2 = lsame_(init, "B", (ftnlen)1, (ftnlen)1) || lsame_(init, "S", (
	    ftnlen)1, (ftnlen)1);

    ml = *m + *l;
    *info = 0;
    *iwarn = 0;
    if (! (init1 || init2 || lsame_(init, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (init1 && *nobr <= 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*l < 0 || init1 && *l == 0) {
	*info = -4;
    } else if (*nsmp < 0 || init1 && *nsmp < (ml + 1 << 1) * *nobr - 1) {
	*info = -5;
    } else if (*n < 0 && ! init1 || (*n == 0 || *n >= *nobr) && init1) {
	*info = -6;
    } else if (*nn < 0) {
	*info = -7;
    } else if (init2 && *itmax1 < 0) {
	*info = -8;
    } else if (*itmax2 < 0) {
	*info = -9;
    } else if (*ldu < max(1,*nsmp)) {
	*info = -12;
    } else if (*ldy < max(1,*nsmp)) {
	*info = -14;
    } else {
	lnol = *l * *nobr - *l;
	mno = *m * *nobr;
	bsn = *nn * (*l + 2) + 1;
	nths = bsn * *l;
	nsml = *nsmp * *l;
	if (*n > 0) {
	    ldac = *n + *l;
	    isad = ldac * (*n + *m);
	    n2 = *n * *n;
	}

/*        Check the workspace size. */

	jwork = 0;
	if (init1) {
/*           Workspace for IB01AD. */
	    jwork = (ml << 1) * *nobr * ((ml << 1) * (*nobr + 1) + 3) + *l * *
		    nobr;
	    if (*n > 0) {
/*              Workspace for IB01BD. */
/* Computing MAX */
/* Computing MAX */
		i__3 = lnol * *n + (*n << 1) + (*m + ml) * *nobr + *l, i__4 = 
			(lnol << 1) * *n + n2 + (*n << 3), i__3 = max(i__3,
			i__4), i__4 = *n + (mno + *n << 2) + 1, i__3 = max(
			i__3,i__4), i__4 = mno + *n * 3 + *l;
		i__1 = (lnol << 1) * *n + (*n << 1), i__2 = lnol * *n + n2 + *
			n * 7, i__1 = max(i__1,i__2), i__2 = *l * *nobr * *n 
			+ max(i__3,i__4);
		iw1 = max(i__1,i__2);
		if (*m > 0) {
/* Computing MAX */
/* Computing 2nd power */
		    i__3 = ldac;
		    i__1 = i__3 * i__3, i__2 = (*m << 2) * ldac + 1;
		    iw2 = *l * *nobr * *n + mno * ldac * (*m * ldac + 1) + 
			    max(i__1,i__2);
		} else {
		    iw2 = 0;
		}
/* Computing MAX */
/* Computing 2nd power */
		i__3 = (ml << 1) * *nobr;
		i__1 = jwork, i__2 = i__3 * i__3 + isad + max(iw1,iw2);
		jwork = max(i__1,i__2);
/*              Workspace for IB01CD. */
/* Computing MAX */
		i__1 = n2 << 1, i__2 = *n << 2;
		iw1 = nsml * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
		i__1 = *n * *l * (*n + 1) + (n2 << 1) + *l * *n, i__2 = *n << 
			2;
		iw2 = *n * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
		i__3 = *n * 5, i__3 = max(i__3,2), i__4 = min(iw1,iw2);
		i__1 = jwork, i__2 = isad + 2 + *n * (*n + 1 + ldac + *m) + 
			max(i__3,i__4);
		jwork = max(i__1,i__2);
/*              Workspace for TF01MX. */
/* Computing MAX */
		i__1 = jwork, i__2 = nsml + isad + ldac + (*n << 1) + *m;
		jwork = max(i__1,i__2);
/*              Workspace for TB01VD. */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
		i__5 = n2 + *n * max(*n,*l) + *n * 6 + min(*n,*l), i__6 = *n *
			 *m;
		i__3 = 1, i__4 = n2 * *l + *n * *l + *n, i__3 = max(i__3,i__4)
			, i__4 = n2 + max(i__5,i__6);
		i__1 = jwork, i__2 = nsml + isad + *n + max(i__3,i__4);
		jwork = max(i__1,i__2);
	    }
	}

	if (init2) {
/*           Workspace for MD03BD (initialization of the nonlinear part). */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	    i__7 = *nn << 1, i__8 = bsn * 5 + 1;
/* Computing 2nd power */
	    i__9 = bsn;
/* Computing MAX */
	    i__10 = *nsmp + (*nn << 1), i__11 = bsn * 5;
	    i__5 = *nsmp * bsn + max(i__7,i__8), i__6 = i__9 * i__9 + bsn + 
		    max(i__10,i__11);
	    i__3 = 4, i__4 = *nsmp + max(i__5,i__6);
	    i__1 = jwork, i__2 = nsml + bsn + max(i__3,i__4);
	    jwork = max(i__1,i__2);
	    if (*n > 0 && ! init1) {
/*              Workspace for TB01VY. */
/* Computing MAX */
		i__1 = jwork, i__2 = nsml + ldac * ((*n << 1) + *m) + (*n << 
			1);
		jwork = max(i__1,i__2);
/*              Workspace for TF01MX. */
		if (*m > 0) {
		    iw1 = *n + *m;
		} else {
		    iw1 = 0;
		}
/* Computing MAX */
		i__1 = jwork, i__2 = nsml + isad + iw1 + ldac + *n;
		jwork = max(i__1,i__2);
	    }
	}

	if (*n >= 0) {

/*           Find the number of parameters. */

	    lths = *n * (ml + 1) + *l * *m;
	    nx = nths + lths;

	    if (*lx < nx) {
		*info = -16;
		i__1 = -(*info);
		xerbla_("IB03BD", &i__1, (ftnlen)6);
		return 0;
	    }

/*           Workspace for MD03BD (whole optimization). */

	    if (*m > 0) {
		iw1 = ldac + *m;
	    } else {
		iw1 = *l;
	    }
/* Computing MAX */
/* Computing MAX */
	    i__3 = *n * ldac;
	    i__1 = *nn << 1, i__2 = isad + (*n << 1) + max(i__3,iw1);
	    iw1 = nsml + max(i__1,i__2);
	    if (*l <= 1 || bsn == 0) {
		iw3 = nx << 2;
		iw2 = iw3 + 1;
	    } else {
/* Computing MAX */
		i__1 = bsn * 3 + 1;
		iw2 = bsn + max(i__1,lths);
		if (*nsmp > bsn) {
/* Computing MAX */
		    i__1 = iw2, i__2 = (lths << 2) + 1;
		    iw2 = max(i__1,i__2);
		    if (*nsmp < bsn << 1) {
/* Computing MAX */
			i__1 = iw2, i__2 = (*nsmp - bsn) * (*l - 1);
			iw2 = max(i__1,i__2);
		    }
		}
		iw3 = lths * bsn + (nx << 1) + (max(bsn,lths) << 1);
	    }
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	    i__7 = nsml + iw1, i__8 = iw2 + nx;
/* Computing MAX */
	    i__9 = nsml + iw1, i__10 = nx + iw3;
	    i__5 = nsml * (bsn + lths) + max(i__7,i__8), i__6 = nx * (bsn + 
		    lths) + nx + max(i__9,i__10);
	    i__3 = 4, i__4 = nsml + max(i__5,i__6);
	    i__1 = jwork, i__2 = nsml + nx + max(i__3,i__4);
	    jwork = max(i__1,i__2);
	}

	if (*ldwork < jwork) {
	    *info = -21;
	    dwork[1] = (doublereal) jwork;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB03BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize the pointers to system matrices and save the possible */
/*     seed for random numbers generation. */

    z__ = 1;
    ac = z__ + nsml;
    dcopy_(&c__4, &dwork[1], &c__1, seed, &c__1);

    wrkopt = 1;

    if (init1) {

/*        Initialize the linear part. */
/*        If N < 0, the order of the system is determined by IB01AD; */
/*        otherwise, the given order will be used. */
/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace:  need:   2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR; */
/*                    prefer: larger. */
/*        Integer workspace:  M+L. (If METH = 'N', (M+L)*NOBR.) */

	ns = *n;
	ir = 1;
	isv = (ml << 1) * *nobr;
	ldr = isv;
	if (lsame_("N", "M", (ftnlen)1, (ftnlen)1)) {
/* Computing MAX */
	    i__1 = ldr, i__2 = mno * 3;
	    ldr = max(i__1,i__2);
	}
	isv = ir + ldr * isv;
	jwork = isv + *l * *nobr;

	i__1 = *ldwork - jwork + 1;
	ib01ad_("M", "F", "N", "O", "N", "N", nobr, m, l, nsmp, &u[u_offset], 
		ldu, &y[y_offset], ldy, n, &dwork[ir], &ldr, &dwork[isv], &
		c_b20, &c_b20, &iwork[1], &dwork[jwork], &i__1, &iwarnl, &
		infol, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);

	if (infol != 0) {
	    *info = infol * 100;
	    return 0;
	}
	if (iwarnl != 0) {
	    *iwarn = iwarnl * 100;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	ircnd = 0;
	if (lsame_("M", "N", (ftnlen)1, (ftnlen)1)) {
	    ircnd = 2;
	    dcopy_(&ircnd, &dwork[jwork + 1], &c__1, rcnd, &c__1);
	}

	if (ns >= 0) {
	    *n = ns;
	} else {

/*           Find the number of parameters. */

	    ldac = *n + *l;
	    isad = ldac * (*n + *m);
	    n2 = *n * *n;
	    lths = *n * (ml + 1) + *l * *m;
	    nx = nths + lths;

	    if (*lx < nx) {
		*lx = nx;
		*info = -16;
		i__1 = -(*info);
		xerbla_("IB03BD", &i__1, (ftnlen)6);
		return 0;
	    }
/*           Workspace for IB01BD. */
/* Computing MAX */
/* Computing MAX */
	    i__3 = lnol * *n + (*n << 1) + (*m + ml) * *nobr + *l, i__4 = (
		    lnol << 1) * *n + n2 + (*n << 3), i__3 = max(i__3,i__4), 
		    i__4 = *n + (mno + *n << 2) + 1, i__3 = max(i__3,i__4), 
		    i__4 = mno + *n * 3 + *l;
	    i__1 = (lnol << 1) * *n + (*n << 1), i__2 = lnol * *n + n2 + *n * 
		    7, i__1 = max(i__1,i__2), i__2 = *l * *nobr * *n + max(
		    i__3,i__4);
	    iw1 = max(i__1,i__2);
	    if (*m > 0) {
/* Computing MAX */
/* Computing 2nd power */
		i__3 = ldac;
		i__1 = i__3 * i__3, i__2 = (*m << 2) * ldac + 1;
		iw2 = *l * *nobr * *n + mno * ldac * (*m * ldac + 1) + max(
			i__1,i__2);
	    } else {
		iw2 = 0;
	    }
	    jwork = isv + isad + max(iw1,iw2);
/*           Workspace for IB01CD. */
/* Computing MAX */
	    i__1 = n2 << 1, i__2 = *n << 2;
	    iw1 = nsml * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
	    i__1 = *n * *l * (*n + 1) + (n2 << 1) + *l * *n, i__2 = *n << 2;
	    iw2 = *n * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
	    i__3 = *n * 5, i__3 = max(i__3,2), i__4 = min(iw1,iw2);
	    i__1 = jwork, i__2 = isad + 2 + *n * (*n + 1 + ldac + *m) + max(
		    i__3,i__4);
	    jwork = max(i__1,i__2);
/*           Workspace for TF01MX. */
/* Computing MAX */
	    i__1 = jwork, i__2 = nsml + isad + ldac + (*n << 1) + *m;
	    jwork = max(i__1,i__2);
/*           Workspace for TB01VD. */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	    i__5 = n2 + *n * max(*n,*l) + *n * 6 + min(*n,*l), i__6 = *n * *m;
	    i__3 = 1, i__4 = n2 * *l + *n * *l + *n, i__3 = max(i__3,i__4), 
		    i__4 = n2 + max(i__5,i__6);
	    i__1 = jwork, i__2 = nsml + isad + *n + max(i__3,i__4);
	    jwork = max(i__1,i__2);
/*           Workspace for MD03BD (whole optimization). */
	    if (*m > 0) {
		iw1 = ldac + *m;
	    } else {
		iw1 = *l;
	    }
/* Computing MAX */
/* Computing MAX */
	    i__3 = *n * ldac;
	    i__1 = *nn << 1, i__2 = isad + (*n << 1) + max(i__3,iw1);
	    iw1 = nsml + max(i__1,i__2);
	    if (*l <= 1 || bsn == 0) {
		iw3 = nx << 2;
		iw2 = iw3 + 1;
	    } else {
/* Computing MAX */
		i__1 = bsn * 3 + 1;
		iw2 = bsn + max(i__1,lths);
		if (*nsmp > bsn) {
/* Computing MAX */
		    i__1 = iw2, i__2 = (lths << 2) + 1;
		    iw2 = max(i__1,i__2);
		    if (*nsmp < bsn << 1) {
/* Computing MAX */
			i__1 = iw2, i__2 = (*nsmp - bsn) * (*l - 1);
			iw2 = max(i__1,i__2);
		    }
		}
		iw3 = lths * bsn + (nx << 1) + (max(bsn,lths) << 1);
	    }
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
	    i__7 = nsml + iw1, i__8 = iw2 + nx;
/* Computing MAX */
	    i__9 = nsml + iw1, i__10 = nx + iw3;
	    i__5 = nsml * (bsn + lths) + max(i__7,i__8), i__6 = nx * (bsn + 
		    lths) + nx + max(i__9,i__10);
	    i__3 = 4, i__4 = nsml + max(i__5,i__6);
	    i__1 = jwork, i__2 = nsml + nx + max(i__3,i__4);
	    jwork = max(i__1,i__2);
	    if (*ldwork < jwork) {
		*info = -21;
		dwork[1] = (doublereal) jwork;
		i__1 = -(*info);
		xerbla_("IB03BD", &i__1, (ftnlen)6);
		return 0;
	    }
	}

	bd = ac + ldac * *n;
	ix = bd + ldac * *m;
	ia = isv;
	ib = ia + ldac * *n;
	iq = ib + ldac * *m;
	if (lsame_("N", "N", (ftnlen)1, (ftnlen)1)) {
	    iry = iq;
	    is = iq;
	    ik = iq;
	    jwork = iq;
	} else {
	    iry = iq + n2;
	    is = iry + *l * *l;
	    ik = is + *n * *l;
	    jwork = ik + *n * *l;
	}

/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace: */
/*          need:  4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) + */
/*                 max( LDW1,LDW2 ), where, */
/*                 LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N, */
/*                              L*NOBR*N + */
/*                              max( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L, */
/*                                   2*(L*NOBR-L)*N+N*N+8*N, */
/*                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) ) */
/*                 LDW2 >= 0,                                  if M = 0; */
/*                 LDW2 >= L*NOBR*N+M*NOBR*(N+L)*(M*(N+L)+1)+ */
/*                         max( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0; */
/*          prefer: larger. */
/*        Integer workspace:  MAX(M*NOBR+N,M*(N+L)). */

	i__1 = *ldwork - jwork + 1;
	ib01bd_("C", "A", "N", nobr, n, m, l, nsmp, &dwork[ir], &ldr, &dwork[
		ia], &ldac, &dwork[ia + *n], &ldac, &dwork[ib], &ldac, &dwork[
		ib + *n], &ldac, &dwork[iq], n, &dwork[iry], l, &dwork[is], n,
		 &dwork[ik], n, &c_b20, &iwork[1], &dwork[jwork], &i__1, 
		bwork, &iwarnl, &infol, (ftnlen)1, (ftnlen)1, (ftnlen)1);

	if (infol == -30) {
	    *info = -21;
	    dwork[1] = dwork[jwork];
	    i__1 = -(*info);
	    xerbla_("IB03BD", &i__1, (ftnlen)6);
	    return 0;
	}
	if (infol != 0) {
	    *info = infol * 100;
	    return 0;
	}
	if (iwarnl != 0) {
	    *iwarn = iwarnl * 100;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	ircndb = 4;
	if (lsame_("N", "K", (ftnlen)1, (ftnlen)1)) {
	    ircndb += 8;
	}
	dcopy_(&ircndb, &dwork[jwork + 1], &c__1, &rcnd[ircnd], &c__1);
	ircnd += ircndb;

/*        Copy the system matrices to the beginning of DWORK, to save */
/*        space, and redefine the pointers. */

	dcopy_(&isad, &dwork[ia], &c__1, &dwork[1], &c__1);
	ia = 1;
	ib = ia + ldac * *n;
	ix0 = ib + ldac * *m;
	iv = ix0 + *n;

/*        Compute the initial condition of the system. On normal exit, */
/*           DWORK(i), i = JWORK+2:JWORK+1+N*N, */
/*           DWORK(j), j = JWORK+2+N*N:JWORK+1+N*N+L*N,  and */
/*           DWORK(k), k = JWORK+2+N*N+L*N:JWORK+1+N*N+L*N+N*M, */
/*        contain the transformed system matrices  At, Ct, and Bt, */
/*        respectively, corresponding to the real Schur form of the */
/*        estimated system state matrix  A. The transformation matrix is */
/*        stored in DWORK(IV:IV+N*N-1). */
/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace: */
/*          need:   (N+L)*(N+M) + N + N*N + 2 + N*( N + M + L ) + */
/*                  max( 5*N, 2, min( LDW1, LDW2 ) ), where, */
/*                  LDW1 = NSMP*L*(N + 1) + 2*N + max( 2*N*N, 4*N), */
/*                  LDW2 = N*(N + 1) + 2*N + */
/*                         max( N*L*(N + 1) + 2*N*N + L*N, 4*N); */
/*          prefer: larger. */
/*        Integer workspace:  N. */

	jwork = iv + n2;
	i__1 = *ldwork - jwork + 1;
	ib01cd_("X needed", "U", "D", n, m, l, nsmp, &dwork[ia], &ldac, &
		dwork[ib], &ldac, &dwork[ia + *n], &ldac, &dwork[ib + *n], &
		ldac, &u[u_offset], ldu, &y[y_offset], ldy, &dwork[ix0], &
		dwork[iv], n, &c_b20, &iwork[1], &dwork[jwork], &i__1, &
		iwarnl, &infol, (ftnlen)8, (ftnlen)1, (ftnlen)1);

	if (infol == -26) {
	    *info = -21;
	    dwork[1] = dwork[jwork];
	    i__1 = -(*info);
	    xerbla_("IB03BD", &i__1, (ftnlen)6);
	    return 0;
	}
	if (infol == 1) {
	    infol = 10;
	}
	if (infol != 0) {
	    *info = infol * 100;
	    return 0;
	}
	if (iwarnl != 0) {
	    *iwarn = iwarnl * 100;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
	++ircnd;
	rcnd[ircnd - 1] = dwork[jwork + 1];

/*        Now, save the system matrices and x0 in the final location. */

	if (iv < ac) {
	    i__1 = isad + *n;
	    dcopy_(&i__1, &dwork[ia], &c__1, &dwork[ac], &c__1);
	} else {
	    i__1 = ac;
	    for (j = ac + isad + *n - 1; j >= i__1; --j) {
		dwork[j] = dwork[ia + j - ac];
/* L10: */
	    }
	}

/*        Compute the output of the linear part. */
/*        Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, */
/*                                                              if M > 0; */
/*                          NSMP*L + (N + L)*N + 2*N + L,       if M = 0; */
/*                   prefer larger. */

	jwork = ix + *n;
	dcopy_(n, &dwork[ix], &c__1, &x[nths + 1], &c__1);
	i__1 = *ldwork - jwork + 1;
	tf01mx_(n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &x[nths 
		+ 1], &dwork[z__], nsmp, &dwork[jwork], &i__1, info);

/*        Convert the state-space representation to output normal form. */
/*        Workspace: */
/*          need:   NSMP*L + (N + L)*(N + M) + N + */
/*                  MAX(1, N*N*L + N*L + N, N*N + */
/*                      MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), N*M)); */
/*          prefer: larger. */

	i__1 = *ldwork - jwork + 1;
	tb01vd_("Apply", n, m, l, &dwork[ac], &ldac, &dwork[bd], &ldac, &
		dwork[ac + *n], &ldac, &dwork[bd + *n], &ldac, &dwork[ix], &x[
		nths + 1], &lths, &dwork[jwork], &i__1, &infol, (ftnlen)5);

	if (infol > 0) {
	    *info = infol + 4;
	    return 0;
	}
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

    }

    lipar = 7;
    iw1 = 0;
    iw2 = 0;
    idiag = ac;

    if (init2) {

/*        Initialize the nonlinear part. */

	if (! init1) {
	    bd = ac + ldac * *n;
	    ix = bd + ldac * *m;

/*           Convert the output normal form to state-space model. */
/*           Workspace: need NSMP*L + (N + L)*(2*N + M) + 2*N. */
/*           (NSMP*L locations are reserved for the output of the linear */
/*           part.) */

	    jwork = ix + *n;
	    i__1 = *ldwork - jwork + 1;
	    tb01vy_("Apply", n, m, l, &x[nths + 1], &lths, &dwork[ac], &ldac, 
		    &dwork[bd], &ldac, &dwork[ac + *n], &ldac, &dwork[bd + *n]
		    , &ldac, &dwork[ix], &dwork[jwork], &i__1, info, (ftnlen)
		    5);

/*           Compute the output of the linear part. */
/*           Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, */
/*                                                              if M > 0; */
/*                             NSMP*L + (N + L)*N + 2*N + L,    if M = 0; */
/*                      prefer larger. */

	    i__1 = *ldwork - jwork + 1;
	    tf01mx_(n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &
		    dwork[ix], &dwork[z__], nsmp, &dwork[jwork], &i__1, info);
	}

/*        Optimize the parameters of the nonlinear part. */
/*        Workspace: */
/*          need   NSMP*L + BSN + */
/*                 MAX( 4, NSMP + */
/*                      MAX( NSMP*BSN + MAX( 2*NN, 5*BSN + 1 ), */
/*                           BSN**2 + BSN + MAX( NSMP + 2*NN, 5*BSN ) )); */
/*          prefer larger. */
/*        Integer workspace:  NN*(L + 2) + 2. */

	work[0] = 0.;
	dcopy_(&c__3, work, &c__0, &work[1], &c__1);

/*        Set the integer parameters needed, including the number of */
/*        neurons. */

	ipar[0] = *nsmp;
	ipar[1] = *l;
	ipar[2] = *nn;
	jwork = idiag + bsn;

	i__1 = *l - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    dcopy_(&c__4, seed, &c__1, &dwork[jwork], &c__1);
	    i__2 = *ldwork - jwork + 1;
	    md03bd_("Random initialization", "I", "E", (U_fp)nf01be_, (U_fp)
		    md03ba_, (U_fp)md03bb_, nsmp, &bsn, itmax1, &c_b63, 
		    nprint, ipar, &lipar, &dwork[z__], nsmp, &y[(i__ + 1) * 
		    y_dim1 + 1], ldy, &x[i__ * bsn + 1], &dwork[idiag], &nfev,
		     &njev, tol1, tol1, &c_b64, &c_b64, &iwork[1], &dwork[
		    jwork], &i__2, &iwarnl, &infol, (ftnlen)21, (ftnlen)1, (
		    ftnlen)1);
	    if (infol != 0) {
		*info = infol * 10;
		return 0;
	    }
	    if (iwarnl < 0) {
		*info = infol;
		*iwarn = iwarnl;
		goto L50;
	    } else if (iwarnl > 0) {
		if (*iwarn > 100) {
/* Computing MAX */
		    i__2 = *iwarn, i__3 = *iwarn / 100 * 100 + iwarnl * 10;
		    *iwarn = max(i__2,i__3);
		} else {
/* Computing MAX */
		    i__2 = *iwarn, i__3 = iwarnl * 10;
		    *iwarn = max(i__2,i__3);
		}
	    }
/* Computing MAX */
	    d__1 = work[0], d__2 = dwork[jwork];
	    work[0] = max(d__1,d__2);
/* Computing MAX */
	    d__1 = work[1], d__2 = dwork[jwork + 1];
	    work[1] = max(d__1,d__2);
/* Computing MAX */
	    d__1 = work[3], d__2 = dwork[jwork + 3];
	    work[3] = max(d__1,d__2);
	    work[2] += dwork[jwork + 2];
	    iw1 = nfev + iw1;
	    iw2 = njev + iw2;
/* L30: */
	}

    }

/*     Main iteration. */
/*     Workspace: */
/*       need   NSMP*L + NX + */
/*              MAX( 4, NSMP*L + */
/*                      MAX( NSMP*L*( BSN + LTHS ) + */
/*                           MAX( NSMP*L + LDW1, LDW2 + NX ), */
/*                                NX*( BSN + LTHS ) + NX + */
/*                                MAX( NSMP*L + LDW1, NX + LDW3 ) ) ), */
/*              LDW0 = MAX( N*(N+L), N+M+L ),    if M > 0; */
/*              LDW0 = MAX( N*(N+L), L ),        if M = 0; */
/*              LDW1 = NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + LDW0); */
/*              LDW2 = 4*NX + 1,  if L <= 1 or BSN = 0; otherwise, */
/*              LDW2 = BSN + MAX(3*BSN+1,LTHS); */
/*              LDW2 = MAX(LDW2, 4*LTHS+1),        if NSMP > BSN; */
/*              LDW2 = MAX(LDW2, (NSMP-BSN)*(L-1)), if BSN < NSMP < 2*BSN; */
/*              LDW3 = 4*NX,                       if L <= 1 or BSN = 0; */
/*              LDW3 = LTHS*BSN + 2*NX + 2*MAX(BSN,LTHS), */
/*                                                 if L > 1 and BSN > 0; */
/*       prefer larger. */
/*     Integer workspace:  NX+L. */

/*     Set the integer parameters describing the Jacobian structure */
/*     and the number of neurons. */

    ipar[0] = lths;
    ipar[1] = *l;
    ipar[2] = *nsmp;
    ipar[3] = bsn;
    ipar[4] = *m;
    ipar[5] = *n;
    ipar[6] = *nn;
    jwork = idiag + nx;

    i__1 = *ldwork - jwork + 1;
    md03bd_("Given initialization", "I", "E", (U_fp)nf01bf_, (U_fp)nf01bs_, (
	    U_fp)nf01bp_, &nsml, &nx, itmax2, &c_b63, nprint, ipar, &lipar, &
	    u[u_offset], ldu, &y[y_offset], ldy, &x[1], &dwork[idiag], &nfev, 
	    &njev, tol2, tol2, &c_b64, &c_b64, &iwork[1], &dwork[jwork], &
	    i__1, &iwarnl, info, (ftnlen)20, (ftnlen)1, (ftnlen)1);
    if (*info != 0) {
	return 0;
    }

    i__1 = nx + *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__ + 3] = iwork[i__];
/* L40: */
    }

L50:
    iwork[1] = iw1 + nfev;
    iwork[2] = iw2 + njev;
    if (iwarnl < 0) {
	*iwarn = iwarnl;
    } else {
	*iwarn += iwarnl;
    }
    dcopy_(&c__4, &dwork[jwork], &c__1, &dwork[1], &c__1);
    if (init2) {
	dcopy_(&c__4, work, &c__1, &dwork[5], &c__1);
    }
    if (init1) {
	iwork[3] = ircnd;
	dcopy_(&ircnd, rcnd, &c__1, &dwork[9], &c__1);
    } else {
	iwork[3] = 0;
    }

    return 0;

/* *** Last line of IB03BD *** */
} /* ib03bd_ */

