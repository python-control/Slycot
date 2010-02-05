      SUBROUTINE IB03BD( INIT, NOBR, M, L, NSMP, N, NN, ITMAX1, ITMAX2,
     $                   NPRINT, U, LDU, Y, LDY, X, LX, TOL1, TOL2,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To compute a set of parameters for approximating a Wiener system
C     in a least-squares sense, using a neural network approach and a
C     MINPACK-like Levenberg-Marquardt algorithm. The Wiener system
C     consists of a linear part and a static nonlinearity, and it is
C     represented as
C
C        x(t+1) = A*x(t) + B*u(t)
C        z(t)   = C*x(t) + D*u(t),
C
C        y(t)   = f(z(t),wb(1:L)),
C
C     where t = 1, 2, ..., NSMP, and f is a nonlinear function,
C     evaluated by the SLICOT Library routine NF01AY. The parameter
C     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ),
C     where theta corresponds to the linear part, and wb(i), i = 1 : L,
C     correspond to the nonlinear part. See SLICOT Library routine
C     NF01AD for further details.
C
C     The sum of squares of the error functions, defined by
C
C        e(t) = y(t) - Y(t),  t = 1, 2, ..., NSMP,
C
C     is minimized, where Y(t) is the measured output vector. The
C     functions and their Jacobian matrices are evaluated by SLICOT
C     Library routine NF01BF (the FCN routine in the call of MD03BD).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     INIT    CHARACTER*1
C             Specifies which parts have to be initialized, as follows:
C             = 'L' : initialize the linear part only, X already
C                     contains an initial approximation of the
C                     nonlinearity;
C             = 'S' : initialize the static nonlinearity only, X
C                     already contains an initial approximation of the
C                     linear part;
C             = 'B' : initialize both linear and nonlinear parts;
C             = 'N' : do not initialize anything, X already contains
C                     an initial approximation.
C             If INIT = 'S' or 'B', the error functions for the
C             nonlinear part, and their Jacobian matrices, are evaluated
C             by SLICOT Library routine NF01BE (used as a second FCN
C             routine in the MD03BD call for the initialization step,
C             see METHOD).
C
C     Input/Output Parameters
C
C     NOBR    (input) INTEGER
C             If INIT = 'L' or 'B', NOBR is the number of block rows, s,
C             in the input and output block Hankel matrices to be
C             processed for estimating the linear part.  NOBR > 0.
C             (In the MOESP theory,  NOBR  should be larger than  n,
C             the estimated dimension of state vector.)
C             This parameter is ignored if INIT is 'S' or 'N'.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     L       (input) INTEGER
C             The number of system outputs.  L >= 0, and L > 0, if
C             INIT = 'L' or 'B'.
C
C     NSMP    (input) INTEGER
C             The number of input and output samples, t.  NSMP >= 0, and
C             NSMP >= 2*(M+L+1)*NOBR - 1, if INIT = 'L' or 'B'.
C
C     N       (input/output) INTEGER
C             The order of the linear part.
C             If INIT = 'L' or 'B', and N < 0 on entry, the order is
C             assumed unknown and it will be found by the routine.
C             Otherwise, the input value will be used. If INIT = 'S'
C             or 'N', N must be non-negative. The values N >= NOBR,
C             or N = 0, are not acceptable if INIT = 'L' or 'B'.
C
C     NN      (input) INTEGER
C             The number of neurons which shall be used to approximate
C             the nonlinear part.  NN >= 0.
C
C     ITMAX1  (input) INTEGER
C             The maximum number of iterations for the initialization of
C             the static nonlinearity.
C             This parameter is ignored if INIT is 'N' or 'L'.
C             Otherwise, ITMAX1 >= 0.
C
C     ITMAX2  (input) INTEGER
C             The maximum number of iterations.  ITMAX2 >= 0.
C
C     NPRINT  (input) INTEGER
C             This parameter enables controlled printing of iterates if
C             it is positive. In this case, FCN is called with IFLAG = 0
C             at the beginning of the first iteration and every NPRINT
C             iterations thereafter and immediately prior to return,
C             and the current error norm is printed. Other intermediate
C             results could be printed by modifying the corresponding
C             FCN routine (NF01BE and/or NF01BF). If NPRINT <= 0, no
C             special calls of FCN with IFLAG = 0 are made.
C
C     U       (input) DOUBLE PRECISION array, dimension (LDU, M)
C             The leading NSMP-by-M part of this array must contain the
C             set of input samples,
C             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ).
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,NSMP).
C
C     Y       (input) DOUBLE PRECISION array, dimension (LDY, L)
C             The leading NSMP-by-L part of this array must contain the
C             set of output samples,
C             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ).
C
C     LDY     INTEGER
C             The leading dimension of array Y.  LDY >= MAX(1,NSMP).
C
C     X       (input/output) DOUBLE PRECISION array dimension (LX)
C             On entry, if INIT = 'L', the leading (NN*(L+2) + 1)*L part
C             of this array must contain the initial parameters for
C             the nonlinear part of the system.
C             On entry, if INIT = 'S', the elements lin1 : lin2 of this
C             array must contain the initial parameters for the linear
C             part of the system, corresponding to the output normal
C             form, computed by SLICOT Library routine TB01VD, where
C                lin1 = (NN*(L+2) + 1)*L + 1;
C                lin2 = (NN*(L+2) + 1)*L + N*(L+M+1) + L*M.
C             On entry, if INIT = 'N', the elements 1 : lin2 of this
C             array must contain the initial parameters for the
C             nonlinear part followed by the initial parameters for the
C             linear part of the system, as specified above.
C             This array need not be set on entry if INIT = 'B'.
C             On exit, the elements 1 : lin2 of this array contain the
C             optimal parameters for the nonlinear part followed by the
C             optimal parameters for the linear part of the system, as
C             specified above.
C
C     LX      (input/output) INTEGER
C             On entry, this parameter must contain the intended length
C             of X. If N >= 0, then LX >= NX := lin2 (see parameter X).
C             If N is unknown (N < 0 on entry), a large enough estimate
C             of N should be used in the formula of lin2.
C             On exit, if N < 0 on entry, but LX is not large enough,
C             then this parameter contains the actual length of X,
C             corresponding to the computed N. Otherwise, its value
C             is unchanged.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If INIT = 'S' or 'B' and TOL1 >= 0, TOL1 is the tolerance
C             which measures the relative error desired in the sum of
C             squares, as well as the relative error desired in the
C             approximate solution, for the initialization step of
C             nonlinear part. Termination occurs when either both the
C             actual and predicted relative reductions in the sum of
C             squares, or the relative error between two consecutive
C             iterates are at most TOL1. If the user sets  TOL1 < 0,
C             then  SQRT(EPS)  is used instead TOL1, where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH).
C             This parameter is ignored if INIT is 'N' or 'L'.
C
C     TOL2    DOUBLE PRECISION
C             If TOL2 >= 0, TOL2 is the tolerance which measures the
C             relative error desired in the sum of squares, as well as
C             the relative error desired in the approximate solution,
C             for the whole optimization process. Termination occurs
C             when either both the actual and predicted relative
C             reductions in the sum of squares, or the relative error
C             between two consecutive iterates are at most TOL2. If the
C             user sets TOL2 < 0, then  SQRT(EPS)  is used instead TOL2.
C             This default value could require many iterations,
C             especially if TOL1 is larger. If INIT = 'S' or 'B', it is
C             advisable that TOL2 be larger than TOL1, and spend more
C             time with cheaper iterations.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX( LIW1, LIW2, LIW3 )), where
C             LIW1 = LIW2 = 0,  if INIT = 'S' or 'N'; otherwise,
C             LIW1 = M+L;
C             LIW2 = MAX(M*NOBR+N,M*(N+L));
C             LIW3 = 3+MAX(NN*(L+2)+2,NX+L), if INIT = 'S' or 'B';
C             LIW3 = 3+NX+L,                 if INIT = 'L' or 'N'.
C             On output, if INFO = 0, IWORK(1) and IWORK(2) return the
C             (total) number of function and Jacobian evaluations,
C             respectively (including the initialization step, if it was
C             performed), and if INIT = 'L' or INIT = 'B', IWORK(3)
C             specifies how many locations of DWORK contain reciprocal
C             condition number estimates (see below); otherwise,
C             IWORK(3) = 0. If INFO = 0, the entries 4 to 3+NX of IWORK
C             define a permutation matrix P such that J*P = Q*R, where
C             J is the final calculated Jacobian, Q is an orthogonal
C             matrix (not stored), and R is upper triangular with
C             diagonal elements of nonincreasing magnitude (possibly
C             for each block column of J). Column j of P is column
C             IWORK(3+j) of the identity matrix. Moreover, the entries
C             4+NX:3+NX+L of this array contain the ranks of the final
C             submatrices S_k (see description of LMPARM in MD03BD).
C
C     DWORK   DOUBLE PRECISION array dimesion (LDWORK)
C             On entry, if desired, and if INIT = 'S' or 'B', the
C             entries DWORK(1:4) are set to initialize the random
C             numbers generator for the nonlinear part parameters (see
C             the description of the argument XINIT of SLICOT Library
C             routine MD03BD); this enables to obtain reproducible
C             results. The same seed is used for all outputs.
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, DWORK(2) returns the residual error norm (the
C             sum of squares), DWORK(3) returns the number of iterations
C             performed, and DWORK(4) returns the final Levenberg
C             factor, for optimizing the parameters of both the linear
C             part and the static nonlinearity part. If INIT = 'S' or
C             INIT = 'B' and INFO = 0, then the elements DWORK(5) to
C             DWORK(8) contain the corresponding four values for the
C             initialization step (see METHOD). (If L > 1, DWORK(8)
C             contains the maximum of the Levenberg factors for all
C             outputs.) If INIT = 'L' or INIT = 'B', and INFO = 0,
C             DWORK(9) to DWORK(8+IWORK(3)) contain reciprocal condition
C             number estimates set by SLICOT Library routines IB01AD,
C             IB01BD, and IB01CD.
C             On exit, if  INFO = -21,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             In the formulas below, N should be taken not larger than
C             NOBR - 1, if N < 0 on entry.
C             LDWORK = MAX( LW1, LW2, LW3, LW4 ), where
C             LW1 = 0, if INIT = 'S' or 'N'; otherwise,
C             LW1 = MAX( 2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR,
C                        4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) +
C                        MAX( LDW1, LDW2 ),
C                        (N+L)*(N+M) + N + N*N + 2 + N*(N+M+L) +
C                        MAX( 5*N, 2, MIN( LDW3, LDW4 ), LDW5, LDW6 ),
C                 where,
C                 LDW1 >= MAX( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N,
C                              L*NOBR*N +
C                              MAX( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L,
C                                   2*(L*NOBR-L)*N+N*N+8*N,
C                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) )
C                 LDW2 >= 0,                                  if M = 0;
C                 LDW2 >= L*NOBR*N + M*NOBR*(N+L)*(M*(N+L)+1) +
C                         MAX( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0;
C                 LDW3 = NSMP*L*(N+1) + 2*N + MAX( 2*N*N, 4*N ),
C                 LDW4 = N*(N+1) + 2*N +
C                        MAX( N*L*(N+1) + 2*N*N + L*N, 4*N );
C                 LDW5 = NSMP*L + (N+L)*(N+M) + 3*N+M+L;
C                 LDW6 = NSMP*L + (N+L)*(N+M) + N +
C                        MAX(1, N*N*L + N*L + N, N*N +
C                            MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L),
C                                N*M));
C             LW2 = LW3 = 0, if INIT = 'L' or 'N'; otherwise,
C             LW2 = NSMP*L + BSN +
C                   MAX( 4, NSMP +
C                           MAX( NSMP*BSN + MAX( 2*NN, 5*BSN + 1 ),
C                                BSN**2 + BSN +
C                                MAX( NSMP + 2*NN, 5*BSN ) ) );
C             LW3 = MAX( LDW7, NSMP*L + (N+L)*(2*N+M) + 2*N );
C                 LDW7 = NSMP*L + (N+L)*(N+M) + 3*N+M+L,  if M > 0;
C                 LDW7 = NSMP*L + (N+L)*N + 2*N+L,        if M = 0;
C             LW4 = NSMP*L + NX +
C                   MAX( 4, NSMP*L +
C                           MAX( NSMP*L*( BSN + LTHS ) +
C                                MAX( NSMP*L + L1, L2 + NX ),
C                                     NX*( BSN + LTHS ) + NX +
C                                     MAX( NSMP*L + L1, NX + L3 ) ) ),
C                  L0 = MAX( N*(N+L), N+M+L ),    if M > 0;
C                  L0 = MAX( N*(N+L), L ),        if M = 0;
C                  L1 = NSMP*L + MAX( 2*NN, (N+L)*(N+M) + 2*N + L0);
C                  L2 = 4*NX + 1,  if L <= 1 or BSN = 0; otherwise,
C                  L2 = BSN + MAX(3*BSN+1,LTHS);
C                  L2 = MAX(L2,4*LTHS+1),         if NSMP > BSN;
C                  L2 = MAX(L2,(NSMP-BSN)*(L-1)), if BSN < NSMP < 2*BSN;
C                  L3 = 4*NX,                     if L <= 1 or BSN = 0;
C                  L3 = LTHS*BSN + 2*NX + 2*MAX(BSN,LTHS),
C                                                 if L > 1 and BSN > 0,
C                  with BSN  = NN*( L + 2 ) + 1,
C                       LTHS = N*( L + M + 1 ) + L*M.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             < 0:  the user set IFLAG = IWARN in (one of) the
C                   subroutine(s) FCN, i.e., NF01BE, if INIT = 'S'
C                   or 'B', and/or NF01BF; this value cannot be returned
C                   without changing the FCN routine(s);
C                   otherwise, IWARN has the value k*100 + j*10 + i,
C                   where k is defined below, i refers to the whole
C                   optimization process, and j refers to the
C                   initialization step (j = 0, if INIT = 'L' or 'N'),
C                   and the possible values for i and j have the
C                   following meaning (where TOL* denotes TOL1 or TOL2,
C                   and similarly for ITMAX*):
C             = 1:  both actual and predicted relative reductions in
C                   the sum of squares are at most TOL*;
C             = 2:  relative error between two consecutive iterates is
C                   at most TOL*;
C             = 3:  conditions for i or j = 1 and i or j = 2 both hold;
C             = 4:  the cosine of the angle between the vector of error
C                   function values and any column of the Jacobian is at
C                   most EPS in absolute value;
C             = 5:  the number of iterations has reached ITMAX* without
C                   satisfying any convergence condition;
C             = 6:  TOL* is too small: no further reduction in the sum
C                   of squares is possible;
C             = 7:  TOL* is too small: no further improvement in the
C                   approximate solution X is possible;
C             = 8:  the vector of function values e is orthogonal to the
C                   columns of the Jacobian to machine precision.
C             The digit k is normally 0, but if INIT = 'L' or 'B', it
C             can have a value in the range 1 to 6 (see IB01AD, IB01BD
C             and IB01CD). In all these cases, the entries DWORK(1:4),
C             DWORK(5:8) (if INIT = 'S' or 'B'), and DWORK(9:8+IWORK(3))
C             (if INIT = 'L' or 'B'), are set as described above.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C                   otherwise, INFO has the value k*100 + j*10 + i,
C                   where k is defined below, i refers to the whole
C                   optimization process, and j refers to the
C                   initialization step (j = 0, if INIT = 'L' or 'N'),
C                   and the possible values for i and j have the
C                   following meaning:
C             = 1:  the routine FCN returned with INFO <> 0 for
C                   IFLAG = 1;
C             = 2:  the routine FCN returned with INFO <> 0 for
C                   IFLAG = 2;
C             = 3:  the routine QRFACT returned with INFO <> 0;
C             = 4:  the routine LMPARM returned with INFO <> 0.
C             In addition, if INIT = 'L' or 'B', i could also be
C             = 5:  if a Lyapunov equation could not be solved;
C             = 6:  if the identified linear system is unstable;
C             = 7:  if the QR algorithm failed on the state matrix
C                   of the identified linear system.
C             QRFACT and LMPARM are generic names for SLICOT Library
C             routines NF01BS and NF01BP, respectively, for the whole
C             optimization process, and MD03BA and MD03BB, respectively,
C             for the initialization step (if INIT = 'S' or 'B').
C             The digit k is normally 0, but if INIT = 'L' or 'B', it
C             can have a value in the range 1 to 10 (see IB01AD/IB01BD).
C
C     METHOD
C
C     If INIT = 'L' or 'B', the linear part of the system is
C     approximated using the combined MOESP and N4SID algorithm. If
C     necessary, this algorithm can also choose the order, but it is
C     advantageous if the order is already known.
C
C     If INIT = 'S' or 'B', the output of the approximated linear part
C     is computed and used to calculate an approximation of the static
C     nonlinearity using the Levenberg-Marquardt algorithm [1,3].
C     This step is referred to as the (nonlinear) initialization step.
C
C     As last step, the Levenberg-Marquardt algorithm is used again to
C     optimize the parameters of the linear part and the static
C     nonlinearity as a whole. Therefore, it is necessary to parametrise
C     the matrices of the linear part. The output normal form [2]
C     parameterisation is used.
C
C     The Jacobian is computed analytically, for the nonlinear part, and
C     numerically, for the linear part.
C
C     REFERENCES
C
C     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E.
C         User's Guide for MINPACK-1.
C         Applied Math. Division, Argonne National Laboratory, Argonne,
C         Illinois, Report ANL-80-74, 1980.
C
C     [2] Peeters, R.L.M., Hanzon, B., and Olivi, M.
C         Balanced realizations of discrete-time stable all-pass
C         systems and the tangential Schur algorithm.
C         Proceedings of the European Control Conference,
C         31 August - 3 September 1999, Karlsruhe, Germany.
C         Session CP-6, Discrete-time Systems, 1999.
C
C     [3] More, J.J.
C         The Levenberg-Marquardt algorithm: implementation and theory.
C         In Watson, G.A. (Ed.), Numerical Analysis, Lecture Notes in
C         Mathematics, vol. 630, Springer-Verlag, Berlin, Heidelberg
C         and New York, pp. 105-116, 1978.
C
C     NUMERICAL ASPECTS
C
C     The Levenberg-Marquardt algorithm described in [3] is scaling
C     invariant and globally convergent to (maybe local) minima.
C     The convergence rate near a local minimum is quadratic, if the
C     Jacobian is computed analytically, and linear, if the Jacobian
C     is computed numerically.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, March, 2002, Apr. 2002, Feb. 2004, March 2005.
C
C     KEYWORDS
C
C     Least-squares approximation, Levenberg-Marquardt algorithm,
C     matrix operations, optimization.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     FACTOR is a scaling factor for variables (see MD03BD).
      DOUBLE PRECISION  FACTOR
      PARAMETER         ( FACTOR = 100.0D0 )
C     Condition estimation and internal scaling of variables are used
C     (see MD03BD).
      CHARACTER         COND, SCALE
      PARAMETER         ( COND = 'E', SCALE = 'I' )
C     Default tolerances are used in MD03BD for measuring the
C     orthogonality between the vector of function values and columns
C     of the Jacobian (GTOL), and for the rank estimations (TOL).
      DOUBLE PRECISION  GTOL, TOL
      PARAMETER         ( GTOL = 0.0D0, TOL = 0.0D0 )
C     For INIT = 'L' or 'B', additional parameters are set:
C     The following six parameters are used in the call of IB01AD;
      CHARACTER         ALG, BATCH, CONCT, CTRL, JOBD, METH
      PARAMETER         ( ALG   = 'Fast QR',     BATCH = 'One batch',
     $                    CONCT = 'Not connect', CTRL  = 'Not confirm',
     $                    JOBD  = 'Not MOESP',   METH  = 'MOESP' )
C     The following three parameters are used in the call of IB01BD;
      CHARACTER         JOB, JOBCK, METHB
      PARAMETER         ( JOB   = 'All matrices',
     $                    JOBCK = 'No Kalman gain',
     $                    METHB = 'Combined MOESP+N4SID' )
C     The following two parameters are used in the call of IB01CD;
      CHARACTER         COMUSE, JOBXD
      PARAMETER         ( COMUSE = 'Use B, D',
     $                    JOBXD  = 'D also' )
C     TOLN controls the estimated order in IB01AD (default value);
      DOUBLE PRECISION  TOLN
      PARAMETER         ( TOLN = -1.0D0 )
C     RCOND controls the rank decisions in IB01AD, IB01BD, and IB01CD
C     (default);
      DOUBLE PRECISION  RCOND
      PARAMETER         ( RCOND = -1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         INIT
      INTEGER           INFO, ITMAX1, ITMAX2, IWARN, L, LDU, LDWORK,
     $                  LDY, LX, M, N, NN, NOBR, NPRINT, NSMP
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), U(LDU, *), X(*), Y(LDY, *)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      INTEGER           AC, BD, BSN, I, IA, IB, IDIAG, IK, INFOL, IQ,
     $                  IR, IRCND, IRCNDB, IRY, IS, ISAD, ISV, IV, IW1,
     $                  IW2, IW3, IWARNL, IX, IX0, J, JWORK, LDAC, LDR,
     $                  LIPAR, LNOL, LTHS, ML, MNO, N2, NFEV, NJEV, NS,
     $                  NSML, NTHS, NX, WRKOPT, Z
      LOGICAL           INIT1, INIT2
C     .. Local Arrays ..
      LOGICAL           BWORK(1)
      INTEGER           IPAR(7)
      DOUBLE PRECISION  RCND(16), SEED(4), WORK(4)
C     .. External Functions ..
      EXTERNAL          LSAME
      LOGICAL           LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY,  IB01AD, IB01BD, IB01CD, MD03BA, MD03BB,
     $                  MD03BD, NF01BE, NF01BF, NF01BP, NF01BS, TB01VD,
     $                  TB01VY, TF01MX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INIT1 = LSAME( INIT, 'B' ) .OR. LSAME( INIT, 'L' )
      INIT2 = LSAME( INIT, 'B' ) .OR. LSAME( INIT, 'S' )
C
      ML    = M + L
      INFO  = 0
      IWARN = 0
      IF ( .NOT.( INIT1 .OR. INIT2 .OR. LSAME( INIT, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( INIT1 .AND. NOBR.LE.0 ) THEN
         INFO = -2
      ELSEIF ( M.LT.0 ) THEN
         INFO = -3
      ELSEIF ( L.LT.0 .OR. ( INIT1 .AND. L.EQ.0 ) ) THEN
         INFO = -4
      ELSEIF ( NSMP.LT.0 .OR.
     $         ( INIT1 .AND. NSMP.LT.2*( ML + 1 )*NOBR - 1 ) ) THEN
         INFO = -5
      ELSEIF ( ( N.LT.0 .AND. .NOT.INIT1 ) .OR.
     $       ( ( N.EQ.0 .OR. N.GE.NOBR ) .AND. INIT1 ) ) THEN
         INFO = -6
      ELSEIF ( NN.LT.0 ) THEN
         INFO = -7
      ELSEIF ( INIT2 .AND. ( ITMAX1.LT.0 ) ) THEN
         INFO = -8
      ELSEIF ( ITMAX2.LT.0 ) THEN
         INFO = -9
      ELSEIF ( LDU.LT.MAX( 1, NSMP ) ) THEN
         INFO = -12
      ELSEIF ( LDY.LT.MAX( 1, NSMP ) ) THEN
         INFO = -14
      ELSE
         LNOL = L*NOBR - L
         MNO  = M*NOBR
         BSN  = NN*( L + 2 ) + 1
         NTHS =  BSN*L
         NSML = NSMP*L
         IF ( N.GT.0 ) THEN
            LDAC = N + L
            ISAD = LDAC*( N + M )
            N2   = N*N
         END IF
C
C        Check the workspace size.
C
         JWORK = 0
         IF ( INIT1 ) THEN
C           Workspace for IB01AD.
            JWORK = 2*ML*NOBR*( 2*ML*( NOBR + 1 ) + 3 ) + L*NOBR
            IF ( N.GT.0 ) THEN
C              Workspace for IB01BD.
               IW1 = MAX( 2*LNOL*N + 2*N, LNOL*N + N2 + 7*N, L*NOBR*N +
     $                    MAX( LNOL*N + 2*N + ( M + ML )*NOBR + L,
     $                         2*LNOL*N + N2 + 8*N, N + 4*( MNO + N ) +
     $                         1, MNO + 3*N + L ) )
               IF ( M.GT.0 ) THEN
                  IW2 = L*NOBR*N + MNO*LDAC*( M*LDAC + 1 ) +
     $                  MAX( LDAC**2, 4*M*LDAC + 1 )
               ELSE
                  IW2 = 0
               END IF
               JWORK = MAX( JWORK,
     $                      ( 2*ML*NOBR )**2 + ISAD + MAX( IW1, IW2 ) )
C              Workspace for IB01CD.
               IW1   = NSML*( N + 1 ) + 2*N + MAX( 2*N2, 4*N )
               IW2   = N*( N + 1 ) + 2*N +
     $                 MAX( N*L*( N + 1 ) + 2*N2 + L*N, 4*N )
               JWORK = MAX( JWORK, ISAD + 2 + N*( N + 1 + LDAC + M ) +
     $                     MAX( 5*N, 2, MIN( IW1, IW2 ) ) )
C              Workspace for TF01MX.
               JWORK = MAX( JWORK, NSML + ISAD + LDAC + 2*N + M )
C              Workspace for TB01VD.
               JWORK = MAX( JWORK, NSML + ISAD + N +
     $                      MAX( 1, N2*L + N*L + N,
     $                           N2 + MAX( N2 + N*MAX( N, L ) +
     $                                     6*N +  MIN( N, L ), N*M ) ) )
            END IF
         END IF
C
         IF ( INIT2 ) THEN
C           Workspace for MD03BD (initialization of the nonlinear part).
            JWORK = MAX( JWORK, NSML + BSN +
     $                   MAX( 4, NSMP +
     $                        MAX( NSMP*BSN + MAX( 2*NN, 5*BSN + 1 ),
     $                             BSN**2 + BSN +
     $                             MAX( NSMP + 2*NN, 5*BSN ) ) ) )
            IF ( N.GT.0 .AND. .NOT.INIT1 ) THEN
C              Workspace for TB01VY.
               JWORK = MAX( JWORK, NSML + LDAC*( 2*N + M ) + 2*N )
C              Workspace for TF01MX.
               IF ( M.GT.0 ) THEN
                  IW1 = N + M
               ELSE
                  IW1 = 0
               END IF
               JWORK = MAX( JWORK, NSML + ISAD + IW1 + LDAC + N )
            END IF
         END IF
C
         IF ( N.GE.0 ) THEN
C
C           Find the number of parameters.
C
            LTHS = N*( ML + 1 ) + L*M
            NX   = NTHS + LTHS
C
            IF ( LX.LT.NX ) THEN
               INFO = -16
               CALL XERBLA( 'IB03BD', -INFO )
               RETURN
            END IF
C
C           Workspace for MD03BD (whole optimization).
C
            IF ( M.GT.0 ) THEN
               IW1 = LDAC + M
            ELSE
               IW1 = L
            END IF
            IW1 = NSML + MAX( 2*NN, ISAD + 2*N + MAX( N*LDAC, IW1 ) )
            IF ( L.LE.1 .OR. BSN.EQ.0 ) THEN
               IW3 = 4*NX
               IW2 = IW3 + 1
            ELSE
               IW2 = BSN + MAX( 3*BSN + 1, LTHS )
               IF ( NSMP.GT.BSN ) THEN
                  IW2 = MAX( IW2, 4*LTHS + 1 )
                  IF ( NSMP.LT.2*BSN )
     $               IW2 = MAX( IW2, ( NSMP - BSN )*( L - 1 ) )
               END IF
               IW3 = LTHS*BSN + 2*NX + 2*MAX( BSN, LTHS )
            END IF
            JWORK = MAX( JWORK, NSML + NX +
     $                   MAX( 4, NSML +
     $                           MAX( NSML*( BSN + LTHS ) +
     $                                MAX( NSML + IW1, IW2 + NX ),
     $                                     NX*( BSN + LTHS ) + NX +
     $                                     MAX( NSML + IW1, NX + IW3 ) )
     $                      ) )
         END IF
C
         IF ( LDWORK.LT.JWORK ) THEN
            INFO = -21
            DWORK(1) = JWORK
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'IB03BD', -INFO )
         RETURN
      END IF
C
C     Initialize the pointers to system matrices and save the possible
C     seed for random numbers generation.
C
      Z  = 1
      AC = Z + NSML
      CALL DCOPY( 4, DWORK, 1, SEED, 1 )
C
      WRKOPT = 1
C
      IF ( INIT1 ) THEN
C
C        Initialize the linear part.
C        If N < 0, the order of the system is determined by IB01AD;
C        otherwise, the given order will be used.
C        The workspace needed is defined for the options set above
C        in the PARAMETER statements.
C        Workspace:  need:   2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR;
C                    prefer: larger.
C        Integer workspace:  M+L. (If METH = 'N', (M+L)*NOBR.)
C
         NS  = N
         IR  = 1
         ISV = 2*ML*NOBR
         LDR = ISV
         IF ( LSAME( JOBD, 'M' ) )
     $      LDR = MAX( LDR, 3*MNO )
         ISV   = IR  + LDR*ISV
         JWORK = ISV + L*NOBR
C
         CALL IB01AD( METH, ALG, JOBD, BATCH, CONCT, CTRL, NOBR, M, L,
     $                NSMP, U, LDU, Y, LDY, N, DWORK(IR), LDR,
     $                DWORK(ISV), RCOND, TOLN, IWORK, DWORK(JWORK),
     $                LDWORK-JWORK+1, IWARNL, INFOL )
C
         IF( INFOL.NE.0 ) THEN
            INFO = 100*INFOL
            RETURN
         END IF
         IF( IWARNL.NE.0 )
     $      IWARN = 100*IWARNL
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IRCND  = 0
         IF ( LSAME( METH, 'N' ) ) THEN
            IRCND = 2
            CALL DCOPY( IRCND, DWORK(JWORK+1), 1, RCND, 1 )
         END IF
C
         IF ( NS.GE.0 ) THEN
            N = NS
         ELSE
C
C           Find the number of parameters.
C
            LDAC = N + L
            ISAD = LDAC*( N + M )
            N2   = N*N
            LTHS = N*( ML + 1 ) + L*M
            NX   = NTHS + LTHS
C
            IF ( LX.LT.NX ) THEN
               LX   = NX
               INFO = -16
               CALL XERBLA( 'IB03BD', -INFO )
               RETURN
            END IF
C           Workspace for IB01BD.
            IW1 = MAX( 2*LNOL*N + 2*N, LNOL*N + N2 + 7*N, L*NOBR*N +
     $                 MAX( LNOL*N + 2*N + ( M + ML )*NOBR + L,
     $                      2*LNOL*N + N2 + 8*N, N + 4*( MNO + N ) + 1,
     $                      MNO + 3*N + L ) )
            IF ( M.GT.0 ) THEN
               IW2 = L*NOBR*N + MNO*LDAC*( M*LDAC + 1 ) +
     $               MAX( LDAC**2, 4*M*LDAC + 1 )
            ELSE
               IW2 = 0
            END IF
            JWORK = ISV + ISAD + MAX( IW1, IW2 )
C           Workspace for IB01CD.
            IW1   = NSML*( N + 1 ) + 2*N + MAX( 2*N2, 4*N )
            IW2   = N*( N + 1 ) + 2*N + MAX( N*L*( N + 1 ) + 2*N2 + L*N,
     $                                       4*N )
            JWORK = MAX( JWORK, ISAD + 2 + N*( N + 1 + LDAC + M ) +
     $                   MAX( 5*N, 2, MIN( IW1, IW2 ) ) )
C           Workspace for TF01MX.
            JWORK = MAX( JWORK, NSML + ISAD + LDAC + 2*N + M )
C           Workspace for TB01VD.
            JWORK = MAX( JWORK, NSML + ISAD + N +
     $                   MAX( 1, N2*L + N*L + N,
     $                        N2 + MAX( N2 + N*MAX( N, L ) +
     $                                  6*N +  MIN( N, L ), N*M ) ) )
C           Workspace for MD03BD (whole optimization).
            IF ( M.GT.0 ) THEN
               IW1 = LDAC + M
            ELSE
               IW1 = L
            END IF
            IW1 = NSML + MAX( 2*NN, ISAD + 2*N + MAX( N*LDAC, IW1 ) )
            IF ( L.LE.1 .OR. BSN.EQ.0 ) THEN
               IW3 = 4*NX
               IW2 = IW3 + 1
            ELSE
               IW2 = BSN + MAX( 3*BSN + 1, LTHS )
               IF ( NSMP.GT.BSN ) THEN
                  IW2 = MAX( IW2, 4*LTHS + 1 )
                  IF ( NSMP.LT.2*BSN )
     $               IW2 = MAX( IW2, ( NSMP - BSN )*( L - 1 ) )
               END IF
               IW3 = LTHS*BSN + 2*NX + 2*MAX( BSN, LTHS )
            END IF
            JWORK = MAX( JWORK, NSML + NX +
     $                   MAX( 4, NSML +
     $                           MAX( NSML*( BSN + LTHS ) +
     $                                MAX( NSML + IW1, IW2 + NX ),
     $                                     NX*( BSN + LTHS ) + NX +
     $                                     MAX( NSML + IW1, NX + IW3 ) )
     $                      ) )
            IF ( LDWORK.LT.JWORK ) THEN
               INFO = -21
               DWORK(1) = JWORK
               CALL XERBLA( 'IB03BD', -INFO )
               RETURN
            END IF
         END IF
C
         BD = AC + LDAC*N
         IX = BD + LDAC*M
         IA = ISV
         IB = IA + LDAC*N
         IQ = IB + LDAC*M
         IF ( LSAME( JOBCK, 'N' ) ) THEN
            IRY   = IQ
            IS    = IQ
            IK    = IQ
            JWORK = IQ
         ELSE
            IRY   = IQ  + N2
            IS    = IRY + L*L
            IK    = IS  + N*L
            JWORK = IK  + N*L
         END IF
C
C        The workspace needed is defined for the options set above
C        in the PARAMETER statements.
C        Workspace:
C          need:  4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) +
C                 max( LDW1,LDW2 ), where,
C                 LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N,
C                              L*NOBR*N +
C                              max( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L,
C                                   2*(L*NOBR-L)*N+N*N+8*N,
C                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) )
C                 LDW2 >= 0,                                  if M = 0;
C                 LDW2 >= L*NOBR*N+M*NOBR*(N+L)*(M*(N+L)+1)+
C                         max( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0;
C          prefer: larger.
C        Integer workspace:  MAX(M*NOBR+N,M*(N+L)).
C
         CALL IB01BD( METHB, JOB, JOBCK, NOBR, N, M, L, NSMP, DWORK(IR),
     $                LDR, DWORK(IA), LDAC, DWORK(IA+N), LDAC,
     $                DWORK(IB), LDAC, DWORK(IB+N), LDAC, DWORK(IQ), N,
     $                DWORK(IRY), L, DWORK(IS), N, DWORK(IK), N, RCOND,
     $                IWORK, DWORK(JWORK), LDWORK-JWORK+1, BWORK,
     $                IWARNL, INFOL )
C
         IF( INFOL.EQ.-30 ) THEN
            INFO = -21
            DWORK(1) = DWORK(JWORK)
            CALL XERBLA( 'IB03BD', -INFO )
            RETURN
         END IF
         IF( INFOL.NE.0 ) THEN
            INFO = 100*INFOL
            RETURN
         END IF
         IF( IWARNL.NE.0 )
     $      IWARN = 100*IWARNL
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IRCNDB = 4
         IF ( LSAME( JOBCK, 'K' ) )
     $      IRCNDB = IRCNDB + 8
         CALL DCOPY( IRCNDB, DWORK(JWORK+1), 1, RCND(IRCND+1), 1 )
         IRCND = IRCND + IRCNDB
C
C        Copy the system matrices to the beginning of DWORK, to save
C        space, and redefine the pointers.
C
         CALL DCOPY( ISAD, DWORK(IA), 1, DWORK, 1 )
         IA  = 1
         IB  = IA  + LDAC*N
         IX0 = IB  + LDAC*M
         IV  = IX0 + N
C
C        Compute the initial condition of the system. On normal exit,
C           DWORK(i), i = JWORK+2:JWORK+1+N*N,
C           DWORK(j), j = JWORK+2+N*N:JWORK+1+N*N+L*N,  and
C           DWORK(k), k = JWORK+2+N*N+L*N:JWORK+1+N*N+L*N+N*M,
C        contain the transformed system matrices  At, Ct, and Bt,
C        respectively, corresponding to the real Schur form of the
C        estimated system state matrix  A. The transformation matrix is
C        stored in DWORK(IV:IV+N*N-1).
C        The workspace needed is defined for the options set above
C        in the PARAMETER statements.
C        Workspace:
C          need:   (N+L)*(N+M) + N + N*N + 2 + N*( N + M + L ) +
C                  max( 5*N, 2, min( LDW1, LDW2 ) ), where,
C                  LDW1 = NSMP*L*(N + 1) + 2*N + max( 2*N*N, 4*N),
C                  LDW2 = N*(N + 1) + 2*N +
C                         max( N*L*(N + 1) + 2*N*N + L*N, 4*N);
C          prefer: larger.
C        Integer workspace:  N.
C
         JWORK = IV + N2
         CALL IB01CD( 'X needed', COMUSE, JOBXD, N, M, L, NSMP,
     $                DWORK(IA), LDAC, DWORK(IB), LDAC, DWORK(IA+N),
     $                LDAC, DWORK(IB+N), LDAC, U, LDU, Y, LDY,
     $                DWORK(IX0), DWORK(IV), N, RCOND, IWORK,
     $                DWORK(JWORK), LDWORK-JWORK+1, IWARNL, INFOL )
C
         IF( INFOL.EQ.-26 ) THEN
            INFO = -21
            DWORK(1) = DWORK(JWORK)
            CALL XERBLA( 'IB03BD', -INFO )
            RETURN
         END IF
         IF( INFOL.EQ.1 )
     $      INFOL = 10
         IF( INFOL.NE.0 ) THEN
            INFO = 100*INFOL
            RETURN
         END IF
         IF( IWARNL.NE.0 )
     $      IWARN = 100*IWARNL
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IRCND  = IRCND + 1
         RCND(IRCND) = DWORK(JWORK+1)
C
C        Now, save the system matrices and x0 in the final location.
C
         IF ( IV.LT.AC ) THEN
            CALL DCOPY( ISAD+N, DWORK(IA), 1, DWORK(AC), 1 )
         ELSE
            DO 10 J = AC + ISAD + N - 1, AC, -1
               DWORK(J) = DWORK(IA+J-AC)
   10       CONTINUE
         END IF
C
C        Compute the output of the linear part.
C        Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L,
C                                                              if M > 0;
C                          NSMP*L + (N + L)*N + 2*N + L,       if M = 0;
C                   prefer larger.
C
         JWORK = IX + N
         CALL DCOPY(  N, DWORK(IX), 1, X(NTHS+1), 1 )
         CALL TF01MX( N, M, L, NSMP, DWORK(AC), LDAC, U, LDU, X(NTHS+1),
     $                DWORK(Z), NSMP, DWORK(JWORK), LDWORK-JWORK+1,
     $                INFO )
C
C        Convert the state-space representation to output normal form.
C        Workspace:
C          need:   NSMP*L + (N + L)*(N + M) + N +
C                  MAX(1, N*N*L + N*L + N, N*N +
C                      MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), N*M));
C          prefer: larger.
C
         CALL TB01VD( 'Apply', N, M, L, DWORK(AC), LDAC, DWORK(BD),
     $                LDAC, DWORK(AC+N), LDAC, DWORK(BD+N), LDAC,
     $                DWORK(IX), X(NTHS+1), LTHS, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFOL )
C
         IF( INFOL.GT.0 ) THEN
            INFO = INFOL + 4
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
      END IF
C
      LIPAR = 7
      IW1   = 0
      IW2   = 0
      IDIAG = AC
C
      IF ( INIT2 ) THEN
C
C        Initialize the nonlinear part.
C
         IF ( .NOT.INIT1 ) THEN
            BD = AC + LDAC*N
            IX = BD + LDAC*M
C
C           Convert the output normal form to state-space model.
C           Workspace: need NSMP*L + (N + L)*(2*N + M) + 2*N.
C           (NSMP*L locations are reserved for the output of the linear
C           part.)
C
            JWORK = IX + N
            CALL TB01VY( 'Apply', N, M, L, X(NTHS+1), LTHS, DWORK(AC),
     $                   LDAC, DWORK(BD), LDAC, DWORK(AC+N), LDAC,
     $                   DWORK(BD+N), LDAC, DWORK(IX), DWORK(JWORK),
     $                   LDWORK-JWORK+1, INFO )
C
C           Compute the output of the linear part.
C           Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L,
C                                                              if M > 0;
C                             NSMP*L + (N + L)*N + 2*N + L,    if M = 0;
C                      prefer larger.
C
            CALL TF01MX( N, M, L, NSMP, DWORK(AC), LDAC, U, LDU,
     $                   DWORK(IX), DWORK(Z), NSMP, DWORK(JWORK),
     $                   LDWORK-JWORK+1, INFO )
         END IF
C
C        Optimize the parameters of the nonlinear part.
C        Workspace:
C          need   NSMP*L + BSN +
C                 MAX( 4, NSMP +
C                      MAX( NSMP*BSN + MAX( 2*NN, 5*BSN + 1 ),
C                           BSN**2 + BSN + MAX( NSMP + 2*NN, 5*BSN ) ));
C          prefer larger.
C        Integer workspace:  NN*(L + 2) + 2.
C
         WORK(1) = ZERO
         CALL DCOPY( 3, WORK(1), 0, WORK(2), 1 )
C
C        Set the integer parameters needed, including the number of
C        neurons.
C
         IPAR(1) = NSMP
         IPAR(2) = L
         IPAR(3) = NN
         JWORK   = IDIAG + BSN
C
         DO 30 I = 0, L - 1
            CALL DCOPY( 4, SEED, 1, DWORK(JWORK), 1 )
            CALL MD03BD( 'Random initialization', SCALE, COND, NF01BE,
     $                   MD03BA, MD03BB, NSMP, BSN, ITMAX1, FACTOR,
     $                   NPRINT, IPAR, LIPAR, DWORK(Z), NSMP, Y(1,I+1),
     $                   LDY, X(I*BSN+1), DWORK(IDIAG), NFEV, NJEV,
     $                   TOL1, TOL1, GTOL, TOL, IWORK, DWORK(JWORK),
     $                   LDWORK-JWORK+1, IWARNL, INFOL )
            IF( INFOL.NE.0 ) THEN
               INFO = 10*INFOL
               RETURN
            END IF
            IF ( IWARNL.LT.0 ) THEN
               INFO  = INFOL
               IWARN = IWARNL
               GO TO 50
            ELSEIF ( IWARNL.GT.0 ) THEN
               IF ( IWARN.GT.100 ) THEN
                  IWARN = MAX( IWARN, ( IWARN/100 )*100 + 10*IWARNL )
               ELSE
                  IWARN = MAX( IWARN, 10*IWARNL )
               END IF
            END IF
            WORK(1) = MAX( WORK(1), DWORK(JWORK) )
            WORK(2) = MAX( WORK(2), DWORK(JWORK+1) )
            WORK(4) = MAX( WORK(4), DWORK(JWORK+3) )
            WORK(3) = WORK(3) + DWORK(JWORK+2)
            IW1     = NFEV + IW1
            IW2     = NJEV + IW2
   30    CONTINUE
C
      END IF
C
C     Main iteration.
C     Workspace:
C       need   NSMP*L + NX +
C              MAX( 4, NSMP*L +
C                      MAX( NSMP*L*( BSN + LTHS ) +
C                           MAX( NSMP*L + LDW1, LDW2 + NX ),
C                                NX*( BSN + LTHS ) + NX +
C                                MAX( NSMP*L + LDW1, NX + LDW3 ) ) ),
C              LDW0 = MAX( N*(N+L), N+M+L ),    if M > 0;
C              LDW0 = MAX( N*(N+L), L ),        if M = 0;
C              LDW1 = NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + LDW0);
C              LDW2 = 4*NX + 1,  if L <= 1 or BSN = 0; otherwise,
C              LDW2 = BSN + MAX(3*BSN+1,LTHS);
C              LDW2 = MAX(LDW2, 4*LTHS+1),        if NSMP > BSN;
C              LDW2 = MAX(LDW2, (NSMP-BSN)*(L-1)), if BSN < NSMP < 2*BSN;
C              LDW3 = 4*NX,                       if L <= 1 or BSN = 0;
C              LDW3 = LTHS*BSN + 2*NX + 2*MAX(BSN,LTHS),
C                                                 if L > 1 and BSN > 0;
C       prefer larger.
C     Integer workspace:  NX+L.
C
C     Set the integer parameters describing the Jacobian structure
C     and the number of neurons.
C
      IPAR(1) = LTHS
      IPAR(2) = L
      IPAR(3) = NSMP
      IPAR(4) = BSN
      IPAR(5) = M
      IPAR(6) = N
      IPAR(7) = NN
      JWORK   = IDIAG + NX
C
      CALL MD03BD( 'Given initialization', SCALE, COND, NF01BF,
     $             NF01BS, NF01BP, NSML, NX, ITMAX2, FACTOR, NPRINT,
     $             IPAR, LIPAR, U, LDU, Y, LDY, X, DWORK(IDIAG), NFEV,
     $             NJEV, TOL2, TOL2, GTOL, TOL, IWORK, DWORK(JWORK),
     $             LDWORK-JWORK+1, IWARNL, INFO )
      IF( INFO.NE.0 )
     $   RETURN
C
      DO 40 I = 1, NX + L
         IWORK(I+3) = IWORK(I)
   40 CONTINUE
C
   50 CONTINUE
      IWORK(1) = IW1 + NFEV
      IWORK(2) = IW2 + NJEV
      IF ( IWARNL.LT.0 ) THEN
         IWARN = IWARNL
      ELSE
         IWARN = IWARN + IWARNL
      END IF
      CALL DCOPY( 4, DWORK(JWORK), 1, DWORK, 1 )
      IF ( INIT2 )
     $   CALL DCOPY( 4, WORK, 1, DWORK(5), 1 )
      IF ( INIT1 ) THEN
         IWORK(3) = IRCND
         CALL DCOPY( IRCND, RCND, 1, DWORK(9), 1 )
      ELSE
         IWORK(3) = 0
      END IF
C
      RETURN
C
C *** Last line of IB03BD ***
      END
