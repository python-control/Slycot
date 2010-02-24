      SUBROUTINE IB03AD( INIT, ALG, STOR, NOBR, M, L, NSMP, N, NN,
     $                   ITMAX1, ITMAX2, NPRINT, U, LDU, Y, LDY, X, LX,
     $                   TOL1, TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     Levenberg-Marquardt algorithm. Conjugate gradients (CG) or
C     Cholesky algorithms are used to solve linear systems of equations.
C     The Wiener system is represented as
C
C        x(t+1) = A*x(t) + B*u(t)
C        z(t)   = C*x(t) + D*u(t),
C
C        y(t)   = f(z(t),wb(1:L)),
C
C     where t = 1, 2, ..., NSMP, and f is a nonlinear function,
C     evaluated by the SLICOT Library routine NF01AY. The parameter
C     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ),
C     where wb(i), i = 1 : L, correspond to the nonlinear part, and
C     theta corresponds to the linear part. See SLICOT Library routine
C     NF01AD for further details.
C
C     The sum of squares of the error functions, defined by
C
C        e(t) = y(t) - Y(t),  t = 1, 2, ..., NSMP,
C
C     is minimized, where Y(t) is the measured output vector. The
C     functions and their Jacobian matrices are evaluated by SLICOT
C     Library routine NF01BB (the FCN routine in the call of MD03AD).
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
C             by SLICOT Library routine NF01BA (used as a second FCN
C             routine in the MD03AD call for the initialization step,
C             see METHOD).
C
C     ALG     CHARACTER*1
C             Specifies the algorithm used for solving the linear
C             systems involving a Jacobian matrix J, as follows:
C             = 'D' :  a direct algorithm, which computes the Cholesky
C                      factor of the matrix J'*J + par*I is used, where
C                      par is the Levenberg factor;
C             = 'I' :  an iterative Conjugate Gradients algorithm, which
C                      only needs the matrix J, is used.
C             In both cases, matrix J is stored in a compressed form.
C
C     STOR    CHARACTER*1
C             If ALG = 'D', specifies the storage scheme for the
C             symmetric matrix J'*J, as follows:
C             = 'F' :  full storage is used;
C             = 'P' :  packed storage is used.
C             The option STOR = 'F' usually ensures a faster execution.
C             This parameter is not relevant if ALG = 'I'.
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
C             FCN routine (NF01BA and/or NF01BB). If NPRINT <= 0, no
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
C             squares, for the initialization step of nonlinear part.
C             Termination occurs when the actual relative reduction in
C             the sum of squares is at most TOL1. In addition, if
C             ALG = 'I', TOL1 also measures the relative residual of
C             the solutions computed by the CG algorithm (for the
C             initialization step). Termination of a CG process occurs
C             when the relative residual is at most TOL1.
C             If the user sets  TOL1 < 0,  then  SQRT(EPS)  is used
C             instead TOL1, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH).
C             This parameter is ignored if INIT is 'N' or 'L'.
C
C     TOL2    DOUBLE PRECISION
C             If TOL2 >= 0, TOL2 is the tolerance which measures the
C             relative error desired in the sum of squares, for the
C             whole optimization process. Termination occurs when the
C             actual relative reduction in the sum of squares is at
C             most TOL2.
C             If ALG = 'I', TOL2 also measures the relative residual of
C             the solutions computed by the CG algorithm (for the whole
C             optimization). Termination of a CG process occurs when the
C             relative residual is at most TOL2.
C             If the user sets  TOL2 < 0,  then  SQRT(EPS)  is used
C             instead TOL2. This default value could require many
C             iterations, especially if TOL1 is larger. If INIT = 'S'
C             or 'B', it is advisable that TOL2 be larger than TOL1,
C             and spend more time with cheaper iterations.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX( 3, LIW1, LIW2 )), where
C             LIW1 = LIW2 = 0,  if INIT = 'S' or 'N'; otherwise,
C             LIW1 = M+L;
C             LIW2 = MAX(M*NOBR+N,M*(N+L)).
C             On output, if INFO = 0, IWORK(1) and IWORK(2) return the
C             (total) number of function and Jacobian evaluations,
C             respectively (including the initialization step, if it was
C             performed), and if INIT = 'L' or INIT = 'B', IWORK(3)
C             specifies how many locations of DWORK contain reciprocal
C             condition number estimates (see below); otherwise,
C             IWORK(3) = 0.
C
C     DWORK   DOUBLE PRECISION array dimesion (LDWORK)
C             On entry, if desired, and if INIT = 'S' or 'B', the
C             entries DWORK(1:4) are set to initialize the random
C             numbers generator for the nonlinear part parameters (see
C             the description of the argument XINIT of SLICOT Library
C             routine MD03AD); this enables to obtain reproducible
C             results. The same seed is used for all outputs.
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, DWORK(2) returns the residual error norm (the
C             sum of squares), DWORK(3) returns the number of iterations
C             performed, DWORK(4) returns the number of conjugate
C             gradients iterations performed, and DWORK(5) returns the
C             final Levenberg factor, for optimizing the parameters of
C             both the linear part and the static nonlinearity part.
C             If INIT = 'S' or INIT = 'B' and INFO = 0, then the
C             elements DWORK(6) to DWORK(10) contain the corresponding
C             five values for the initialization step (see METHOD).
C             (If L > 1, DWORK(10) contains the maximum of the Levenberg
C             factors for all outputs.) If INIT = 'L' or INIT = 'B', and
C             INFO = 0, DWORK(11) to DWORK(10+IWORK(3)) contain
C             reciprocal condition number estimates set by SLICOT
C             Library routines IB01AD, IB01BD, and IB01CD.
C             On exit, if  INFO = -23,  DWORK(1)  returns the minimum
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
C             LW2 = NSMP*L +
C                   MAX( 5, NSMP + 2*BSN + NSMP*BSN +
C                           MAX( 2*NN + BSN, LDW7 ) );
C                 LDW7 = BSN*BSN,       if ALG = 'D' and STOR = 'F';
C                 LDW7 = BSN*(BSN+1)/2, if ALG = 'D' and STOR = 'P';
C                 LDW7 = 3*BSN + NSMP,  if ALG = 'I';
C             LW3 = MAX( LDW8, NSMP*L + (N+L)*(2*N+M) + 2*N );
C                 LDW8 = NSMP*L + (N+L)*(N+M) + 3*N+M+L,  if M > 0;
C                 LDW8 = NSMP*L + (N+L)*N + 2*N+L,        if M = 0;
C             LW4 = MAX( 5, NSMP*L + 2*NX + NSMP*L*( BSN + LTHS ) +
C                           MAX( L1 + NX, NSMP*L + L1, L2 ) ),
C                  L0 = MAX( N*(N+L), N+M+L ),    if M > 0;
C                  L0 = MAX( N*(N+L), L ),        if M = 0;
C                  L1 = NSMP*L + MAX( 2*NN, (N+L)*(N+M) + 2*N + L0);
C                  L2 = NX*NX,          if ALG = 'D' and STOR = 'F';
C                  L2 = NX*(NX+1)/2,    if ALG = 'D' and STOR = 'P';
C                  L2 = 3*NX + NSMP*L,  if ALG = 'I',
C                  with BSN  = NN*( L + 2 ) + 1,
C                       LTHS = N*( L + M + 1 ) + L*M.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             < 0:  the user set IFLAG = IWARN in (one of) the
C                   subroutine(s) FCN, i.e., NF01BA, if INIT = 'S'
C                   or 'B', and/or NF01BB; this value cannot be returned
C                   without changing the FCN routine(s);
C                   otherwise, IWARN has the value k*100 + j*10 + i,
C                   where k is defined below, i refers to the whole
C                   optimization process, and j refers to the
C                   initialization step (j = 0, if INIT = 'L' or 'N'),
C                   and the possible values for i and j have the
C                   following meaning (where TOL* denotes TOL1 or TOL2,
C                   and similarly for ITMAX*):
C             = 1:  the number of iterations has reached ITMAX* without
C                   satisfying the convergence condition;
C             = 2:  if alg = 'I' and in an iteration of the Levenberg-
C                   Marquardt algorithm, the CG algorithm finished
C                   after 3*NX iterations (or 3*(lin1-1) iterations, for
C                   the initialization phase), without achieving the
C                   precision required in the call;
C             = 3:  the cosine of the angle between the vector of error
C                   function values and any column of the Jacobian is at
C                   most FACTOR*EPS in absolute value (FACTOR = 100);
C             = 4:  TOL* is too small: no further reduction in the sum
C                   of squares is possible.
C             The digit k is normally 0, but if INIT = 'L' or 'B', it
C             can have a value in the range 1 to 6 (see IB01AD, IB01BD
C             and IB01CD). In all these cases, the entries DWORK(1:5),
C             DWORK(6:10) (if INIT = 'S' or 'B'), and
C             DWORK(11:10+IWORK(3)) (if INIT = 'L' or 'B'), are set as
C             described above.
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
C             = 3:  ALG = 'D' and SLICOT Library routines MB02XD or
C                   NF01BU (or NF01BV, if INIT = 'S' or 'B') or
C                   ALG = 'I' and SLICOT Library routines MB02WD or
C                   NF01BW (or NF01BX, if INIT = 'S' or 'B') returned
C                   with INFO <> 0.
C             In addition, if INIT = 'L' or 'B', i could also be
C             = 4:  if a Lyapunov equation could not be solved;
C             = 5:  if the identified linear system is unstable;
C             = 6:  if the QR algorithm failed on the state matrix
C                   of the identified linear system.
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
C     nonlinearity using the Levenberg-Marquardt algorithm [1].
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
C     [1] Kelley, C.T.
C         Iterative Methods for Optimization.
C         Society for Industrial and Applied Mathematics (SIAM),
C         Philadelphia (Pa.), 1999.
C
C     [2] Peeters, R.L.M., Hanzon, B., and Olivi, M.
C         Balanced realizations of discrete-time stable all-pass
C         systems and the tangential Schur algorithm.
C         Proceedings of the European Control Conference,
C         31 August - 3 September 1999, Karlsruhe, Germany.
C         Session CP-6, Discrete-time Systems, 1999.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Mar. 2002, Apr. 2002, Feb. 2004, March 2005, Nov. 2005.
C
C     KEYWORDS
C
C     Conjugate gradients, least-squares approximation,
C     Levenberg-Marquardt algorithm, matrix operations, optimization.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     The upper triangular part is used in MD03AD;
      CHARACTER         UPLO
      PARAMETER         ( UPLO = 'U' )
C     For INIT = 'L' or 'B', additional parameters are set:
C     The following six parameters are used in the call of IB01AD;
      CHARACTER         IALG, BATCH, CONCT, CTRL, JOBD, METH
      PARAMETER         ( IALG  = 'Fast QR',     BATCH = 'One batch',
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
      CHARACTER         ALG, INIT, STOR
      INTEGER           INFO, ITMAX1, ITMAX2, IWARN, L, LDU, LDWORK,
     $                  LDY, LX, M, N, NN, NOBR, NPRINT, NSMP
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), U(LDU, *), X(*), Y(LDY, *)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      INTEGER           AC, BD, BSN, I, IA, IB, IK, INFOL, IQ, IR,
     $                  IRCND, IRCNDB, IRY, IS, ISAD, ISV, IV, IW1, IW2,
     $                  IWARNL, IX, IX0, J, JWORK, LDAC, LDR, LIPAR,
     $                  LNOL, LTHS, ML, MNO, N2, NFEV, NJEV, NS, NSML,
     $                  NTHS, NX, WRKOPT, Z
      LOGICAL           CHOL, FULL, INIT1, INIT2
C     .. Local Arrays ..
      LOGICAL           BWORK(1)
      INTEGER           IPAR(7)
      DOUBLE PRECISION  RCND(16), SEED(4), WORK(5)
C     .. External Functions ..
      EXTERNAL          LSAME
      LOGICAL           LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY,  IB01AD, IB01BD, IB01CD, MD03AD, NF01BA,
     $                  NF01BB, NF01BU, NF01BV, NF01BW, NF01BX, TB01VD,
     $                  TB01VY, TF01MX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      CHOL  = LSAME( ALG,  'D' )
      FULL  = LSAME( STOR, 'F' )
      INIT1 = LSAME( INIT, 'B' ) .OR. LSAME( INIT, 'L' )
      INIT2 = LSAME( INIT, 'B' ) .OR. LSAME( INIT, 'S' )
C
      ML    = M + L
      INFO  = 0
      IWARN = 0
      IF ( .NOT.( INIT1 .OR. INIT2 .OR. LSAME( INIT, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( CHOL .OR. LSAME( ALG,  'I' ) ) ) THEN
         INFO = -2
      ELSEIF ( CHOL .AND. .NOT.( FULL .OR. LSAME( STOR, 'P' ) ) ) THEN
         INFO = -3
      ELSEIF ( INIT1 .AND. NOBR.LE.0 ) THEN
         INFO = -4
      ELSEIF ( M.LT.0 ) THEN
         INFO = -5
      ELSEIF ( L.LT.0 .OR. ( INIT1 .AND. L.EQ.0 ) ) THEN
         INFO = -6
      ELSEIF ( NSMP.LT.0 .OR.
     $         ( INIT1 .AND. NSMP.LT.2*( ML + 1 )*NOBR - 1 ) ) THEN
         INFO = -7
      ELSEIF ( ( N.LT.0 .AND. .NOT.INIT1 ) .OR.
     $       ( ( N.EQ.0 .OR. N.GE.NOBR ) .AND. INIT1 ) ) THEN
         INFO = -8
      ELSEIF ( NN.LT.0 ) THEN
         INFO = -9
      ELSEIF ( INIT2 .AND. ( ITMAX1.LT.0 ) ) THEN
         INFO = -10
      ELSEIF ( ITMAX2.LT.0 ) THEN
         INFO = -11
      ELSEIF ( LDU.LT.MAX( 1, NSMP ) ) THEN
         INFO = -14
      ELSEIF ( LDY.LT.MAX( 1, NSMP ) ) THEN
         INFO = -16
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
C           Workspace for MD03AD (initialization of the nonlinear part).
            IF ( CHOL ) THEN
               IF ( FULL ) THEN
                  IW1 = BSN**2
               ELSE
                  IW1 = ( BSN*( BSN + 1 ) )/2
               END IF
            ELSE
               IW1 = 3*BSN + NSMP
            END IF
            JWORK = MAX( JWORK, NSML +
     $                   MAX( 5, NSMP + 2*BSN + NSMP*BSN +
     $                        MAX( 2*NN + BSN, IW1 ) ) )
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
               INFO = -18
               CALL XERBLA( 'IB03AD', -INFO )
               RETURN
            END IF
C
C           Workspace for MD03AD (whole optimization).
C
            IF ( M.GT.0 ) THEN
               IW1 = LDAC + M
            ELSE
               IW1 = L
            END IF
            IW1 = NSML + MAX( 2*NN, ISAD + 2*N + MAX( N*LDAC, IW1 ) )
            IF ( CHOL ) THEN
               IF ( FULL ) THEN
                  IW2 = NX**2
               ELSE
                  IW2 = ( NX*( NX + 1 ) )/2
               END IF
            ELSE
               IW2 = 3*NX + NSML
            END IF
            JWORK = MAX( JWORK,
     $                   5, NSML + 2*NX + NSML*( BSN + LTHS ) +
     $                      MAX( IW1 + NX, NSML + IW1, IW2 ) )
         END IF
C
         IF ( LDWORK.LT.JWORK ) THEN
            INFO = -23
            DWORK(1) = JWORK
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'IB03AD', -INFO )
         RETURN
      ENDIF
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
         CALL IB01AD( METH, IALG, JOBD, BATCH, CONCT, CTRL, NOBR, M, L,
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
               INFO = -18
               CALL XERBLA( 'IB03AD', -INFO )
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
C           Workspace for MD03AD (whole optimization).
            IF ( M.GT.0 ) THEN
               IW1 = LDAC + M
            ELSE
               IW1 = L
            END IF
            IW1 = NSML + MAX( 2*NN, ISAD + 2*N + MAX( N*LDAC, IW1 ) )
            IF ( CHOL ) THEN
               IF ( FULL ) THEN
                  IW2 = NX**2
               ELSE
                  IW2 = ( NX*( NX + 1 ) )/2
               END IF
            ELSE
               IW2 = 3*NX + NSML
            END IF
            JWORK = MAX( JWORK,
     $                   5, NSML + 2*NX + NSML*( BSN + LTHS ) +
     $                      MAX( IW1 + NX, NSML + IW1, IW2 ) )
            IF ( LDWORK.LT.JWORK ) THEN
               INFO = -23
               DWORK(1) = JWORK
               CALL XERBLA( 'IB03AD', -INFO )
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
            INFO = -23
            DWORK(1) = DWORK(JWORK)
            CALL XERBLA( 'IB03AD', -INFO )
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
            INFO = -23
            DWORK(1) = DWORK(JWORK)
            CALL XERBLA( 'IB03AD', -INFO )
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
            DO 5 J = AC + ISAD + N - 1, AC, -1
               DWORK(J) = DWORK(IA+J-AC)
    5       CONTINUE
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
            INFO = INFOL + 3
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
      END IF
C
      LIPAR = 7
      IW1   = 0
      IW2   = 0
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
C          need   NSMP*L +
C                 MAX( 5, NSMP + 2*BSN + NSMP*BSN +
C                         MAX( 2*NN + BSN, DW( sol ) ) ),
C                 where, if ALG = 'D',
C                      DW( sol ) = BSN*BSN,        if STOR = 'F';
C                      DW( sol ) = BSN*(BSN+1)/2,  if STOR = 'P';
C                 and  DW( sol ) = 3*BSN + NSMP,   if ALG  = 'I';
C          prefer larger.
C
         JWORK   = AC
         WORK(1) = ZERO
         CALL DCOPY( 4, WORK(1), 0, WORK(2), 1 )
C
C        Set the integer parameters needed, including the number of
C        neurons.
C
         IPAR(1) = NSMP
         IPAR(2) = L
         IPAR(3) = NN
C
         DO 10 I = 0, L - 1
            CALL DCOPY( 4, SEED, 1, DWORK(JWORK), 1 )
            IF ( CHOL ) THEN
               CALL MD03AD( 'Random initialization', ALG, STOR, UPLO,
     $                      NF01BA, NF01BV, NSMP, BSN, ITMAX1, NPRINT,
     $                      IPAR, LIPAR, DWORK(Z), NSMP, Y(1,I+1), LDY,
     $                      X(I*BSN+1), NFEV, NJEV, TOL1, TOL1,
     $                      DWORK(JWORK), LDWORK-JWORK+1, IWARNL,
     $                      INFOL )
            ELSE
               CALL MD03AD( 'Random initialization', ALG, STOR, UPLO,
     $                      NF01BA, NF01BX, NSMP, BSN, ITMAX1, NPRINT,
     $                      IPAR, LIPAR, DWORK(Z), NSMP, Y(1,I+1), LDY,
     $                      X(I*BSN+1), NFEV, NJEV, TOL1, TOL1,
     $                      DWORK(JWORK), LDWORK-JWORK+1, IWARNL,
     $                      INFOL )
            END IF
C
            IF( INFOL.NE.0 ) THEN
               INFO = 10*INFOL
               RETURN
            END IF
            IF ( IWARNL.LT.0 ) THEN
               INFO  = INFOL
               IWARN = IWARNL
               GO TO 20
            ELSEIF ( IWARNL.GT.0 ) THEN
               IF ( IWARN.GT.100 ) THEN
                  IWARN = MAX( IWARN, ( IWARN/100 )*100 + 10*IWARNL )
               ELSE
                  IWARN = MAX( IWARN, 10*IWARNL )
               END IF
            END IF
            WORK(1) = MAX( WORK(1), DWORK(JWORK) )
            WORK(2) = MAX( WORK(2), DWORK(JWORK+1) )
            WORK(5) = MAX( WORK(5), DWORK(JWORK+4) )
            WORK(3) = WORK(3) + DWORK(JWORK+2)
            WORK(4) = WORK(4) + DWORK(JWORK+3)
            IW1     = NFEV + IW1
            IW2     = NJEV + IW2
   10    CONTINUE
C
      ENDIF
C
C     Main iteration.
C     Workspace: need   MAX( 5, NFUN + 2*NX + NFUN*( BSN + LTHS ) +
C                            MAX( LDW1 + NX, NFUN + LDW1, DW( sol ) ) ),
C                       where NFUN = NSMP*L, and
C                       LDW1 = NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N +
C                                          MAX( N*(N + L), N + M + L )),
C                                                              if M > 0,
C                       LDW1 = NFUN + MAX( 2*NN, (N + L)*N + 2*N +
C                                          MAX( N*(N + L), L ) ),
C                                                              if M = 0;
C                       if ALG = 'D',
C                             DW( sol ) = NX*NX,        if STOR = 'F';
C                             DW( sol ) = NX*(NX+1)/2,  if STOR = 'P';
C                       and   DW( sol ) = 3*NX + NFUN,  if ALG  = 'I',
C                       and DW( f ) is the workspace needed by the
C                       subroutine f;
C                prefer larger.
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
C
      IF ( CHOL ) THEN
         CALL MD03AD( 'Given initialization', ALG, STOR, UPLO, NF01BB,
     $                NF01BU, NSML, NX, ITMAX2, NPRINT, IPAR, LIPAR,
     $                U, LDU, Y, LDY, X, NFEV, NJEV, TOL2, TOL2,
     $                DWORK, LDWORK, IWARNL, INFO )
      ELSE
         CALL MD03AD( 'Given initialization', ALG, STOR, UPLO, NF01BB,
     $                NF01BW, NSML, NX, ITMAX2, NPRINT, IPAR, LIPAR,
     $                U, LDU, Y, LDY, X, NFEV, NJEV, TOL2, TOL2,
     $                DWORK, LDWORK, IWARNL, INFO )
      END IF
C
      IF( INFO.NE.0 )
     $   RETURN
C
   20 CONTINUE
      IWORK(1) = IW1 + NFEV
      IWORK(2) = IW2 + NJEV
      IF ( IWARNL.LT.0 ) THEN
         IWARN = IWARNL
      ELSE
         IWARN = IWARN + IWARNL
      END IF
      IF ( INIT2 )
     $   CALL DCOPY( 5, WORK, 1, DWORK(6), 1 )
      IF ( INIT1 ) THEN
         IWORK(3) = IRCND
         CALL DCOPY( IRCND, RCND, 1, DWORK(11), 1 )
      ELSE
         IWORK(3) = 0
      END IF
      RETURN
C
C *** Last line of IB03AD ***
      END
