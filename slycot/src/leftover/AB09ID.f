      SUBROUTINE AB09ID( DICO, JOBC, JOBO, JOB, WEIGHT, EQUIL, ORDSEL,
     $                   N, M, P, NV, PV, NW, MW, NR, ALPHA, ALPHAC,
     $                   ALPHAO, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                   AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $                   NS, HSV, TOL1, TOL2, IWORK, DWORK, LDWORK,
     $                   IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for an original
C     state-space representation (A,B,C,D) by using the frequency
C     weighted square-root or balancing-free square-root
C     Balance & Truncate (B&T) or Singular Perturbation Approximation
C     (SPA) model reduction methods. The algorithm tries to minimize
C     the norm of the frequency-weighted error
C
C           ||V*(G-Gr)*W||
C
C     where G and Gr are the transfer-function matrices of the original
C     and reduced order models, respectively, and V and W are
C     frequency-weighting transfer-function matrices. V and W must not
C     have poles on the imaginary axis for a continuous-time
C     system or on the unit circle for a discrete-time system.
C     If G is unstable, only the ALPHA-stable part of G is reduced.
C     In case of possible pole-zero cancellations in V*G and/or G*W,
C     the absolute values of parameters ALPHAO and/or ALPHAC must be
C     different from 1.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOBC    CHARACTER*1
C             Specifies the choice of frequency-weighted controllability
C             Grammian as follows:
C             = 'S': choice corresponding to a combination method [4]
C                    of the approaches of Enns [1] and Lin-Chiu [2,3];
C             = 'E': choice corresponding to the stability enhanced
C                    modified combination method of [4].
C
C     JOBO    CHARACTER*1
C             Specifies the choice of frequency-weighted observability
C             Grammian as follows:
C             = 'S': choice corresponding to a combination method [4]
C                    of the approaches of Enns [1] and Lin-Chiu [2,3];
C             = 'E': choice corresponding to the stability enhanced
C                    modified combination method of [4].
C
C     JOB     CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root Balance & Truncate method;
C             = 'F':  use the balancing-free square-root
C                     Balance & Truncate method;
C             = 'S':  use the square-root Singular Perturbation
C                     Approximation method;
C             = 'P':  use the balancing-free square-root
C                     Singular Perturbation Approximation method.
C
C     WEIGHT  CHARACTER*1
C             Specifies the type of frequency weighting, as follows:
C             = 'N':  no weightings are used (V = I, W = I);
C             = 'L':  only left weighting V is used (W = I);
C             = 'R':  only right weighting W is used (V = I);
C             = 'B':  both left and right weightings V and W are used.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e., the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NV      (input) INTEGER
C             The order of the matrix AV. Also the number of rows of
C             the matrix BV and the number of columns of the matrix CV.
C             NV represents the dimension of the state vector of the
C             system with the transfer-function matrix V.  NV >= 0.
C
C     PV      (input) INTEGER
C             The number of rows of the matrices CV and DV.  PV >= 0.
C             PV represents the dimension of the output vector of the
C             system with the transfer-function matrix V.
C
C     NW      (input) INTEGER
C             The order of the matrix AW. Also the number of rows of
C             the matrix BW and the number of columns of the matrix CW.
C             NW represents the dimension of the state vector of the
C             system with the transfer-function matrix W.  NW >= 0.
C
C     MW      (input) INTEGER
C             The number of columns of the matrices BW and DW.  MW >= 0.
C             MW represents the dimension of the input vector of the
C             system with the transfer-function matrix W.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of the
C             resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. For a system with NU ALPHA-unstable
C             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N),
C             NR is set as follows: if ORDSEL = 'F', NR is equal to
C             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order
C             on entry, NMIN is the number of frequency-weighted Hankel
C             singular values greater than NS*EPS*S1, EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH)
C             and S1 is the largest Hankel singular value (computed
C             in HSV(1)); NR can be further reduced to ensure
C             HSV(NR-NU) > HSV(NR+1-NU);
C             if ORDSEL = 'A', NR is the sum of NU and the number of
C             Hankel singular values greater than MAX(TOL1,NS*EPS*S1).
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the ALPHA-stability boundary for the eigenvalues
C             of the state dynamics matrix A. For a continuous-time
C             system (DICO = 'C'), ALPHA <= 0 is the boundary value for
C             the real parts of eigenvalues, while for a discrete-time
C             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the
C             boundary value for the moduli of eigenvalues.
C             The ALPHA-stability domain does not include the boundary.
C
C     ALPHAC  (input) DOUBLE PRECISION
C             Combination method parameter for defining the
C             frequency-weighted controllability Grammian (see METHOD);
C             ABS(ALPHAC) <= 1.
C
C     ALPHAO  (input) DOUBLE PRECISION
C             Combination method parameter for defining the
C             frequency-weighted observability Grammian (see METHOD);
C             ABS(ALPHAO) <= 1.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the
C             reduced order system.
C             The resulting A has a block-diagonal form with two blocks.
C             For a system with NU ALPHA-unstable eigenvalues and
C             NS ALPHA-stable eigenvalues (NU+NS = N), the leading
C             NU-by-NU block contains the unreduced part of A
C             corresponding to ALPHA-unstable eigenvalues.
C             The trailing (NR+NS-N)-by-(NR+NS-N) block contains
C             the reduced part of A corresponding to ALPHA-stable
C             eigenvalues.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, if INFO = 0, the leading NR-by-M part of this
C             array contains the input/state matrix Br of the reduced
C             order system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, if INFO = 0, the leading P-by-NR part of this
C             array contains the state/output matrix Cr of the reduced
C             order system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the original input/output matrix D.
C             On exit, if INFO = 0, the leading P-by-M part of this
C             array contains the input/output matrix Dr of the reduced
C             order system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV)
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV
C             part of this array must contain the state matrix AV of
C             the system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading NVR-by-NVR part of this array
C             contains the state matrix of a minimal realization of V
C             in a real Schur form. NVR is returned in IWORK(2).
C             AV is not referenced if WEIGHT = 'R' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDAV    INTEGER
C             The leading dimension of array AV.
C             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDAV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P)
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part
C             of this array must contain the input matrix BV of the
C             system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading NVR-by-P part of this array contains
C             the input matrix of a minimal realization of V.
C             BV is not referenced if WEIGHT = 'R' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDBV    INTEGER
C             The leading dimension of array BV.
C             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDBV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV)
C             On entry, if WEIGHT = 'L' or 'B', the leading PV-by-NV
C             part of this array must contain the output matrix CV of
C             the system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading PV-by-NVR part of this array
C             contains the output matrix of a minimal realization of V.
C             CV is not referenced if WEIGHT = 'R' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDCV    INTEGER
C             The leading dimension of array CV.
C             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B';
C             LDCV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P)
C             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this
C             array must contain the feedthrough matrix DV of the system
C             with the transfer-function matrix V.
C             DV is not referenced if WEIGHT = 'R' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDDV    INTEGER
C             The leading dimension of array DV.
C             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B';
C             LDDV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW
C             part of this array must contain the state matrix AW of
C             the system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading NWR-by-NWR part of this array
C             contains the state matrix of a minimal realization of W
C             in a real Schur form. NWR is returned in IWORK(3).
C             AW is not referenced if WEIGHT = 'L' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDAW    INTEGER
C             The leading dimension of array AW.
C             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDAW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,MW)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-MW
C             part of this array must contain the input matrix BW of the
C             system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading NWR-by-MW part of this array
C             contains the input matrix of a minimal realization of W.
C             BW is not referenced if WEIGHT = 'L' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDBW    INTEGER
C             The leading dimension of array BW.
C             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDBW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part
C             of this array must contain the output matrix CW of the
C             system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and
C             INFO = 0, the leading M-by-NWR part of this array contains
C             the output matrix of a minimal realization of W.
C             CW is not referenced if WEIGHT = 'L' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDCW    INTEGER
C             The leading dimension of array CW.
C             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDCW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW)
C             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this
C             array must contain the feedthrough matrix DW of the system
C             with the transfer-function matrix W.
C             DW is not referenced if WEIGHT = 'L' or 'N',
C             or MIN(N,M,P) = 0.
C
C     LDDW    INTEGER
C             The leading dimension of array DW.
C             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDDW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of this array contain
C             the frequency-weighted Hankel singular values, ordered
C             decreasingly, of the ALPHA-stable part of the original
C             system.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*S1, where c is a constant in the
C             interval [0.00001,0.001], and S1 is the largest
C             frequency-weighted Hankel singular value of the
C             ALPHA-stable part of the original system (computed
C             in HSV(1)).
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = NS*EPS*S1, where NS is the number of
C             ALPHA-stable eigenvalues of A and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the ALPHA-stable part of the given system.
C             The recommended value is TOL2 = NS*EPS*S1.
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension
C             ( MAX( 3, LIWRK1, LIWRK2, LIWRK3 ) ), where
C             LIWRK1 = 0,             if JOB = 'B';
C             LIWRK1 = N,             if JOB = 'F';
C             LIWRK1 = 2*N,           if JOB = 'S' or 'P';
C             LIWRK2 = 0,             if WEIGHT = 'R' or 'N' or  NV = 0;
C             LIWRK2 = NV+MAX(P,PV),  if WEIGHT = 'L' or 'B' and NV > 0;
C             LIWRK3 = 0,             if WEIGHT = 'L' or 'N' or  NW = 0;
C             LIWRK3 = NW+MAX(M,MW),  if WEIGHT = 'R' or 'B' and NW > 0.
C             On exit, if INFO = 0, IWORK(1) contains the order of a
C             minimal realization of the stable part of the system,
C             IWORK(2) and IWORK(3) contain the actual orders
C             of the state space realizations of V and W, respectively.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( LMINL, LMINR, LRCF,
C                            2*N*N + MAX( 1, LLEFT, LRIGHT, 2*N*N+5*N,
C                                         N*MAX(M,P) ) ),
C             where
C             LMINL  = 0, if WEIGHT = 'R' or 'N' or NV = 0; otherwise,
C             LMINL  = MAX(LLCF,NV+MAX(NV,3*P))           if P =  PV;
C             LMINL  = MAX(P,PV)*(2*NV+MAX(P,PV))+
C                      MAX(LLCF,NV+MAX(NV,3*P,3*PV))      if P <> PV;
C             LRCF   = 0, and
C             LMINR  = 0, if WEIGHT = 'L' or 'N' or NW = 0; otherwise,
C             LMINR  = NW+MAX(NW,3*M)                     if M =  MW;
C             LMINR  = 2*NW*MAX(M,MW)+NW+MAX(NW,3*M,3*MW) if M <> MW;
C             LLCF   = PV*(NV+PV)+PV*NV+MAX(NV*(NV+5), PV*(PV+2),
C                                           4*PV, 4*P);
C             LRCF   = MW*(NW+MW)+MAX(NW*(NW+5),MW*(MW+2),4*MW,4*M)
C             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5)
C                              if WEIGHT = 'L' or 'B' and PV > 0;
C             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0;
C             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5)
C                              if WEIGHT = 'R' or 'B' and MW > 0;
C             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than NSMIN, the sum of the order of the
C                   ALPHA-unstable part and the order of a minimal
C                   realization of the ALPHA-stable part of the given
C                   system; in this case, the resulting NR is set equal
C                   to NSMIN;
C             = 2:  with ORDSEL = 'F', the selected order NR corresponds
C                   to repeated singular values for the ALPHA-stable
C                   part, which are neither all included nor all
C                   excluded from the reduced model; in this case, the
C                   resulting NR is automatically decreased to exclude
C                   all repeated singular values;
C             = 3:  with ORDSEL = 'F', the selected order NR is less
C                   than the order of the ALPHA-unstable part of the
C                   given system; in this case NR is set equal to the
C                   order of the ALPHA-unstable part.
C             = 10+K:  K violations of the numerical stability condition
C                   occured during the assignment of eigenvalues in the
C                   SLICOT Library routines SB08CD and/or SB08DD.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the computation of the ordered real Schur form of A
C                   failed;
C             = 2:  the separation of the ALPHA-stable/unstable
C                   diagonal blocks failed because of very close
C                   eigenvalues;
C             = 3:  the reduction to a real Schur form of the state
C                   matrix of a minimal realization of V failed;
C             = 4:  a failure was detected during the ordering of the
C                   real Schur form of the state matrix of a minimal
C                   realization of V or in the iterative process to
C                   compute a left coprime factorization with inner
C                   denominator;
C             = 5:  if DICO = 'C' and the matrix AV has an observable
C                   eigenvalue on the imaginary axis, or DICO = 'D' and
C                   AV has an observable eigenvalue on the unit circle;
C             = 6:  the reduction to a real Schur form of the state
C                   matrix of a minimal realization of W failed;
C             = 7:  a failure was detected during the ordering of the
C                   real Schur form of the state matrix of a minimal
C                   realization of W or in the iterative process to
C                   compute a right coprime factorization with inner
C                   denominator;
C             = 8:  if DICO = 'C' and the matrix AW has a controllable
C                   eigenvalue on the imaginary axis, or DICO = 'D' and
C                   AW has a controllable eigenvalue on the unit circle;
C             = 9:  the computation of eigenvalues failed;
C             = 10: the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let G be the transfer-function matrix of the original
C     linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t),                          (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09ID determines
C     the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                      (2)
C
C     such that the corresponding transfer-function matrix Gr minimizes
C     the norm of the frequency-weighted error
C
C             V*(G-Gr)*W,                                    (3)
C
C     where V and W are transfer-function matrices without poles on the
C     imaginary axis in continuous-time case or on the unit circle in
C     discrete-time case.
C
C     The following procedure is used to reduce G:
C
C     1) Decompose additively G, of order N, as
C
C          G = G1 + G2,
C
C        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and
C        G2 = (A2,B2,C2,0), of order NU, has only ALPHA-unstable poles.
C
C     2) Compute for G1 a B&T or SPA frequency-weighted approximation
C        G1r of order NR-NU using the combination method or the
C        modified combination method of [4].
C
C     3) Assemble the reduced model Gr as
C
C           Gr = G1r + G2.
C
C     For the frequency-weighted reduction of the ALPHA-stable part,
C     several methods described in [4] can be employed in conjunction
C     with the combination method and modified combination method
C     proposed in [4].
C
C     If JOB = 'B', the square-root B&T method is used.
C     If JOB = 'F', the balancing-free square-root version of the
C     B&T method is used.
C     If JOB = 'S', the square-root version of the SPA method is used.
C     If JOB = 'P', the balancing-free square-root version of the
C     SPA method is used.
C
C     For each of these methods, left and right truncation matrices
C     are determined using the Cholesky factors of an input
C     frequency-weighted controllability Grammian P and an output
C     frequency-weighted observability Grammian Q.
C     P and Q are computed from the controllability Grammian Pi of G*W
C     and the observability Grammian Qo of V*G. Using special
C     realizations of G*W and V*G, Pi and Qo are computed in the
C     partitioned forms
C
C           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) ,
C                ( P12' P22 )               ( Q12' Q22 )
C
C     where P11 and Q11 are the leading N-by-N parts of Pi and Qo,
C     respectively. Let P0 and Q0 be non-negative definite matrices
C     defined below
C                                        -1
C            P0 = P11 - ALPHAC**2*P12*P22 *P21 ,
C                                        -1
C            Q0 = Q11 - ALPHAO**2*Q12*Q22 *Q21.
C
C     The frequency-weighted controllability and observability
C     Grammians, P and Q, respectively, are defined as follows:
C     P = P0 if JOBC = 'S' (standard combination method [4]);
C     P = P1 >= P0 if JOBC = 'E', where P1 is the controllability
C     Grammian defined to enforce stability for a modified combination
C     method of [4];
C     Q = Q0 if JOBO = 'S' (standard combination method [4]);
C     Q = Q1 >= Q0 if JOBO = 'E', where Q1 is the observability
C     Grammian defined to enforce stability for a modified combination
C     method of [4].
C
C     If JOBC = JOBO = 'S' and ALPHAC = ALPHAO = 0, the choice of
C     Grammians corresponds to the method of Enns [1], while if
C     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds
C     to the method of Lin and Chiu [2,3].
C
C     If JOBC = 'S' and ALPHAC = 1, no pole-zero cancellations must
C     occur in G*W. If JOBO = 'S' and ALPHAO = 1, no pole-zero
C     cancellations must occur in V*G. The presence of pole-zero
C     cancellations leads to meaningless results and must be avoided.
C
C     The frequency-weighted Hankel singular values HSV(1), ....,
C     HSV(N) are computed as the square roots of the eigenvalues
C     of the product P*Q.
C
C     REFERENCES
C
C     [1] Enns, D.
C         Model reduction with balanced realizations: An error bound
C         and a frequency weighted generalization.
C         Proc. 23-th CDC, Las Vegas, pp. 127-132, 1984.
C
C     [2] Lin, C.-A. and Chiu, T.-Y.
C         Model reduction via frequency-weighted balanced realization.
C         Control Theory and Advanced Technology, vol. 8,
C         pp. 341-351, 1992.
C
C     [3] Sreeram, V., Anderson, B.D.O and Madievski, A.G.
C         New results on frequency weighted balanced reduction
C         technique.
C         Proc. ACC, Seattle, Washington, pp. 4004-4009, 1995.
C
C     [4] Varga, A. and Anderson, B.D.O.
C         Square-root balancing-free methods for the frequency-weighted
C         balancing related model reduction.
C         (report in preparation)
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root
C     techniques.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000.
C     D. Sima, University of Bucharest, August 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C              Sep. 2001.
C
C     KEYWORDS
C
C     Frequency weighting, model reduction, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  C100, ONE, ZERO
      PARAMETER         ( C100 = 100.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOB, JOBC, JOBO, ORDSEL, WEIGHT
      INTEGER           INFO, IWARN, LDA, LDAV, LDAW, LDB, LDBV, LDBW,
     $                  LDC, LDCV, LDCW, LDD, LDDV, LDDW, LDWORK, M, MW,
     $                  N, NR, NS, NV, NW, P, PV
      DOUBLE PRECISION  ALPHA, ALPHAC, ALPHAO, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AV(LDAV,*), AW(LDAW,*),
     $                  B(LDB,*), BV(LDBV,*), BW(LDBW,*),
     $                  C(LDC,*), CV(LDCV,*), CW(LDCW,*),
     $                  D(LDD,*), DV(LDDV,*), DW(LDDW,*), DWORK(*),
     $                  HSV(*)
C     .. Local Scalars ..
      LOGICAL           BAL, BTA, DISCR, FIXORD, FRWGHT, LEFTW, RIGHTW,
     $                  SCALE, SPA
      INTEGER           IERR, IWARNL, KBR, KBV, KBW, KCR, KCV, KCW, KDR,
     $                  KDV, KI, KL, KT, KTI, KU, KW, LCF, LDW, LW, NMR,
     $                  NN, NNQ, NNR, NNV, NNW, NRA, NU, NU1, NVR, NWR,
     $                  PPV, WRKOPT
      DOUBLE PRECISION  ALPWRK, MAXRED, SCALEC, SCALEO
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09IX, AB09IY, DLACPY, SB08CD, SB08DD, TB01ID,
     $                  TB01KD, TB01PD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      BTA    = LSAME( JOB,    'B' ) .OR. LSAME( JOB, 'F' )
      SPA    = LSAME( JOB,    'S' ) .OR. LSAME( JOB, 'P' )
      BAL    = LSAME( JOB,    'B' ) .OR. LSAME( JOB, 'S' )
      SCALE  = LSAME( EQUIL,  'S' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LEFTW  = LSAME( WEIGHT, 'L' ) .OR. LSAME( WEIGHT, 'B' )
      RIGHTW = LSAME( WEIGHT, 'R' ) .OR. LSAME( WEIGHT, 'B' )
      FRWGHT = LEFTW .OR. RIGHTW
C
      LW  = 1
      NN  = N*N
      NNV = N + NV
      NNW = N + NW
      PPV = MAX( P, PV )
C
      IF( LEFTW .AND. PV.GT.0 ) THEN
         LW = MAX( LW, NNV*( NNV + MAX( NNV, PV ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( P + 5 ) )
      END IF
C
      IF( RIGHTW .AND. MW.GT.0 ) THEN
         LW = MAX( LW, NNW*( NNW + MAX( NNW, MW ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( M + 5 ) )
      END IF
      LW = 2*NN + MAX( LW, 2*NN + 5*N, N*MAX( M, P ) )
C
      IF( LEFTW .AND. NV.GT.0 ) THEN
         LCF = PV*( NV + PV ) + PV*NV +
     $         MAX( NV*( NV + 5 ), PV*( PV + 2 ), 4*PPV )
         IF( PV.EQ.P ) THEN
            LW = MAX( LW, LCF, NV + MAX( NV, 3*P ) )
         ELSE
            LW = MAX( LW, PPV*( 2*NV + PPV ) +
     $                    MAX( LCF, NV + MAX( NV, 3*PPV ) ) )
         END IF
      END IF
C
      IF( RIGHTW .AND. NW.GT.0 ) THEN
         IF( MW.EQ.M ) THEN
            LW = MAX( LW, NW + MAX( NW, 3*M ) )
         ELSE
            LW = MAX( LW, 2*NW*MAX( M, MW ) +
     $                    NW + MAX( NW, 3*M, 3*MW ) )
         END IF
         LW = MAX( LW, MW*( NW + MW ) +
     $             MAX( NW*( NW + 5 ), MW*( MW + 2 ), 4*MW, 4*M ) )
      END IF
C
C     Check the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( JOBC, 'S' ) .OR. LSAME( JOBC, 'E' ) ) )
     $     THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( JOBO, 'S' ) .OR. LSAME( JOBO, 'E' ) ) )
     $     THEN
         INFO = -3
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( .NOT. ( SCALE  .OR. LSAME( EQUIL,  'N' ) ) ) THEN
         INFO = -6
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -7
      ELSE IF( N.LT.0 ) THEN
         INFO = -8
      ELSE IF( M.LT.0 ) THEN
         INFO = -9
      ELSE IF( P.LT.0 ) THEN
         INFO = -10
      ELSE IF( NV.LT.0 ) THEN
         INFO = -11
      ELSE IF( PV.LT.0 ) THEN
         INFO = -12
      ELSE IF( NW.LT.0 ) THEN
         INFO = -13
      ELSE IF( MW.LT.0 ) THEN
         INFO = -14
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -15
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -16
      ELSE IF( ABS( ALPHAC ).GT.ONE  ) THEN
         INFO = -17
      ELSE IF( ABS( ALPHAO ).GT.ONE  ) THEN
         INFO = -18
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -20
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -22
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -24
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -26
      ELSE IF( LDAV.LT.1 .OR. ( LEFTW  .AND. LDAV.LT.NV ) ) THEN
         INFO = -28
      ELSE IF( LDBV.LT.1 .OR. ( LEFTW  .AND. LDBV.LT.NV ) ) THEN
         INFO = -30
      ELSE IF( LDCV.LT.1 .OR. ( LEFTW  .AND. LDCV.LT.PV ) ) THEN
         INFO = -32
      ELSE IF( LDDV.LT.1 .OR. ( LEFTW  .AND. LDDV.LT.PV ) ) THEN
         INFO = -34
      ELSE IF( LDAW.LT.1 .OR. ( RIGHTW .AND. LDAW.LT.NW ) ) THEN
         INFO = -36
      ELSE IF( LDBW.LT.1 .OR. ( RIGHTW .AND. LDBW.LT.NW ) ) THEN
         INFO = -38
      ELSE IF( LDCW.LT.1 .OR. ( RIGHTW .AND. LDCW.LT.M  ) ) THEN
         INFO = -40
      ELSE IF( LDDW.LT.1 .OR. ( RIGHTW .AND. LDDW.LT.M  ) ) THEN
         INFO = -42
      ELSE IF( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 ) THEN
         INFO = -46
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -49
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09ID', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         NS = 0
         IWORK(1) = 0
         IWORK(2) = NV
         IWORK(3) = NW
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( SCALE ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C        Workspace: N.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
      END IF
C
C     Correct the value of ALPHA to ensure stability.
C
      ALPWRK = ALPHA
      IF( DISCR ) THEN
         IF( ALPHA.EQ.ONE ) ALPWRK = ONE - SQRT( DLAMCH( 'E' ) )
      ELSE
         IF( ALPHA.EQ.ZERO ) ALPWRK = -SQRT( DLAMCH( 'E' ) )
      END IF
C
C     Allocate working storage.
C
      KU = 1
      KL = KU + NN
      KI = KL + N
      KW = KI + N
C
C     Reduce A to a block-diagonal real Schur form, with the
C     ALPHA-unstable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation, A <- inv(T)*A*T, and
C     apply the transformation to B and C: B <- inv(T)*B and C <- C*T.
C
C     Workspace needed:      N*(N+2);
C     Additional workspace:  need   3*N;
C                            prefer larger.
C
      CALL TB01KD( DICO, 'Unstable', 'General', N, M, P, ALPWRK, A, LDA,
     $             B, LDB, C, LDC, NU, DWORK(KU), N, DWORK(KL),
     $             DWORK(KI), DWORK(KW), LDWORK-KW+1, IERR )
C
      IF( IERR.NE.0 ) THEN
         IF( IERR.NE.3 ) THEN
            INFO = 1
         ELSE
            INFO = 2
         END IF
         RETURN
      END IF
C
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Determine NRA, the desired order for the reduction of stable part.
C
      IWARNL = 0
      NS = N - NU
      IF( FIXORD ) THEN
         NRA = MAX( 0, NR-NU )
         IF( NR.LT.NU )
     $      IWARNL = 3
      ELSE
         NRA = 0
      END IF
C
C     Finish if only unstable part is present.
C
      IF( NS.EQ.0 ) THEN
         NR = NU
         DWORK(1) = WRKOPT
         IWORK(1) = 0
         IWORK(2) = NV
         IWORK(3) = NW
         RETURN
      END IF
C
      NVR = NV
      IF( LEFTW .AND. NV.GT.0 ) THEN
C
C        Compute a left-coprime factorization with inner denominator
C        of a minimal realization of V. The resulting AV is in
C        real Schur form.
C        Workspace needed:   real  LV+MAX( 1, LCF,
C                                          NV + MAX( NV, 3*P, 3*PV ) ),
C                                  where
C                                  LV = 0 if P = PV and
C                                  LV = MAX(P,PV)*(2*NV+MAX(P,PV))
C                                         otherwise;
C                                  LCF = PV*(NV+PV) +
C                                        MAX( 1, PV*NV + MAX( NV*(NV+5),
C                                             PV*(PV+2),4*PV,4*P ) );
C                                  prefer larger;
C                          integer NV + MAX(P,PV).
C
         IF( P.EQ.PV ) THEN
            KW = 1
            CALL TB01PD( 'Minimal', 'Scale', NV, P, PV, AV, LDAV,
     $                   BV, LDBV, CV, LDCV, NVR, ZERO,
     $                   IWORK, DWORK, LDWORK, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
            KBR = 1
            KDR = KBR + PV*NVR
            KW  = KDR + PV*PV
            CALL SB08CD( DICO, NVR, P, PV, AV, LDAV, BV, LDBV, CV, LDCV,
     $                   DV, LDDV, NNQ, NNR, DWORK(KBR), MAX( 1, NVR ),
     $                   DWORK(KDR), PV, ZERO, DWORK(KW), LDWORK-KW+1,
     $                   IWARN, IERR )
         ELSE
            LDW = MAX( P, PV )
            KBV = 1
            KCV = KBV + NV*LDW
            KW  = KCV + NV*LDW
            CALL DLACPY( 'Full', NV, P, BV, LDBV, DWORK(KBV), NV )
            CALL DLACPY( 'Full', PV, NV, CV, LDCV, DWORK(KCV), LDW )
            CALL TB01PD( 'Minimal', 'Scale', NV, P, PV, AV, LDAV,
     $                   DWORK(KBV), NV, DWORK(KCV), LDW, NVR, ZERO,
     $                   IWORK, DWORK(KW), LDWORK-KW+1, INFO )
            KDV = KW
            KBR = KDV + LDW*LDW
            KDR = KBR + PV*NVR
            KW  = KDR + PV*PV
            CALL DLACPY( 'Full', PV, P, DV, LDDV, DWORK(KDV), LDW )
            CALL SB08CD( DICO, NVR, P, PV, AV, LDAV, DWORK(KBV), NV,
     $                   DWORK(KCV), LDW, DWORK(KDV), LDW, NNQ, NNR,
     $                   DWORK(KBR), MAX( 1, NVR ), DWORK(KDR), PV,
     $                   ZERO, DWORK(KW), LDWORK-KW+1, IWARN, IERR )
            CALL DLACPY( 'Full', NVR, P, DWORK(KBV), NV, BV, LDBV )
            CALL DLACPY( 'Full', PV, NVR, DWORK(KCV), LDW, CV, LDCV )
            CALL DLACPY( 'Full', PV, P, DWORK(KDV), LDW, DV, LDDV )
         END IF
         IF( IERR.NE.0 ) THEN
            INFO = IERR + 2
            RETURN
         END IF
         NVR = NNQ
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         IF( IWARN.GT.0 )
     $      IWARN = 10 + IWARN
      END IF
C
      NWR = NW
      IF( RIGHTW .AND. NW.GT.0 ) THEN
C
C        Compute a minimal realization of W.
C        Workspace needed:   real  LW+MAX(1, NW + MAX(NW, 3*M, 3*MW));
C                                  where
C                                  LW = 0,              if M = MW and
C                                  LW = 2*NW*MAX(M,MW), otherwise;
C                                  prefer larger;
C                          integer NW + MAX(M,MW).
C
         IF( M.EQ.MW ) THEN
            KW = 1
            CALL TB01PD( 'Minimal', 'Scale', NW, MW, M, AW, LDAW,
     $                   BW, LDBW, CW, LDCW, NWR, ZERO, IWORK, DWORK,
     $                   LDWORK, INFO )
         ELSE
            LDW = MAX( M, MW )
            KBW = 1
            KCW = KBW + NW*LDW
            KW  = KCW + NW*LDW
            CALL DLACPY( 'Full', NW, MW, BW, LDBW, DWORK(KBW), NW )
            CALL DLACPY( 'Full', M, NW, CW, LDCW, DWORK(KCW), LDW )
            CALL TB01PD( 'Minimal', 'Scale', NW, MW, M, AW, LDAW,
     $                   DWORK(KBW), NW, DWORK(KCW), LDW, NWR, ZERO,
     $                   IWORK, DWORK(KW), LDWORK-KW+1, INFO )
            CALL DLACPY( 'Full', NWR, MW, DWORK(KBW), NW, BW, LDBW )
            CALL DLACPY( 'Full', M, NWR, DWORK(KCW), LDW, CW, LDCW )
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      END IF
C
      IF( RIGHTW .AND. NWR.GT.0 ) THEN
C
C        Compute a right-coprime factorization with inner denominator
C        of the minimal realization of W. The resulting AW is in
C        real Schur form.
C
C        Workspace needed:  MW*(NW+MW) +
C                           MAX( 1, NW*(NW+5), MW*(MW+2), 4*MW, 4*M );
C                           prefer larger.
C
         LDW = MAX( 1, MW )
         KCR = 1
         KDR = KCR + NWR*LDW
         KW  = KDR + MW*LDW
         CALL SB08DD( DICO, NWR, MW, M, AW, LDAW, BW, LDBW, CW, LDCW,
     $                DW, LDDW, NNQ, NNR, DWORK(KCR), LDW, DWORK(KDR),
     $                LDW, ZERO, DWORK(KW), LDWORK-KW+1, IWARN, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = IERR + 5
            RETURN
         END IF
         NWR = NNQ
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         IF( IWARN.GT.0 )
     $      IWARN = 10 + IWARN
      END IF
C
      NU1 = NU + 1
C
C     Allocate working storage.
C
      KT  = 1
      KTI = KT  + NN
      KW  = KTI + NN
C
C     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R
C     of the controllability and observability Grammians, respectively.
C     Real workspace:    need  2*N*N + MAX( 1, LLEFT, LRIGHT ),
C             where
C             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5)
C                              if WEIGHT = 'L' or 'B' and PV > 0;
C             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0;
C             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5)
C                              if WEIGHT = 'R' or 'B' and MW > 0;
C             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0.
C                        prefer larger.
C
      CALL AB09IY( DICO, JOBC, JOBO, WEIGHT, NS, M, P, NVR, PV, NWR,
     $             MW, ALPHAC, ALPHAO, A(NU1,NU1), LDA, B(NU1,1), LDB,
     $             C(1,NU1), LDC, AV, LDAV, BV, LDBV, CV, LDCV,
     $             DV, LDDV, AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $             SCALEC, SCALEO, DWORK(KTI), N, DWORK(KT), N,
     $             DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 9
         RETURN
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Compute a BTA or SPA of the stable part.
C     Real workspace:  need  2*N*N + MAX( 1, 2*N*N+5*N, N*MAX(M,P) ).
C
      CALL AB09IX( DICO, JOB, 'Schur', ORDSEL, NS, M, P, NRA,
     $             SCALEC, SCALEO, A(NU1,NU1), LDA, B(NU1,1), LDB,
     $             C(1,NU1), LDC, D, LDD, DWORK(KTI), N, DWORK(KT), N,
     $             NMR, HSV, TOL1, TOL2, IWORK, DWORK(KW), LDWORK-KW+1,
     $             IWARN, IERR )
      IWARN = MAX( IWARN, IWARNL )
      IF( IERR.NE.0 ) THEN
         INFO = 10
         RETURN
      END IF
      NR = NRA + NU
C
      DWORK(1) = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      IWORK(1) = NMR
      IWORK(2) = NVR
      IWORK(3) = NWR
C
      RETURN
C *** Last line of AB09ID ***
      END
