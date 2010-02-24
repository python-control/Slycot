      SUBROUTINE AB09JD( JOBV, JOBW, JOBINV, DICO, EQUIL, ORDSEL,
     $                   N, NV, NW, M, P, NR, ALPHA, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, AV, LDAV, BV, LDBV,
     $                   CV, LDCV, DV, LDDV, AW, LDAW, BW, LDBW,
     $                   CW, LDCW, DW, LDDW, NS, HSV, TOL1, TOL2,
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for an original
C     state-space representation (A,B,C,D) by using the frequency
C     weighted optimal Hankel-norm approximation method.
C     The Hankel norm of the weighted error
C
C           op(V)*(G-Gr)*op(W)
C
C     is minimized, where G and Gr are the transfer-function matrices
C     of the original and reduced systems, respectively, V and W are
C     invertible transfer-function matrices representing the left and
C     right frequency weights, and op(X) denotes X, inv(X), conj(X) or
C     conj(inv(X)). V and W are specified by their state space
C     realizations (AV,BV,CV,DV) and (AW,BW,CW,DW), respectively.
C     When minimizing ||V*(G-Gr)*W||, V and W must be antistable.
C     When minimizing inv(V)*(G-Gr)*inv(W), V and W must have only
C     antistable zeros.
C     When minimizing conj(V)*(G-Gr)*conj(W), V and W must be stable.
C     When minimizing conj(inv(V))*(G-Gr)*conj(inv(W)), V and W must
C     be minimum-phase.
C     If the original system is unstable, then the frequency weighted
C     Hankel-norm approximation is computed only for the
C     ALPHA-stable part of the system.
C
C     For a transfer-function matrix G, conj(G) denotes the conjugate
C     of G given by G'(-s) for a continuous-time system or G'(1/z)
C     for a discrete-time system.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBV    CHARACTER*1
C             Specifies the left frequency-weighting as follows:
C             = 'N':  V = I;
C             = 'V':  op(V) = V;
C             = 'I':  op(V) = inv(V);
C             = 'C':  op(V) = conj(V);
C             = 'R':  op(V) = conj(inv(V)).
C
C     JOBW    CHARACTER*1
C             Specifies the right frequency-weighting as follows:
C             = 'N':  W = I;
C             = 'W':  op(W) = W;
C             = 'I':  op(W) = inv(W);
C             = 'C':  op(W) = conj(W);
C             = 'R':  op(W) = conj(inv(W)).
C
C     JOBINV  CHARACTER*1
C             Specifies the computational approach to be used as
C             follows:
C             = 'N':  use the inverse free descriptor system approach;
C             = 'I':  use the inversion based standard approach;
C             = 'A':  switch automatically to the inverse free
C                     descriptor approach in case of badly conditioned
C                     feedthrough matrices in V or W (see METHOD).
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
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
C     NV      (input) INTEGER
C             The order of the realization of the left frequency
C             weighting V, i.e., the order of the matrix AV.  NV >= 0.
C
C     NW      (input) INTEGER
C             The order of the realization of the right frequency
C             weighting W, i.e., the order of the matrix AW.  NW >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of
C             the resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. For a system with NU ALPHA-unstable
C             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N),
C             NR is set as follows: if ORDSEL = 'F', NR is equal to
C             NU+MIN(MAX(0,NR-NU-KR+1),NMIN), where KR is the
C             multiplicity of the Hankel singular value HSV(NR-NU+1),
C             NR is the desired order on entry, and NMIN is the order
C             of a minimal realization of the ALPHA-stable part of the
C             given system; NMIN is determined as the number of Hankel
C             singular values greater than NS*EPS*HNORM(As,Bs,Cs), where
C             EPS is the machine precision (see LAPACK Library Routine
C             DLAMCH) and HNORM(As,Bs,Cs) is the Hankel norm of the
C             ALPHA-stable part of the weighted system (computed in
C             HSV(1));
C             if ORDSEL = 'A', NR is the sum of NU and the number of
C             Hankel singular values greater than
C             MAX(TOL1,NS*EPS*HNORM(As,Bs,Cs)).
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
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the
C             reduced order system in a real Schur form.
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
C             On entry, if JOBV <> 'N', the leading NV-by-NV part of
C             this array must contain the state matrix AV of a state
C             space realization of the left frequency weighting V.
C             On exit, if JOBV <> 'N', and INFO = 0, the leading
C             NV-by-NV part of this array contains the real Schur form
C             of AV.
C             AV is not referenced if JOBV = 'N'.
C
C     LDAV    INTEGER
C             The leading dimension of the array AV.
C             LDAV >= MAX(1,NV), if JOBV <> 'N';
C             LDAV >= 1,         if JOBV =  'N'.
C
C     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P)
C             On entry, if JOBV <> 'N', the leading NV-by-P part of
C             this array must contain the input matrix BV of a state
C             space realization of the left frequency weighting V.
C             On exit, if JOBV <> 'N', and INFO = 0, the leading
C             NV-by-P part of this array contains the transformed
C             input matrix BV corresponding to the transformed AV.
C             BV is not referenced if JOBV = 'N'.
C
C     LDBV    INTEGER
C             The leading dimension of the array BV.
C             LDBV >= MAX(1,NV), if JOBV <> 'N';
C             LDBV >= 1,         if JOBV =  'N'.
C
C     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV)
C             On entry, if JOBV <> 'N', the leading P-by-NV part of
C             this array must contain the output matrix CV of a state
C             space realization of the left frequency weighting V.
C             On exit, if JOBV <> 'N', and INFO = 0, the leading
C             P-by-NV part of this array contains the transformed output
C             matrix CV corresponding to the transformed AV.
C             CV is not referenced if JOBV = 'N'.
C
C     LDCV    INTEGER
C             The leading dimension of the array CV.
C             LDCV >= MAX(1,P), if JOBV <> 'N';
C             LDCV >= 1,        if JOBV =  'N'.
C
C     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P)
C             If JOBV <> 'N', the leading P-by-P part of this array
C             must contain the feedthrough matrix DV of a state space
C             realization of the left frequency weighting V.
C             DV is not referenced if JOBV = 'N'.
C
C     LDDV    INTEGER
C             The leading dimension of the array DV.
C             LDDV >= MAX(1,P), if JOBV <> 'N';
C             LDDV >= 1,        if JOBV =  'N'.
C
C     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW)
C             On entry, if JOBW <> 'N', the leading NW-by-NW part of
C             this array must contain the state matrix AW of a state
C             space realization of the right frequency weighting W.
C             On exit, if JOBW <> 'N', and INFO = 0, the leading
C             NW-by-NW part of this array contains the real Schur form
C             of AW.
C             AW is not referenced if JOBW = 'N'.
C
C     LDAW    INTEGER
C             The leading dimension of the array AW.
C             LDAW >= MAX(1,NW), if JOBW <> 'N';
C             LDAW >= 1,         if JOBW =  'N'.
C
C     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M)
C             On entry, if JOBW <> 'N', the leading NW-by-M part of
C             this array must contain the input matrix BW of a state
C             space realization of the right frequency weighting W.
C             On exit, if JOBW <> 'N', and INFO = 0, the leading
C             NW-by-M part of this array contains the transformed
C             input matrix BW corresponding to the transformed AW.
C             BW is not referenced if JOBW = 'N'.
C
C     LDBW    INTEGER
C             The leading dimension of the array BW.
C             LDBW >= MAX(1,NW), if JOBW <> 'N';
C             LDBW >= 1,         if JOBW =  'N'.
C
C     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW)
C             On entry, if JOBW <> 'N', the leading M-by-NW part of
C             this array must contain the output matrix CW of a state
C             space realization of the right frequency weighting W.
C             On exit, if JOBW <> 'N', and INFO = 0, the leading
C             M-by-NW part of this array contains the transformed output
C             matrix CW corresponding to the transformed AW.
C             CW is not referenced if JOBW = 'N'.
C
C     LDCW    INTEGER
C             The leading dimension of the array CW.
C             LDCW >= MAX(1,M), if JOBW <> 'N';
C             LDCW >= 1,        if JOBW =  'N'.
C
C     DW      (input) DOUBLE PRECISION array, dimension (LDDW,M)
C             If JOBW <> 'N', the leading M-by-M part of this array
C             must contain the feedthrough matrix DW of a state space
C             realization of the right frequency weighting W.
C             DW is not referenced if JOBW = 'N'.
C
C     LDDW    INTEGER
C             The leading dimension of the array DW.
C             LDDW >= MAX(1,M), if JOBW <> 'N';
C             LDDW >= 1,        if JOBW =  'N'.
C
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of this array contain
C             the Hankel singular values, ordered decreasingly, of the
C             projection G1s of op(V)*G1*op(W) (see METHOD), where G1
C             is the ALPHA-stable part of the original system.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(G1s), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(G1s) is the
C             Hankel-norm of the projection G1s of op(V)*G1*op(W)
C             (see METHOD), computed in HSV(1).
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = NS*EPS*HNORM(G1s), where NS is the number of
C             ALPHA-stable eigenvalues of A and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C             TOL1 < 1.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the ALPHA-stable part of the given system.
C             The recommended value is TOL2 = NS*EPS*HNORM(G1s).
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C             TOL2 < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = MAX(1,M,c,d),    if DICO = 'C',
C             LIWORK = MAX(1,N,M,c,d),  if DICO = 'D', where
C                c = 0,                          if JOBV =  'N',
C                c = MAX(2*P,NV+P+N+6,2*NV+P+2), if JOBV <> 'N',
C                d = 0,                          if JOBW =  'N',
C                d = MAX(2*M,NW+M+N+6,2*NW+M+2), if JOBW <> 'N'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( LDW1, LDW2, LDW3, LDW4 ), where
C             for NVP = NV+P and NWM = NW+M we have
C             LDW1 = 0 if JOBV =  'N' and
C             LDW1 = 2*NVP*(NVP+P) + P*P +
C                    MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ),
C                          NVP*N + MAX( NVP*N+N*N, P*N, P*M ) )
C                      if JOBV <> 'N',
C             LDW2 = 0 if JOBW =  'N' and
C             LDW2 = 2*NWM*(NWM+M) + M*M +
C                    MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ),
C                          NWM*N + MAX( NWM*N+N*N, M*N, P*M ) )
C                      if JOBW <> 'N',
C             LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2,
C             LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) +
C                    MAX( 3*M+1, MIN(N,M)+P ).
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
C                   system. In this case, the resulting NR is set equal
C                   to NSMIN.
C             = 2:  with ORDSEL = 'F', the selected order NR is less
C                   than the order of the ALPHA-unstable part of the
C                   given system. In this case NR is set equal to the
C                   order of the ALPHA-unstable part.
C
C     Error Indicator
C
C     INFO    INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  the computation of the ordered real Schur form of A
C                    failed;
C             =  2:  the separation of the ALPHA-stable/unstable
C                    diagonal blocks failed because of very close
C                    eigenvalues;
C             =  3:  the reduction of AV to a real Schur form failed;
C             =  4:  the reduction of AW to a real Schur form failed;
C             =  5:  the reduction to generalized Schur form of the
C                    descriptor pair corresponding to the inverse of V
C                    failed;
C             =  6:  the reduction to generalized Schur form of the
C                    descriptor pair corresponding to the inverse of W
C                    failed;
C             =  7:  the computation of Hankel singular values failed;
C             =  8:  the computation of stable projection in the
C                    Hankel-norm approximation algorithm failed;
C             =  9:  the order of computed stable projection in the
C                    Hankel-norm approximation algorithm differs
C                    from the order of Hankel-norm approximation;
C             = 10:  the reduction of AV-BV*inv(DV)*CV to a
C                    real Schur form failed;
C             = 11:  the reduction of AW-BW*inv(DW)*CW to a
C                    real Schur form failed;
C             = 12:  the solution of the Sylvester equation failed
C                    because the poles of V (if JOBV = 'V') or of
C                    conj(V) (if JOBV = 'C') are not distinct from
C                    the poles of G1 (see METHOD);
C             = 13:  the solution of the Sylvester equation failed
C                    because the poles of W (if JOBW = 'W') or of
C                    conj(W) (if JOBW = 'C') are not distinct from
C                    the poles of G1 (see METHOD);
C             = 14:  the solution of the Sylvester equation failed
C                    because the zeros of V (if JOBV = 'I') or of
C                    conj(V) (if JOBV = 'R') are not distinct from
C                    the poles of G1sr (see METHOD);
C             = 15:  the solution of the Sylvester equation failed
C                    because the zeros of W (if JOBW = 'I') or of
C                    conj(W) (if JOBW = 'R') are not distinct from
C                    the poles of G1sr (see METHOD);
C             = 16:  the solution of the generalized Sylvester system
C                    failed because the zeros of V (if JOBV = 'I') or
C                    of conj(V) (if JOBV = 'R') are not distinct from
C                    the poles of G1sr (see METHOD);
C             = 17:  the solution of the generalized Sylvester system
C                    failed because the zeros of W (if JOBW = 'I') or
C                    of conj(W) (if JOBW = 'R') are not distinct from
C                    the poles of G1sr (see METHOD);
C             = 18:  op(V) is not antistable;
C             = 19:  op(W) is not antistable;
C             = 20:  V is not invertible;
C             = 21:  W is not invertible.
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
C     for a discrete-time system. The subroutine AB09JD determines
C     the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                      (2)
C
C     such that the corresponding transfer-function matrix Gr minimizes
C     the Hankel-norm of the frequency-weighted error
C
C             op(V)*(G-Gr)*op(W).                            (3)
C
C     For minimizing (3) with op(V) = V and op(W) = W, V and W are
C     assumed to have poles distinct from those of G, while with
C     op(V) = conj(V) and op(W) = conj(W), conj(V) and conj(W) are
C     assumed to have poles distinct from those of G. For minimizing (3)
C     with op(V) = inv(V) and op(W) = inv(W), V and W are assumed to
C     have zeros distinct from the poles of G, while with
C     op(V) = conj(inv(V)) and op(W) = conj(inv(W)), conj(V) and conj(W)
C     are assumed to have zeros distinct from the poles of G.
C
C     Note: conj(G) = G'(-s) for a continuous-time system and
C           conj(G) = G'(1/z) for a discrete-time system.
C
C     The following procedure is used to reduce G (see [1]):
C
C     1) Decompose additively G as
C
C          G = G1 + G2,
C
C        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and
C        G2 = (A2,B2,C2,0) has only ALPHA-unstable poles.
C
C     2) Compute G1s, the projection of op(V)*G1*op(W) containing the
C        poles of G1, using explicit formulas [4] or the inverse-free
C        descriptor system formulas of [5].
C
C     3) Determine G1sr, the optimal Hankel-norm approximation of G1s,
C        of order r.
C
C     4) Compute G1r, the projection of inv(op(V))*G1sr*inv(op(W))
C        containing the poles of G1sr, using explicit formulas [4]
C        or the inverse-free descriptor system formulas of [5].
C
C     5) Assemble the reduced model Gr as
C
C           Gr = G1r + G2.
C
C     To reduce the weighted ALPHA-stable part G1s at step 3, the
C     optimal Hankel-norm approximation method of [2], based on the
C     square-root balancing projection formulas of [3], is employed.
C
C     The optimal weighted approximation error satisfies
C
C          HNORM[op(V)*(G-Gr)*op(W)] >= S(r+1),
C
C     where S(r+1) is the (r+1)-th Hankel singular value of G1s, the
C     transfer-function matrix computed at step 2 of the above
C     procedure, and HNORM(.) denotes the Hankel-norm.
C
C     REFERENCES
C
C     [1] Latham, G.A. and Anderson, B.D.O.
C         Frequency-weighted optimal Hankel-norm approximation of stable
C         transfer functions.
C         Systems & Control Letters, Vol. 5, pp. 229-236, 1985.
C
C     [2] Glover, K.
C         All optimal Hankel norm approximation of linear
C         multivariable systems and their L-infinity error bounds.
C         Int. J. Control, Vol. 36, pp. 1145-1193, 1984.
C
C     [3] Tombs, M.S. and Postlethwaite, I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [4] Varga, A.
C         Explicit formulas for an efficient implementation
C         of the frequency-weighting model reduction approach.
C         Proc. 1993 European Control Conference, Groningen, NL,
C         pp. 693-696, 1993.
C
C     [5] Varga, A.
C         Efficient and numerically reliable implementation of the
C         frequency-weighted Hankel-norm approximation model reduction
C         approach.
C         Proc. 2001 ECC, Porto, Portugal, 2001.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on an accuracy enhancing square-root
C     technique.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2001.
C     D. Sima, University of Bucharest, April 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
C
C     REVISIONS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C     March 2005.
C
C     KEYWORDS
C
C     Frequency weighting, model reduction, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  C100, ONE, P0001, ZERO
      PARAMETER         ( C100 = 100.0D0, ONE = 1.0D0, P0001 = 0.0001D0,
     $                    ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOBINV, JOBV, JOBW, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDAV, LDAW, LDB, LDBV, LDBW,
     $                  LDC, LDCV, LDCW, LDD, LDDV, LDDW, LDWORK, M, N,
     $                  NR, NS, NV, NW, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AV(LDAV,*), AW(LDAW,*),
     $                  B(LDB,*), BV(LDBV,*), BW(LDBW,*),
     $                  C(LDC,*), CV(LDCV,*), CW(LDCW,*),
     $                  D(LDD,*), DV(LDDV,*), DW(LDDW,*), DWORK(*),
     $                  HSV(*)
C     .. Local Scalars ..
      CHARACTER         JOBVL, JOBWL
      LOGICAL           AUTOM, CONJV, CONJW, DISCR, FIXORD, FRWGHT,
     $                  INVFR, LEFTI, LEFTW, RIGHTI, RIGHTW
      INTEGER           IERR, IWARNL, KAV, KAW, KBV, KBW, KCV, KCW, KDV,
     $                  KDW, KEV, KEW, KI, KL, KU, KW, LDABV, LDABW,
     $                  LDCDV, LDCDW, LW, NRA, NU, NU1, NVP, NWM, RANK
      DOUBLE PRECISION  ALPWRK, MAXRED, RCOND, SQREPS, TOL, WRKOPT
C     .. Local Arrays ..
      DOUBLE PRECISION  TEMP(1)
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07ND, AB08MD, AB09CX, AB09JV, AB09JW, AG07BD,
     $                  DLACPY, TB01ID, TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LEFTI  = LSAME( JOBV, 'I' ) .OR. LSAME( JOBV, 'R' )
      LEFTW  = LSAME( JOBV, 'V' ) .OR. LSAME( JOBV, 'C' ) .OR. LEFTI
      CONJV  = LSAME( JOBV, 'C' ) .OR. LSAME( JOBV, 'R' )
      RIGHTI = LSAME( JOBW, 'I' ) .OR. LSAME( JOBW, 'R' )
      RIGHTW = LSAME( JOBW, 'W' ) .OR. LSAME( JOBW, 'C' ) .OR. RIGHTI
      CONJW  = LSAME( JOBW, 'C' ) .OR. LSAME( JOBW, 'R' )
      FRWGHT = LEFTW .OR. RIGHTW
      INVFR  = LSAME( JOBINV, 'N' )
      AUTOM  = LSAME( JOBINV, 'A' )
C
      LW = 1
      IF( LEFTW ) THEN
         NVP = NV + P
         LW  = MAX( LW, 2*NVP*( NVP + P ) + P*P +
     $              MAX( 2*NVP*NVP + MAX( 11*NVP + 16, P*NVP ),
     $                   NVP*N + MAX( NVP*N+N*N, P*N, P*M ) ) )
      END IF
      IF( RIGHTW ) THEN
         NWM = NW + M
         LW  = MAX( LW, 2*NWM*( NWM + M ) + M*M +
     $              MAX( 2*NWM*NWM + MAX( 11*NWM + 16, M*NWM ),
     $                   NWM*N + MAX( NWM*N+N*N, M*N, P*M ) ) )
      END IF
      LW = MAX( LW, N*( 2*N + MAX( N, M, P ) + 5 ) + ( N*( N + 1 ) )/2 )
      LW = MAX( LW, N*( M + P + 2 ) + 2*M*P + MIN( N, M ) +
     $                             MAX ( 3*M + 1, MIN( N, M ) + P ) )
C
C     Check the input scalar arguments.
C
      IF( .NOT. ( LSAME( JOBV, 'N' ) .OR. LEFTW ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( JOBW, 'N' ) .OR. RIGHTW ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( INVFR .OR. AUTOM .OR. LSAME( JOBINV, 'I' ) ) )
     $   THEN
         INFO = -3
      ELSE IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( NV.LT.0 ) THEN
         INFO = -8
      ELSE IF( NW.LT.0 ) THEN
         INFO = -9
      ELSE IF( M.LT.0 ) THEN
         INFO = -10
      ELSE IF( P.LT.0 ) THEN
         INFO = -11
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -12
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -13
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -19
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -21
      ELSE IF( LDAV.LT.1 .OR. ( LEFTW  .AND. LDAV.LT.NV ) ) THEN
         INFO = -23
      ELSE IF( LDBV.LT.1 .OR. ( LEFTW  .AND. LDBV.LT.NV ) ) THEN
         INFO = -25
      ELSE IF( LDCV.LT.1 .OR. ( LEFTW  .AND. LDCV.LT.P  ) ) THEN
         INFO = -27
      ELSE IF( LDDV.LT.1 .OR. ( LEFTW  .AND. LDDV.LT.P  ) ) THEN
         INFO = -29
      ELSE IF( LDAW.LT.1 .OR. ( RIGHTW .AND. LDAW.LT.NW ) ) THEN
         INFO = -31
      ELSE IF( LDBW.LT.1 .OR. ( RIGHTW .AND. LDBW.LT.NW ) ) THEN
         INFO = -33
      ELSE IF( LDCW.LT.1 .OR. ( RIGHTW .AND. LDCW.LT.M  ) ) THEN
         INFO = -35
      ELSE IF( LDDW.LT.1 .OR. ( RIGHTW .AND. LDDW.LT.M  ) ) THEN
         INFO = -37
      ELSE IF( TOL1.GE.ONE ) THEN
         INFO = -40
      ELSE IF( ( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 )
     $      .OR. TOL2.GE.ONE ) THEN
         INFO = -41
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -44
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09JD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         NS = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D,  B <- inv(D)*B  and  C <- C*D,  where D is a
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
      SQREPS = SQRT( DLAMCH( 'E' ) )
      IF( DISCR ) THEN
         IF( ALPHA.EQ.ONE ) ALPWRK = ONE - SQREPS
      ELSE
         IF( ALPHA.EQ.ZERO ) ALPWRK = -SQREPS
      END IF
C
C     Allocate working storage.
C
      KU = 1
      KL = KU + N*N
      KI = KL + N
      KW = KI + N
C
C     Compute an additive decomposition G = G1 + G2, where G1
C     is the ALPHA-stable projection of G.
C
C     Reduce A to a block-diagonal real Schur form, with the NU-th order
C     ALPHA-unstable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation A <- inv(T)*A*T and
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
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
      IWARNL = 0
C
      NS = N - NU
      IF( FIXORD ) THEN
         NRA = MAX( 0, NR-NU )
         IF( NR.LT.NU )
     $      IWARNL = 2
      ELSE
         NRA = 0
      END IF
C
C     Finish if only unstable part is present.
C
      IF( NS.EQ.0 ) THEN
         NR = NU
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      NU1 = NU + 1
      IF( CONJV ) THEN
         JOBVL = 'C'
      ELSE
         JOBVL = 'V'
      END IF
      IF( CONJW ) THEN
         JOBWL = 'C'
      ELSE
         JOBWL = 'W'
      END IF
      IF( LEFTW ) THEN
C
C        Check if V is invertible.
C        Real workspace:    need   (NV+P)**2 + MAX( P + MAX(3*P,NV),
C                                  MIN(P+1,NV) + MAX(3*(P+1),NV+P) );
C                           prefer larger.
C        Integer workspace: need   2*NV+P+2.
C
         TOL = ZERO
         CALL AB08MD( 'S', NV, P, P, AV, LDAV, BV, LDBV, CV, LDCV,
     $                DV, LDDV, RANK, TOL, IWORK, DWORK, LDWORK,
     $                IERR )
         IF( RANK.NE.P ) THEN
            INFO = 20
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, DWORK(1) )
C
         IF( LEFTI ) THEN
            IF( INVFR ) THEN
               IERR = 1
            ELSE
C
C              Allocate storage for a standard inverse of V.
C              Workspace: need  NV*(NV+2*P) + P*P.
C
               KAV = 1
               KBV = KAV + NV*NV
               KCV = KBV + NV*P
               KDV = KCV + P*NV
               KW  = KDV + P*P
C
               LDABV = MAX( NV, 1 )
               LDCDV = P
               CALL DLACPY( 'Full', NV, NV, AV, LDAV,
     $                      DWORK(KAV), LDABV )
               CALL DLACPY( 'Full', NV, P,  BV, LDBV,
     $                      DWORK(KBV), LDABV )
               CALL DLACPY( 'Full', P,  NV, CV, LDCV,
     $                      DWORK(KCV), LDCDV )
               CALL DLACPY( 'Full', P,  P,  DV, LDDV,
     $                      DWORK(KDV), LDCDV )
C
C              Compute the standard inverse of V.
C              Additional real workspace:   need   MAX(1,4*P);
C                                           prefer larger.
C              Integer workspace:           need   2*P.
C
               CALL AB07ND( NV, P, DWORK(KAV), LDABV, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV,
     $                      RCOND, IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
C              Check if inversion is accurate.
C
               IF( AUTOM ) THEN
                  IF( IERR.EQ.0 .AND. RCOND.LE.P0001  ) IERR = 1
               ELSE
                  IF( IERR.EQ.0 .AND. RCOND.LE.SQREPS ) IERR = 1
               END IF
               IF( IERR.NE.0 .AND. NV.EQ.0 ) THEN
                  INFO = 20
                  RETURN
               END IF
            END IF
C
            IF( IERR.NE.0 ) THEN
C
C              Allocate storage for a descriptor inverse of V.
C
               KAV = 1
               KEV = KAV + NVP*NVP
               KBV = KEV + NVP*NVP
               KCV = KBV + NVP*P
               KDV = KCV + P*NVP
               KW  = KDV + P*P
C
               LDABV = MAX( NVP, 1 )
               LDCDV = P
C
C              DV is singular or ill-conditioned.
C              Form a descriptor inverse of V.
C              Workspace: need  2*(NV+P)*(NV+2*P) + P*P.
C
               CALL AG07BD( 'I', NV, P, AV, LDAV, TEMP, 1, BV, LDBV,
     $                      CV, LDCV, DV, LDDV, DWORK(KAV), LDABV,
     $                      DWORK(KEV), LDABV, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV, IERR )
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using descriptor inverse of V
C              of order NVP = NV + P.
C              Additional real workspace: need
C                 MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ),
C                      NVP*N + MAX( NVP*N+N*N, P*N, P*M ) );
C                 prefer larger.
C              Integer workspace: need NVP+N+6.
C
               CALL AB09JV( JOBVL, DICO, 'G', 'C', NS, M, P, NVP, P,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD,
     $                      DWORK(KAV), LDABV, DWORK(KEV), LDABV,
     $                      DWORK(KBV), LDABV, DWORK(KCV), LDCDV,
     $                      DWORK(KDV), LDCDV, IWORK, DWORK(KW),
     $                      LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 5
                  ELSE IF( IERR.EQ.2 ) THEN
                     INFO = 16
                  ELSE IF( IERR.EQ.4 ) THEN
                     INFO = 18
                  END IF
                  RETURN
               END IF
            ELSE
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using explicit inverse of V.
C              Additional real workspace: need
C                 MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                      a = 0,    if DICO = 'C' or  JOBVL = 'V',
C                      a = 2*NV, if DICO = 'D' and JOBVL = 'C';
C                 prefer larger.
C
               CALL AB09JV( JOBVL, DICO, 'I', 'C', NS, M, P, NV, P,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAV), LDABV,
     $                      TEMP, 1, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV, IWORK,
     $                      DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 10
                  ELSE IF( IERR.EQ.3 ) THEN
                     INFO = 14
                  ELSE IF( IERR.EQ.4 ) THEN
                     INFO = 18
                  END IF
                  RETURN
               END IF
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW - 1 ) )
         ELSE
C
C           Compute the projection of V*G1 or conj(V)*G1 containing the
C           poles of G.
C
C           Workspace need:
C           real   MAX( 1, NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                       a = 0,    if DICO = 'C' or  JOBVL = 'V',
C                       a = 2*NV, if DICO = 'D' and JOBVL = 'C';
C           prefer larger.
C
            CALL AB09JV( JOBVL, DICO, 'I', 'C', NS, M, P, NV, P,
     $                   A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                   C(1,NU1), LDC, D, LDD, AV, LDAV,
     $                   TEMP, 1, BV, LDBV, CV, LDCV, DV, LDDV, IWORK,
     $                   DWORK, LDWORK, IERR )
            IF( IERR.NE.0 ) THEN
               IF( IERR.EQ.1 ) THEN
                  INFO = 3
               ELSE IF( IERR.EQ.3 ) THEN
                  INFO = 12
               ELSE IF( IERR.EQ.4 ) THEN
                  INFO = 18
               END IF
               RETURN
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(1) )
         END IF
      END IF
C
      IF( RIGHTW ) THEN
C
C        Check if W is invertible.
C        Real workspace:    need   (NW+M)**2 + MAX( M + MAX(3*M,NW),
C                                  MIN(M+1,NW) + MAX(3*(M+1),NW+M) );
C                           prefer larger.
C        Integer workspace: need   2*NW+M+2.
C
         TOL = ZERO
         CALL AB08MD( 'S', NW, M, M, AW, LDAW, BW, LDBW, CW, LDCW,
     $                DW, LDDW, RANK, TOL, IWORK, DWORK, LDWORK,
     $                IERR )
         IF( RANK.NE.M ) THEN
            INFO = 21
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, DWORK(1) )
C
         IF( RIGHTI ) THEN
            IF( INVFR ) THEN
               IERR = 1
            ELSE
C
C              Allocate storage for a standard inverse of W.
C              Workspace: need  NW*(NW+2*M) + M*M.
C
               KAW = 1
               KBW = KAW + NW*NW
               KCW = KBW + NW*M
               KDW = KCW + M*NW
               KW  = KDW + M*M
C
               LDABW = MAX( NW, 1 )
               LDCDW = M
               CALL DLACPY( 'Full', NW, NW, AW, LDAW,
     $                      DWORK(KAW), LDABW )
               CALL DLACPY( 'Full', NW, M,  BW, LDBW,
     $                      DWORK(KBW), LDABW )
               CALL DLACPY( 'Full', M,  NW, CW, LDCW,
     $                      DWORK(KCW), LDCDW )
               CALL DLACPY( 'Full', M,  M,  DW, LDDW,
     $                      DWORK(KDW), LDCDW )
C
C              Compute the standard inverse of W.
C              Additional real workspace:   need   MAX(1,4*M);
C                                           prefer larger.
C              Integer workspace:           need   2*M.
C
               CALL AB07ND( NW, M, DWORK(KAW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      RCOND, IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
C              Check if inversion is accurate.
C
               IF( AUTOM ) THEN
                  IF( IERR.EQ.0 .AND. RCOND.LE.P0001  ) IERR = 1
               ELSE
                  IF( IERR.EQ.0 .AND. RCOND.LE.SQREPS ) IERR = 1
               END IF
               IF( IERR.NE.0 .AND. NW.EQ.0 ) THEN
                  INFO = 21
                  RETURN
               END IF
            END IF
C
            IF( IERR.NE.0 ) THEN
C
C              Allocate storage for a descriptor inverse of W.
C
               KAW = 1
               KEW = KAW + NWM*NWM
               KBW = KEW + NWM*NWM
               KCW = KBW + NWM*M
               KDW = KCW + M*NWM
               KW  = KDW + M*M
C
               LDABW = MAX( NWM, 1 )
               LDCDW = M
C
C              DW is singular or ill-conditioned.
C              Form the descriptor inverse of W.
C              Workspace: need  2*(NW+M)*(NW+2*M) + M*M.
C
               CALL AG07BD( 'I', NW, M, AW, LDAW, TEMP, 1, BW, LDBW,
     $                      CW, LDCW, DW, LDDW, DWORK(KAW), LDABW,
     $                      DWORK(KEW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW, IERR )
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using descriptor inverse of W
C              of order NWM = NW + M.
C              Additional real workspace: need
C                 MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ),
C                      NWM*N + MAX( NWM*N+N*N, M*N, P*M ) );
C                 prefer larger.
C              Integer workspace: need NWM+N+6.
C
               CALL AB09JW( JOBWL, DICO, 'G', 'C', NS, M, P, NWM, M,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAW), LDABW,
     $                      DWORK(KEW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 6
                  ELSE IF( IERR.EQ.2 ) THEN
                     INFO = 17
                  ELSE IF( IERR.EQ.4 ) THEN
                     INFO = 19
                  END IF
                  RETURN
               END IF
            ELSE
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using explicit inverse of W.
C              Additional real workspace: need
C                 MAX( NW*(NW+5), NW*N + MAX( a, M*N, P*M ) )
C                      a = 0,    if DICO = 'C' or  JOBWL = 'W',
C                      a = 2*NW, if DICO = 'D' and JOBWL = 'C';
C                 prefer larger.
C
               CALL AB09JW( JOBWL, DICO, 'I', 'C', NS, M, P, NW, M,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAW), LDABW,
     $                      TEMP, 1, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 11
                  ELSE IF( IERR.EQ.3 ) THEN
                     INFO = 15
                  ELSE IF( IERR.EQ.4 ) THEN
                     INFO = 19
                  END IF
                  RETURN
               END IF
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW - 1 ) )
         ELSE
C
C           Compute the projection G1s of V*G1*W or conj(V)*G1*conj(W)
C           containing the poles of G.
C
C           Workspace need:
C           real   MAX( 1, NW*(NW+5), NW*N + MAX( b, M*N, P*M ) )
C                    b = 0,    if DICO = 'C' or  JOBWL = 'W',
C                    b = 2*NW, if DICO = 'D' and JOBWL = 'C';
C           prefer larger.
C
            CALL AB09JW( JOBWL, DICO, 'I', 'C', NS, M, P, NW, M,
     $                   A(NU1,NU1), LDA, B(NU1,1), LDB, C(1,NU1), LDC,
     $                   D, LDD, AW, LDAW, TEMP, 1, BW, LDBW, CW, LDCW,
     $                   DW, LDDW, IWORK, DWORK, LDWORK, IERR )
            IF( IERR.NE.0 ) THEN
               IF( IERR.EQ.1 ) THEN
                  INFO = 4
               ELSE IF( IERR.EQ.3 ) THEN
                  INFO = 13
               ELSE IF( IERR.EQ.4 ) THEN
                  INFO = 19
               END IF
               RETURN
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(1) )
         END IF
      END IF
C
C     Determine a reduced order approximation G1sr of G1s using the
C     Hankel-norm approximation method. The resulting A(NU1:N,NU1:N)
C     is further in a real Schur form.
C
C     Workspace: need   MAX( LDW3, LDW4 ),
C                LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2,
C                LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) +
C                       MAX( 3*M+1, MIN(N,M)+P );
C                prefer larger.
C
      CALL AB09CX( DICO, ORDSEL, NS, M, P, NRA, A(NU1,NU1), LDA,
     $             B(NU1,1), LDB, C(1,NU1), LDC, D, LDD, HSV, TOL1,
     $             TOL2, IWORK, DWORK, LDWORK, IWARN, IERR )
C
      IF( IERR.NE.0 ) THEN
C
C        Set INFO = 7, 8 or 9.
C
         INFO = IERR + 5
         RETURN
      END IF
C
      IWARN  = MAX( IWARNL, IWARN )
      WRKOPT = MAX( WRKOPT, DWORK(1) )
C
      IF( LEFTW ) THEN
         IF( .NOT.LEFTI ) THEN
            IF( INVFR ) THEN
               IERR = 1
            ELSE
C
C              Allocate storage for a standard inverse of V.
C              Workspace: need  NV*(NV+2*P) + P*P.
C
               KAV = 1
               KBV = KAV + NV*NV
               KCV = KBV + NV*P
               KDV = KCV + P*NV
               KW  = KDV + P*P
C
               LDABV = MAX( NV, 1 )
               LDCDV = P
               CALL DLACPY( 'Full', NV, NV, AV, LDAV,
     $                      DWORK(KAV), LDABV )
               CALL DLACPY( 'Full', NV, P,  BV, LDBV,
     $                      DWORK(KBV), LDABV )
               CALL DLACPY( 'Full', P,  NV, CV, LDCV,
     $                      DWORK(KCV), LDCDV )
               CALL DLACPY( 'Full', P,  P,  DV, LDDV,
     $                      DWORK(KDV), LDCDV )
C
C              Compute the standard inverse of V.
C              Additional real workspace:   need   MAX(1,4*P);
C                                           prefer larger.
C              Integer workspace:           need   2*P.
C
               CALL AB07ND( NV, P, DWORK(KAV), LDABV, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV,
     $                      RCOND, IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
C              Check if inversion is accurate.
C
               IF( AUTOM ) THEN
                  IF( IERR.EQ.0 .AND. RCOND.LE.P0001  ) IERR = 1
               ELSE
                  IF( IERR.EQ.0 .AND. RCOND.LE.SQREPS ) IERR = 1
               END IF
               IF( IERR.NE.0 .AND. NV.EQ.0 ) THEN
                  INFO = 20
                  RETURN
               END IF
            END IF
C
            IF( IERR.NE.0 ) THEN
C
C              Allocate storage for a descriptor inverse of V.
C
               KAV = 1
               KEV = KAV + NVP*NVP
               KBV = KEV + NVP*NVP
               KCV = KBV + NVP*P
               KDV = KCV + P*NVP
               KW  = KDV + P*P
C
               LDABV = MAX( NVP, 1 )
               LDCDV = P
C
C              DV is singular or ill-conditioned.
C              Form a descriptor inverse of V.
C              Workspace: need  2*(NV+P)*(NV+2*P) + P*P.
C
               CALL AG07BD( 'I', NV, P, AV, LDAV, TEMP, 1, BV, LDBV,
     $                      CV, LDCV, DV, LDDV, DWORK(KAV), LDABV,
     $                      DWORK(KEV), LDABV, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV, IERR )
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using descriptor inverse of V
C              of order NVP = NV + P.
C              Additional real workspace: need
C                 MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ),
C                      NVP*N + MAX( NVP*N+N*N, P*N, P*M ) );
C                 prefer larger.
C              Integer workspace: need NVP+N+6.
C
               CALL AB09JV( JOBVL, DICO, 'G', 'N', NRA, M, P, NVP, P,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD,
     $                      DWORK(KAV), LDABV, DWORK(KEV), LDABV,
     $                      DWORK(KBV), LDABV, DWORK(KCV), LDCDV,
     $                      DWORK(KDV), LDCDV, IWORK, DWORK(KW),
     $                      LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 5
                  ELSE IF( IERR.EQ.2 ) THEN
                     INFO = 16
                  END IF
                  RETURN
               END IF
            ELSE
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using explicit inverse of V.
C              Additional real workspace: need
C                 MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                      a = 0,    if DICO = 'C' or  JOBVL = 'V',
C                      a = 2*NV, if DICO = 'D' and JOBVL = 'C';
C                 prefer larger.
C
               CALL AB09JV( JOBVL, DICO, 'I', 'N', NRA, M, P, NV, P,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAV), LDABV,
     $                      TEMP, 1, DWORK(KBV), LDABV,
     $                      DWORK(KCV), LDCDV, DWORK(KDV), LDCDV, IWORK,
     $                      DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 10
                  ELSE IF( IERR.EQ.3 ) THEN
                     INFO = 14
                  END IF
                  RETURN
               END IF
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW - 1 ) )
         ELSE
C
C           Compute the projection of V*G1sr or conj(V)*G1sr containing
C           the poles of G.
C
C           Workspace need:
C           real    MAX( 1, NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                        a = 0,    if DICO = 'C' or  JOBVL = 'V',
C                        a = 2*NV, if DICO = 'D' and JOBVL = 'C';
C           prefer larger.
C
            CALL AB09JV( JOBVL, DICO, 'I', 'N', NRA, M, P, NV, P,
     $                   A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                   C(1,NU1), LDC, D, LDD, AV, LDAV,
     $                   TEMP, 1, BV, LDBV, CV, LDCV, DV, LDDV, IWORK,
     $                   DWORK, LDWORK, IERR )
            IF( IERR.NE.0 ) THEN
               IF( IERR.EQ.1 ) THEN
                  INFO = 3
               ELSE IF( IERR.EQ.3 ) THEN
                  INFO = 12
               END IF
               RETURN
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(1) )
         END IF
      END IF
C
      IF( RIGHTW ) THEN
         IF( .NOT.RIGHTI ) THEN
            IF( INVFR ) THEN
               IERR = 1
            ELSE
C
C              Allocate storage for a standard inverse of W.
C              Workspace: need  NW*(NW+2*M) + M*M.
C
               KAW = 1
               KBW = KAW + NW*NW
               KCW = KBW + NW*M
               KDW = KCW + M*NW
               KW  = KDW + M*M
C
               LDABW = MAX( NW, 1 )
               LDCDW = M
               CALL DLACPY( 'Full', NW, NW, AW, LDAW,
     $                      DWORK(KAW), LDABW )
               CALL DLACPY( 'Full', NW, M,  BW, LDBW,
     $                      DWORK(KBW), LDABW )
               CALL DLACPY( 'Full', M,  NW, CW, LDCW,
     $                      DWORK(KCW), LDCDW )
               CALL DLACPY( 'Full', M,  M,  DW, LDDW,
     $                      DWORK(KDW), LDCDW )
C
C              Compute the standard inverse of W.
C              Additional real workspace:   need   MAX(1,4*M);
C                                           prefer larger.
C              Integer workspace:           need   2*M.
C
               CALL AB07ND( NW, M, DWORK(KAW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      RCOND, IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
C              Check if inversion is accurate.
C
               IF( AUTOM ) THEN
                  IF( IERR.EQ.0 .AND. RCOND.LE.P0001  ) IERR = 1
               ELSE
                  IF( IERR.EQ.0 .AND. RCOND.LE.SQREPS ) IERR = 1
               END IF
               IF( IERR.NE.0 .AND. NW.EQ.0 ) THEN
                  INFO = 21
                  RETURN
               END IF
            END IF
C
            IF( IERR.NE.0 ) THEN
C
C              Allocate storage for a descriptor inverse of W.
C
               KAW = 1
               KEW = KAW + NWM*NWM
               KBW = KEW + NWM*NWM
               KCW = KBW + NWM*M
               KDW = KCW + M*NWM
               KW  = KDW + M*M
C
               LDABW = MAX( NWM, 1 )
               LDCDW = M
C
C              DW is singular or ill-conditioned.
C              Form the descriptor inverse of W.
C              Workspace: need  2*(NW+M)*(NW+2*M) + M*M.
C
               CALL AG07BD( 'I', NW, M, AW, LDAW, TEMP, 1, BW, LDBW,
     $                      CW, LDCW, DW, LDDW, DWORK(KAW), LDABW,
     $                      DWORK(KEW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW, IERR )
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using descriptor inverse of W
C              of order NWM = NW + M.
C              Additional real workspace: need
C                 MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ),
C                      NWM*N + MAX( NWM*N+N*N, M*N, P*M ) );
C                 prefer larger.
C              Integer workspace: need NWM+N+6.
C
               CALL AB09JW( JOBWL, DICO, 'G', 'N', NRA, M, P, NWM, M,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAW), LDABW,
     $                      DWORK(KEW), LDABW, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 6
                  ELSE IF( IERR.EQ.2 ) THEN
                     INFO = 17
                  END IF
                  RETURN
               END IF
            ELSE
C
C              Compute the projection containing the poles of weighted
C              reduced ALPHA-stable part using explicit inverse of W.
C              Additional real workspace: need
C                 MAX( NW*(NW+5), NW*N + MAX( a, M*N, P*M ) )
C                      a = 0,    if DICO = 'C' or  JOBWL = 'W',
C                      a = 2*NW, if DICO = 'D' and JOBWL = 'C';
C                 prefer larger.
C
               CALL AB09JW( JOBWL, DICO, 'I', 'N', NRA, M, P, NW, M,
     $                      A(NU1,NU1), LDA, B(NU1,1), LDB,
     $                      C(1,NU1), LDC, D, LDD, DWORK(KAW), LDABW,
     $                      TEMP, 1, DWORK(KBW), LDABW,
     $                      DWORK(KCW), LDCDW, DWORK(KDW), LDCDW,
     $                      IWORK, DWORK(KW), LDWORK-KW+1, IERR )
               IF( IERR.NE.0 ) THEN
                  IF( IERR.EQ.1 ) THEN
                     INFO = 11
                  ELSE IF( IERR.EQ.3 ) THEN
                     INFO = 15
                  END IF
                  RETURN
               END IF
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW - 1 ) )
         ELSE
C
C           Compute the projection G1r of V*G1sr*W or
C           conj(V)*G1sr*conj(W) containing the poles of G.
C
C           Workspace need:
C           real   MAX( 1, NW*(NW+5), NW*N + MAX( b, M*N, P*M ) )
C                    b = 0,    if DICO = 'C' or  JOBWL = 'W',
C                    b = 2*NW, if DICO = 'D' and JOBWL = 'C';
C           prefer larger.
C
            CALL AB09JW( JOBWL, DICO, 'I', 'N', NRA, M, P, NW, M,
     $                   A(NU1,NU1), LDA, B(NU1,1), LDB, C(1,NU1), LDC,
     $                   D, LDD, AW, LDAW, TEMP, 1, BW, LDBW, CW, LDCW,
     $                   DW, LDDW, IWORK, DWORK, LDWORK, IERR )
C
            IF( IERR.NE.0 ) THEN
               IF( IERR.EQ.1 ) THEN
                  INFO = 4
               ELSE IF( IERR.EQ.3 ) THEN
                  INFO = 13
               END IF
               RETURN
            END IF
C
            WRKOPT = MAX( WRKOPT, DWORK(1) )
         END IF
      END IF
C
      NR = NRA + NU
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09JD ***
      END
