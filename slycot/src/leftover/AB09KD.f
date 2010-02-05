      SUBROUTINE AB09KD( JOB, DICO, WEIGHT, EQUIL, ORDSEL, N, NV, NW, M,
     $                   P, NR, ALPHA, A, LDA, B, LDB, C, LDC, D, LDD,
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
C     weighted optimal Hankel-norm approximation method.
C     The Hankel norm of the weighted error
C
C           V*(G-Gr)*W    or    conj(V)*(G-Gr)*conj(W)
C
C     is minimized, where G and Gr are the transfer-function matrices
C     of the original and reduced systems, respectively, and V and W
C     are the transfer-function matrices of the left and right frequency
C     weights, specified by their state space realizations (AV,BV,CV,DV)
C     and (AW,BW,CW,DW), respectively. When minimizing the weighted
C     error V*(G-Gr)*W, V and W must be antistable transfer-function
C     matrices. When minimizing conj(V)*(G-Gr)*conj(W), V and W must be
C     stable transfer-function matrices.
C     Additionally, V and W must be invertible transfer-function
C     matrices, with the feedthrough matrices DV and DW invertible.
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
C     JOB     CHARACTER*1
C             Specifies the frequency-weighting problem as follows:
C             = 'N':  solve min||V*(G-Gr)*W||_H;
C             = 'C':  solve min||conj(V)*(G-Gr)*conj(W)||_H.
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
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
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV
C             part of this array must contain the state matrix AV of a
C             state space realization of the left frequency weighting V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             NV-by-NV part of this array contains a real Schur form
C             of the state matrix of a state space realization of the
C             inverse of V.
C             AV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDAV    INTEGER
C             The leading dimension of the array AV.
C             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDAV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P)
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part
C             of this array must contain the input matrix BV of a state
C             space realization of the left frequency weighting V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             NV-by-P part of this array contains the input matrix of a
C             state space realization of the inverse of V.
C             BV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDBV    INTEGER
C             The leading dimension of the array BV.
C             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDBV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV)
C             On entry, if WEIGHT = 'L' or 'B', the leading P-by-NV part
C             of this array must contain the output matrix CV of a state
C             space realization of the left frequency weighting V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             P-by-NV part of this array contains the output matrix of a
C             state space realization of the inverse of V.
C             CV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDCV    INTEGER
C             The leading dimension of the array CV.
C             LDCV >= MAX(1,P), if WEIGHT = 'L' or 'B';
C             LDCV >= 1,        if WEIGHT = 'R' or 'N'.
C
C     DV      (input/output) DOUBLE PRECISION array, dimension (LDDV,P)
C             On entry, if WEIGHT = 'L' or 'B', the leading P-by-P part
C             of this array must contain the feedthrough matrix DV of a
C             state space realization of the left frequency weighting V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             P-by-P part of this array contains the feedthrough matrix
C             of a state space realization of the inverse of V.
C             DV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDDV    INTEGER
C             The leading dimension of the array DV.
C             LDDV >= MAX(1,P), if WEIGHT = 'L' or 'B';
C             LDDV >= 1,        if WEIGHT = 'R' or 'N'.
C
C     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW
C             part of this array must contain the state matrix AW of
C             a state space realization of the right frequency
C             weighting W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             NW-by-NW part of this array contains a real Schur form of
C             the state matrix of a state space realization of the
C             inverse of W.
C             AW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDAW    INTEGER
C             The leading dimension of the array AW.
C             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDAW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-M part
C             of this array must contain the input matrix BW of a state
C             space realization of the right frequency weighting W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             NW-by-M part of this array contains the input matrix of a
C             state space realization of the inverse of W.
C             BW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDBW    INTEGER
C             The leading dimension of the array BW.
C             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDBW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part
C             of this array must contain the output matrix CW of a state
C             space realization of the right frequency weighting W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             M-by-NW part of this array contains the output matrix of a
C             state space realization of the inverse of W.
C             CW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDCW    INTEGER
C             The leading dimension of the array CW.
C             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDCW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     DW      (input/output) DOUBLE PRECISION array, dimension (LDDW,M)
C             On entry, if WEIGHT = 'R' or 'B', the leading M-by-M part
C             of this array must contain the feedthrough matrix DW of
C             a state space realization of the right frequency
C             weighting W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             M-by-M part of this array contains the feedthrough matrix
C             of a state space realization of the inverse of W.
C             DW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDDW    INTEGER
C             The leading dimension of the array DW.
C             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDDW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of this array contain
C             the Hankel singular values, ordered decreasingly, of the
C             ALPHA-stable part of the weighted original system.
C             HSV(1) is the Hankel norm of the ALPHA-stable weighted
C             subsystem.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(As,Bs,Cs), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(As,Bs,Cs) is the
C             Hankel-norm of the ALPHA-stable part of the weighted
C             original system (computed in HSV(1)).
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = NS*EPS*HNORM(As,Bs,Cs), where NS is the number of
C             ALPHA-stable eigenvalues of A and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the ALPHA-stable part of the given system.
C             The recommended value is TOL2 = NS*EPS*HNORM(As,Bs,Cs).
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = MAX(1,M,c),      if DICO = 'C',
C             LIWORK = MAX(1,N,M,c),    if DICO = 'D',
C             where  c = 0,             if WEIGHT = 'N',
C                    c = 2*P,           if WEIGHT = 'L',
C                    c = 2*M,           if WEIGHT = 'R',
C                    c = MAX(2*M,2*P),  if WEIGHT = 'B'.
C             On exit, if INFO = 0, IWORK(1) contains NMIN, the order of
C             the computed minimal realization.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( LDW1, LDW2, LDW3, LDW4 ), where
C             LDW1 = 0 if WEIGHT = 'R' or 'N' and
C             LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                    if WEIGHT = 'L' or WEIGHT = 'B',
C             LDW2 = 0 if WEIGHT = 'L' or 'N' and
C             LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) )
C                    if WEIGHT = 'R' or WEIGHT = 'B', with
C                a = 0,    b = 0,     if DICO = 'C' or  JOB = 'N',
C                a = 2*NV, b = 2*NW,  if DICO = 'D' and JOB = 'C';
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
C                   system; in this case, the resulting NR is set equal
C                   to NSMIN;
C             = 2:  with ORDSEL = 'F', the selected order NR is less
C                   than the order of the ALPHA-unstable part of the
C                   given system; in this case NR is set equal to the
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
C             =  3:  the reduction of AV or AV-BV*inv(DV)*CV to a
C                    real Schur form failed;
C             =  4:  the reduction of AW or AW-BW*inv(DW)*CW to a
C                    real Schur form failed;
C             =  5:  JOB = 'N' and AV is not antistable, or
C                    JOB = 'C' and AV is not stable;
C             =  6:  JOB = 'N' and AW is not antistable, or
C                    JOB = 'C' and AW is not stable;
C             =  7:  the computation of Hankel singular values failed;
C             =  8:  the computation of stable projection in the
C                    Hankel-norm approximation algorithm failed;
C             =  9:  the order of computed stable projection in the
C                    Hankel-norm approximation algorithm differs
C                    from the order of Hankel-norm approximation;
C             = 10:  DV is singular;
C             = 11:  DW is singular;
C             = 12:  the solution of the Sylvester equation failed
C                    because the zeros of V (if JOB = 'N') or of conj(V)
C                    (if JOB = 'C') are not distinct from the poles
C                    of G1sr (see METHOD);
C             = 13:  the solution of the Sylvester equation failed
C                    because the zeros of W (if JOB = 'N') or of conj(W)
C                    (if JOB = 'C') are not distinct from the poles
C                    of G1sr (see METHOD).
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
C     for a discrete-time system. The subroutine AB09KD determines
C     the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                      (2)
C
C     such that the corresponding transfer-function matrix Gr minimizes
C     the Hankel-norm of the frequency-weighted error
C
C             V*(G-Gr)*W,                                    (3)
C     or
C             conj(V)*(G-Gr)*conj(W).                        (4)
C
C     For minimizing (3), V and W are assumed to be antistable, while
C     for minimizing (4), V and W are assumed to be stable transfer-
C     function matrices.
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
C     2) Compute G1s, the stable projection of V*G1*W or
C        conj(V)*G1*conj(W), using explicit formulas [4].
C
C     3) Determine G1sr, the optimal Hankel-norm approximation of G1s
C        of order r.
C
C     4) Compute G1r, the stable projection of either inv(V)*G1sr*inv(W)
C        or conj(inv(V))*G1sr*conj(inv(W)), using explicit formulas [4].
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
C          HNORM[V*(G-Gr)*W] = S(r+1),
C     or
C          HNORM[conj(V)*(G-Gr)*conj(W)] = S(r+1),
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
C     [3] Tombs M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [4] Varga A.
C         Explicit formulas for an efficient implementation
C         of the frequency-weighting model reduction approach.
C         Proc. 1993 European Control Conference, Groningen, NL,
C         pp. 693-696, 1993.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on an accuracy enhancing square-root
C     technique.
C                                         3
C     The algorithms require less than 30N  floating point operations.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, April 2000.
C     D. Sima, University of Bucharest, May 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2000.
C     Based on the RASP routines SFRLW, SFRLW1, SFRRW and SFRRW1,
C     by A. Varga, 1992.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000.
C              Oct. 2001, March 2005.
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
      CHARACTER         DICO, EQUIL, JOB, ORDSEL, WEIGHT
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
      LOGICAL           CONJS, DISCR, FIXORD, FRWGHT, LEFTW, RIGHTW
      INTEGER           IA, IB, IERR, IWARNL, KI, KL, KU, KW, LW, NMIN,
     $                  NRA, NU, NU1
      DOUBLE PRECISION  ALPWRK, MAXRED, RCOND, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07ND, AB09CX, AB09KX, TB01ID, TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      CONJS  = LSAME( JOB,    'C' )
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LEFTW  = LSAME( WEIGHT, 'L' ) .OR. LSAME( WEIGHT, 'B' )
      RIGHTW = LSAME( WEIGHT, 'R' ) .OR. LSAME( WEIGHT, 'B' )
      FRWGHT = LEFTW .OR. RIGHTW
C
      IF ( DISCR .AND. CONJS ) THEN
         IA = 2*NV
         IB = 2*NW
      ELSE
         IA = 0
         IB = 0
      END IF
      LW = 1
      IF( LEFTW )
     $   LW = MAX( LW, NV*(NV+5), NV*N + MAX( IA, P*N, P*M ) )
      IF( RIGHTW )
     $   LW = MAX( LW, MAX( NW*(NW+5), NW*N + MAX( IB, M*N, P*M ) ) )
      LW = MAX( LW, N*( 2*N + MAX( N, M, P ) + 5 ) + ( N*( N + 1 ) )/2 )
      LW = MAX( LW, N*( M + P + 2 ) + 2*M*P + MIN( N, M ) +
     $                             MAX ( 3*M + 1, MIN( N, M ) + P ) )
C
C     Check the input scalar arguments.
C
      IF( .NOT. ( LSAME( JOB, 'N' ) .OR. CONJS ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( NV.LT.0 ) THEN
         INFO = -7
      ELSE IF( NW.LT.0 ) THEN
         INFO = -8
      ELSE IF( M.LT.0 ) THEN
         INFO = -9
      ELSE IF( P.LT.0 ) THEN
         INFO = -10
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -11
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -12
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -18
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -20
      ELSE IF( LDAV.LT.1 .OR. ( LEFTW  .AND. LDAV.LT.NV ) ) THEN
         INFO = -22
      ELSE IF( LDBV.LT.1 .OR. ( LEFTW  .AND. LDBV.LT.NV ) ) THEN
         INFO = -24
      ELSE IF( LDCV.LT.1 .OR. ( LEFTW  .AND. LDCV.LT.P  ) ) THEN
         INFO = -26
      ELSE IF( LDDV.LT.1 .OR. ( LEFTW  .AND. LDDV.LT.P  ) ) THEN
         INFO = -28
      ELSE IF( LDAW.LT.1 .OR. ( RIGHTW .AND. LDAW.LT.NW ) ) THEN
         INFO = -30
      ELSE IF( LDBW.LT.1 .OR. ( RIGHTW .AND. LDBW.LT.NW ) ) THEN
         INFO = -32
      ELSE IF( LDCW.LT.1 .OR. ( RIGHTW .AND. LDCW.LT.M  ) ) THEN
         INFO = -34
      ELSE IF( LDDW.LT.1 .OR. ( RIGHTW .AND. LDDW.LT.M  ) ) THEN
         INFO = -36
      ELSE IF( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 ) THEN
         INFO = -40
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -43
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09KD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         NS = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
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
      KL = KU + N*N
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
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
C     Compute the stable projection of the weighted ALPHA-stable part.
C
C     Workspace: need   MAX( 1, LDW1, LDW2 ),
C                LDW1 = 0 if WEIGHT = 'R' or 'N' and
C                LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                       if WEIGHT = 'L' or 'B',
C                LDW2 = 0 if WEIGHT = 'L' or 'N' and
C                LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) )
C                       if WEIGHT = 'R' or 'B',
C                where  a = 0,    b = 0,    if DICO = 'C' or  JOB = 'N',
C                       a = 2*NV, b = 2*NW, if DICO = 'D' and JOB = 'C';
C                prefer larger.
C
      NS = N - NU
C
C     Finish if only unstable part is present.
C
      IF( NS.EQ.0 ) THEN
         NR = NU
         IWORK(1) = 0
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      NU1 = NU + 1
      IF( FRWGHT ) THEN
         CALL AB09KX( JOB, DICO, WEIGHT, NS, NV, NW, M, P, A(NU1,NU1),
     $                LDA, B(NU1,1), LDB, C(1,NU1), LDC, D, LDD,
     $                AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $                DWORK, LDWORK, IWARNL, IERR )
C
         IF( IERR.NE.0 ) THEN
C
C           Note: Only IERR = 1 or IERR = 2 are possible.
C           Set INFO to 3 or 4.
C
            INFO = IERR + 2
            RETURN
         END IF
C
         IF( IWARNL.NE.0 ) THEN
C
C           Stability/antistability of V and W are compulsory.
C
            IF( IWARNL.EQ.1 .OR. IWARNL.EQ.3 ) THEN
               INFO = 5
            ELSE
               INFO = 6
            END IF
            RETURN
         END IF
C
         DWORK(1) = MAX( WRKOPT, DWORK(1) )
      END IF
C
C     Determine a reduced order approximation of the ALPHA-stable part.
C
C     Workspace: need   MAX( LDW3, LDW4 ),
C                LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2,
C                LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) +
C                       MAX( 3*M+1, MIN(N,M)+P );
C                prefer larger.
C
      IWARNL = 0
      IF( FIXORD ) THEN
         NRA = MAX( 0, NR - NU )
         IF( NRA.EQ.0 )
     $      IWARNL = 2
      ELSE
         NRA = 0
      END IF
      CALL AB09CX( DICO, ORDSEL, NS, M, P, NRA, A(NU1,NU1), LDA,
     $             B(NU1,1), LDB, C(1,NU1), LDC, D, LDD, HSV, TOL1,
     $             TOL2, IWORK, DWORK, LDWORK, IWARN, IERR )
C
      IWARN = MAX( IWARN, IWARNL )
      IF( IERR.NE.0 ) THEN
C
C        Set INFO = 7, 8 or 9.
C
         INFO = IERR + 5
         RETURN
      END IF
C
      WRKOPT = MAX( WRKOPT, DWORK(1) )
      NMIN = IWORK(1)
C
C     Compute the state space realizations of the inverses of V and W.
C
C     Integer workspace: need   c,
C     Real workspace:    need   MAX(1,2*c),
C                        where  c = 0,             if WEIGHT = 'N',
C                               c = 2*P,           if WEIGHT = 'L',
C                               c = 2*M,           if WEIGHT = 'R',
C                               c = MAX(2*M,2*P),  if WEIGHT = 'B'.
C
      IF( LEFTW ) THEN
         CALL AB07ND( NV, P, AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                RCOND, IWORK, DWORK, LDWORK, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 10
            RETURN
         END IF
      END IF
      IF( RIGHTW ) THEN
         CALL AB07ND( NW, M, AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $                RCOND, IWORK, DWORK, LDWORK, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 11
            RETURN
         END IF
      END IF
C
      WRKOPT = MAX( WRKOPT, DWORK(1) )
C
C     Compute the stable projection of weighted reduced ALPHA-stable
C     part.
C
      IF( FRWGHT ) THEN
         CALL AB09KX( JOB, DICO, WEIGHT, NRA, NV, NW, M, P, A(NU1,NU1),
     $                LDA, B(NU1,1), LDB, C(1,NU1), LDC, D, LDD,
     $                AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $                DWORK, LDWORK, IWARNL, IERR )
C
         IF( IERR.NE.0 ) THEN
            IF( IERR.LE.2 ) THEN
C
C              Set INFO to 3 or 4.
C
               INFO = IERR + 2
            ELSE
C
C              Set INFO to 12 or 13.
C
               INFO = IERR + 9
            END IF
            RETURN
         END IF
      END IF
C
      NR = NRA + NU
      IWORK(1) = NMIN
      DWORK(1) = MAX( WRKOPT, DWORK(1) )
C
      RETURN
C *** Last line of AB09KD ***
      END
