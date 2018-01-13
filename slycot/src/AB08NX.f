      SUBROUTINE AB08NX( N, M, P, RO, SIGMA, SVLMAX, ABCD, LDABCD,
     $                   NINFZ, INFZ, KRONL, MU, NU, NKROL, TOL, IWORK,
     $                   DWORK, LDWORK, INFO )
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
C     To extract from the (N+P)-by-(M+N) system
C                  ( B  A )
C                  ( D  C )
C     an (NU+MU)-by-(M+NU) "reduced" system
C                  ( B' A')
C                  ( D' C')
C     having the same transmission zeros but with D' of full row rank.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     RO      (input/output) INTEGER
C             On entry,
C             = P     for the original system;
C             = MAX(P-M, 0) for the pertransposed system.
C             On exit, RO contains the last computed rank.
C
C     SIGMA   (input/output) INTEGER
C             On entry,
C             = 0  for the original system;
C             = M  for the pertransposed system.
C             On exit, SIGMA contains the last computed value sigma in
C             the algorithm.
C
C     SVLMAX  (input) DOUBLE PRECISION
C             During each reduction step, the rank-revealing QR
C             factorization of a matrix stops when the estimated minimum
C             singular value is smaller than TOL * MAX(SVLMAX,EMSV),
C             where EMSV is the estimated maximum singular value.
C             SVLMAX >= 0.
C
C     ABCD    (input/output) DOUBLE PRECISION array, dimension
C             (LDABCD,M+N)
C             On entry, the leading (N+P)-by-(M+N) part of this array
C             must contain the compound input matrix of the system.
C             On exit, the leading (NU+MU)-by-(M+NU) part of this array
C             contains the reduced compound input matrix of the system.
C
C     LDABCD  INTEGER
C             The leading dimension of array ABCD.
C             LDABCD >= MAX(1,N+P).
C
C     NINFZ   (input/output) INTEGER
C             On entry, the currently computed number of infinite zeros.
C             It should be initialized to zero on the first call.
C             NINFZ >= 0.
C             On exit, the number of infinite zeros.
C
C     INFZ    (input/output) INTEGER array, dimension (N)
C             On entry, INFZ(i) must contain the current number of
C             infinite zeros of degree i, where i = 1,2,...,N, found in
C             the previous call(s) of the routine. It should be
C             initialized to zero on the first call.
C             On exit, INFZ(i) contains the number of infinite zeros of
C             degree i, where i = 1,2,...,N.
C
C     KRONL   (input/output) INTEGER array, dimension (N+1)
C             On entry, this array must contain the currently computed
C             left Kronecker (row) indices found in the previous call(s)
C             of the routine. It should be initialized to zero on the
C             first call.
C             On exit, the leading NKROL elements of this array contain
C             the left Kronecker (row) indices.
C
C     MU      (output) INTEGER
C             The normal rank of the transfer function matrix of the
C             original system.
C
C     NU      (output) INTEGER
C             The dimension of the reduced system matrix and the number
C             of (finite) invariant zeros if D' is invertible.
C
C     NKROL   (output) INTEGER
C             The number of left Kronecker indices.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             NOTE that when SVLMAX > 0, the estimated ranks could be
C             less than those defined above (see SVLMAX).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX(M,P))
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
C                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ).
C             For optimum performance LDWORK should be larger.
C
C             If LDWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             DWORK array, returns this value as the first entry of
C             the DWORK array, and no error message related to LDWORK
C             is issued by XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB08BZ by F. Svaricek.
C
C     REVISIONS
C
C     V. Sima, Oct. 1997; Feb. 1998, Jan. 2009, Apr. 2009.
C     A. Varga, May 1999; May 2001.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, multivariable
C     system, orthogonal transformation, structural invariant.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDABCD, LDWORK, M, MU, N, NINFZ, NKROL,
     $                  NU, P, RO, SIGMA
      DOUBLE PRECISION  SVLMAX, TOL
C     .. Array Arguments ..
      INTEGER           INFZ(*), IWORK(*), KRONL(*)
      DOUBLE PRECISION  ABCD(LDABCD,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LQUERY
      INTEGER           I1, IK, IROW, ITAU, IZ, JWORK, MM1, MNTAU, MNU,
     $                  MPM, NB, NP, RANK, RO1, TAU, WRKOPT
      DOUBLE PRECISION  T
C     .. Local Arrays ..
      DOUBLE PRECISION  SVAL(3)
C     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV
C     .. External Subroutines ..
      EXTERNAL          DLAPMT, DLARFG, DLASET, SLCT_DLATZM, DORMQR,
     $                  DORMRQ, MB03OY, MB03PY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      NP   = N + P
      MPM  = MIN( P, M )
      INFO = 0
      LQUERY = ( LDWORK.EQ.-1 )
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( RO.NE.P .AND. RO.NE.MAX( P-M, 0 ) ) THEN
         INFO = -4
      ELSE IF( SIGMA.NE.0 .AND. SIGMA.NE.M ) THEN
         INFO = -5
      ELSE IF( SVLMAX.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDABCD.LT.MAX( 1, NP ) ) THEN
         INFO = -8
      ELSE IF( NINFZ.LT.0 ) THEN
         INFO = -9
      ELSE
         JWORK = MAX( 1,      MPM + MAX( 3*M - 1, N ),
     $                MIN( P, N ) + MAX( 3*P - 1, NP, N+M ) )
         IF( LQUERY ) THEN
            IF( M.GT.0 ) THEN
               NB = MIN( 64, ILAENV( 1, 'DORMQR', 'LT', P, N, MPM,
     $                               -1 ) )
               WRKOPT = MAX( JWORK, MPM + MAX( 1, N )*NB )
            ELSE
               WRKOPT = JWORK
            END IF
            NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RT', NP, N, MIN( P, N ),
     $                            -1 ) )
            WRKOPT = MAX( WRKOPT, MIN( P, N ) + MAX( 1, NP )*NB )
            NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'LN', N, M+N,
     $                            MIN( P, N ), -1 ) )
            WRKOPT = MAX( WRKOPT, MIN( P, N ) + MAX( 1, M+N )*NB )
         ELSE IF( LDWORK.LT.JWORK ) THEN
            INFO = -18
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB08NX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      MU = P
      NU = N
C
      IZ = 0
      IK = 1
      MM1 = M + 1
      ITAU = 1
      NKROL = 0
      WRKOPT = 1
C
C     Main reduction loop:
C
C            M   NU                  M     NU
C      NU  [ B   A ]           NU  [ B     A ]
C      MU  [ D   C ]  -->    SIGMA [ RD   C1 ]   (SIGMA = rank(D) =
C                             TAU  [ 0    C2 ]    row size of RD)
C
C                                    M   NU-RO  RO
C                            NU-RO [ B1   A11  A12 ]
C                     -->      RO  [ B2   A21  A22 ]  (RO = rank(C2) =
C                            SIGMA [ RD   C11  C12 ]   col size of LC)
C                             TAU  [ 0     0   LC  ]
C
C                                      M   NU-RO
C                            NU-RO [ B1   A11 ]     NU := NU - RO
C                                  [----------]     MU := RO + SIGMA
C                     -->      RO  [ B2   A21 ]      D := [B2;RD]
C                            SIGMA [ RD   C11 ]      C := [A21;C11]
C
   20 IF ( MU.EQ.0 )
     $   GO TO 80
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      RO1 = RO
      MNU = M + NU
      IF ( M.GT.0 ) THEN
         IF ( SIGMA.NE.0 ) THEN
            IROW = NU + 1
C
C           Compress rows of D.  First exploit triangular shape.
C           Workspace: need   M+N-1.
C
            DO 40 I1 = 1, SIGMA
               CALL DLARFG( RO+1, ABCD(IROW,I1), ABCD(IROW+1,I1), 1, T )
               CALL SLCT_DLATZM( 'L', RO+1, MNU-I1, ABCD(IROW+1,I1),
     $                          1, T,
     $                      ABCD(IROW,I1+1), ABCD(IROW+1,I1+1), LDABCD,
     $                      DWORK )
               IROW = IROW + 1
   40       CONTINUE
            CALL DLASET( 'Lower', RO+SIGMA-1, SIGMA, ZERO, ZERO,
     $                   ABCD(NU+2,1), LDABCD )
         END IF
C
C        Continue with Householder with column pivoting.
C
C        The rank of D is the number of (estimated) singular values
C        that are greater than TOL * MAX(SVLMAX,EMSV). This number
C        includes the singular values of the first SIGMA columns.
C        Integer workspace: need   M;
C        Workspace: need   min(RO1,M) + 3*M - 1.  RO1 <= P.
C
         IF ( SIGMA.LT.M ) THEN
            JWORK = ITAU + MIN( RO1, M )
            I1    = SIGMA + 1
            IROW  = NU + I1
            CALL MB03OY( RO1, M-SIGMA, ABCD(IROW,I1), LDABCD, TOL,
     $                   SVLMAX, RANK, SVAL, IWORK, DWORK(ITAU),
     $                   DWORK(JWORK), INFO )
            WRKOPT = MAX( WRKOPT, JWORK + 3*M - 2 )
C
C           Apply the column permutations to matrices B and part of D.
C
            CALL DLAPMT( .TRUE., NU+SIGMA, M-SIGMA, ABCD(1,I1), LDABCD,
     $                   IWORK )
C
            IF ( RANK.GT.0 ) THEN
C
C              Apply the Householder transformations to the submatrix C.
C              Workspace: need   min(RO1,M) + NU;
C                         prefer min(RO1,M) + NU*NB.
C
               CALL DORMQR( 'Left', 'Transpose', RO1, NU, RANK,
     $                      ABCD(IROW,I1), LDABCD, DWORK(ITAU),
     $                      ABCD(IROW,MM1), LDABCD, DWORK(JWORK),
     $                      LDWORK-JWORK+1, INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
               IF ( RO1.GT.1 )
     $            CALL DLASET( 'Lower', RO1-1, MIN( RO1-1, RANK ), ZERO,
     $                         ZERO, ABCD(IROW+1,I1), LDABCD )
               RO1 = RO1 - RANK
            END IF
         END IF
      END IF
C
      TAU = RO1
      SIGMA = MU - TAU
C
C     Determination of the orders of the infinite zeros.
C
      IF ( IZ.GT.0 ) THEN
         INFZ(IZ) = INFZ(IZ) + RO - TAU
         NINFZ = NINFZ + IZ*( RO - TAU )
      END IF
      IF ( RO1.EQ.0 )
     $   GO TO 80
      IZ = IZ + 1
C
      IF ( NU.LE.0 ) THEN
         MU = SIGMA
         NU = 0
         RO = 0
      ELSE
C
C        Compress the columns of C2 using RQ factorization with row
C        pivoting, P * C2 = R * Q.
C
         I1 = NU + SIGMA + 1
         MNTAU = MIN( TAU, NU )
         JWORK = ITAU + MNTAU
C
C        The rank of C2 is the number of (estimated) singular values
C        greater than TOL * MAX(SVLMAX,EMSV).
C        Integer Workspace: need TAU;
C        Workspace: need min(TAU,NU) + 3*TAU - 1.
C
         CALL MB03PY( TAU, NU, ABCD(I1,MM1), LDABCD, TOL, SVLMAX, RANK,
     $                SVAL, IWORK, DWORK(ITAU), DWORK(JWORK), INFO )
         WRKOPT = MAX( WRKOPT, JWORK + 3*TAU - 1 )
         IF ( RANK.GT.0 ) THEN
            IROW = I1 + TAU - RANK
C
C           Apply Q' to the first NU columns of [A; C1] from the right.
C           Workspace: need   min(TAU,NU) + NU + SIGMA; SIGMA <= P;
C                      prefer min(TAU,NU) + (NU  + SIGMA)*NB.
C
            CALL DORMRQ( 'Right', 'Transpose', I1-1, NU, RANK,
     $                   ABCD(IROW,MM1), LDABCD, DWORK(MNTAU-RANK+1),
     $                   ABCD(1,MM1), LDABCD, DWORK(JWORK),
     $                   LDWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C           Apply Q to the first NU rows and M + NU columns of [ B  A ]
C           from the left.
C           Workspace: need   min(TAU,NU) + M + NU;
C                      prefer min(TAU,NU) + (M + NU)*NB.
C
            CALL DORMRQ( 'Left', 'NoTranspose', NU, MNU, RANK,
     $                   ABCD(IROW,MM1), LDABCD, DWORK(MNTAU-RANK+1),
     $                   ABCD, LDABCD, DWORK(JWORK), LDWORK-JWORK+1,
     $                   INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
            CALL DLASET( 'Full', RANK, NU-RANK, ZERO, ZERO,
     $                   ABCD(IROW,MM1), LDABCD )
            IF ( RANK.GT.1 )
     $         CALL DLASET( 'Lower', RANK-1, RANK-1, ZERO, ZERO,
     $                      ABCD(IROW+1,MM1+NU-RANK), LDABCD )
         END IF
C
         RO = RANK
      END IF
C
C     Determine the left Kronecker indices (row indices).
C
      KRONL(IK) = KRONL(IK) + TAU - RO
      NKROL = NKROL + KRONL(IK)
      IK = IK + 1
C
C     C and D are updated to [A21 ; C11] and [B2 ; RD].
C
      NU = NU - RO
      MU = SIGMA + RO
      IF ( RO.NE.0 )
     $   GO TO 20
C
   80 CONTINUE
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of AB08NX ***
      END
