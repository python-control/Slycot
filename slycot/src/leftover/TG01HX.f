      SUBROUTINE TG01HX( COMPQ, COMPZ, L, N, M, P, N1, LBE, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, NR,
     $                   NRBLCK, RTAU, TOL, IWORK, DWORK, INFO )
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
C     Given the descriptor system (A-lambda*E,B,C) with the system
C     matrices A, E and B of the form
C
C            ( A1 X1 )        ( E1 Y1 )        ( B1 )
C        A = (       ) ,  E = (       ) ,  B = (    ) ,
C            ( 0  X2 )        ( 0  Y2 )        ( 0  )
C
C     where
C          - B is an L-by-M matrix, with B1 an N1-by-M  submatrix
C          - A is an L-by-N matrix, with A1 an N1-by-N1 submatrix
C          - E is an L-by-N matrix, with E1 an N1-by-N1 submatrix
C              with LBE nonzero sub-diagonals,
C     this routine reduces the pair (A1-lambda*E1,B1) to the form
C
C     Qc'*[A1-lambda*E1 B1]*diag(Zc,I) =
C
C                              ( Bc Ac-lambda*Ec      *         )
C                              (                                ) ,
C                              ( 0     0         Anc-lambda*Enc )
C
C     where:
C     1) the pencil ( Bc Ac-lambda*Ec ) has full row rank NR for
C        all finite lambda and is in a staircase form with
C                           _      _          _        _
C                         ( A1,0   A1,1  ...  A1,k-1   A1,k  )
C                         (        _          _        _     )
C             ( Bc Ac ) = (  0     A2,1  ...  A2,k-1   A2,k  ) ,  (1)
C                         (              ...  _        _     )
C                         (  0       0   ...  Ak,k-1   Ak,k  )
C
C                           _          _        _
C                         ( E1,1  ...  E1,k-1   E1,k  )
C                         (            _        _     )
C               Ec      = (   0   ...  E2,k-1   E2,k  ) ,         (2)
C                         (       ...           _     )
C                         (   0   ...    0      Ek,k  )
C               _
C         where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank
C                                       _
C         matrix (with rtau(0) = M) and Ei,i is an rtau(i)-by-rtau(i)
C         upper triangular matrix.
C
C      2) the pencil Anc-lambda*Enc is regular of order N1-NR with Enc
C         upper triangular; this pencil contains the uncontrollable
C         finite eigenvalues of the pencil (A1-lambda*E1).
C
C     The transformations are applied to the whole matrices A, E, B
C     and C. The left and/or right orthogonal transformations Qc and Zc
C     performed to reduce the pencil S(lambda) can be optionally
C     accumulated in the matrices Q and Z, respectivelly.
C
C     The reduced order descriptor system (Ac-lambda*Ec,Bc,Cc) has no
C     uncontrollable finite eigenvalues and has the same
C     transfer-function matrix as the original system (A-lambda*E,B,C).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPQ   CHARACTER*1
C             = 'N':  do not compute Q;
C             = 'I':  Q is initialized to the unit matrix, and the
C                     orthogonal matrix Q is returned;
C             = 'U':  Q must contain an orthogonal matrix Q1 on entry,
C                     and the product Q1*Q is returned.
C
C     COMPZ   CHARACTER*1
C             = 'N':  do not compute Z;
C             = 'I':  Z is initialized to the unit matrix, and the
C                     orthogonal matrix Z is returned;
C             = 'U':  Z must contain an orthogonal matrix Z1 on entry,
C                     and the product Z1*Z is returned.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of descriptor state equations; also the number
C             of rows of matrices A, E and B.  L >= 0.
C
C     N       (input) INTEGER
C             The dimension of the descriptor state vector; also the
C             number of columns of matrices A, E and C.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of descriptor system input vector; also the
C             number of columns of matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of descriptor system output; also the
C             number of rows of matrix C.  P >= 0.
C
C     N1      (input) INTEGER
C             The order of subsystem (A1-lambda*E1,B1,C1) to be reduced.
C             MIN(L,N) >= N1 >= 0.
C
C     LBE     (input) INTEGER
C             The number of nonzero sub-diagonals of submatrix E1.
C             MAX(0,N1-1) >= LBE >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the L-by-N state matrix A in the partitioned
C             form
C                      ( A1 X1 )
C                  A = (       ) ,
C                      ( 0  X2 )
C
C             where A1 is N1-by-N1.
C             On exit, the leading L-by-N part of this array contains
C             the transformed state matrix,
C
C                                  ( Ac  *   * )
C                       Qc'*A*Zc = ( 0  Anc  * ) ,
C                                  ( 0   0   * )
C
C             where Ac is NR-by-NR and Anc is (N1-NR)-by-(N1-NR).
C             The matrix ( Bc Ac ) is in the controlability
C             staircase form (1).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the L-by-N descriptor matrix E in the partitioned
C             form
C                      ( E1 Y1 )
C                  E = (       ) ,
C                      ( 0  Y2 )
C
C             where E1 is N1-by-N1 matrix with LBE nonzero
C             sub-diagonals.
C             On exit, the leading L-by-N part of this array contains
C             the transformed descriptor matrix
C
C                                  ( Ec  *   * )
C                       Qc'*E*Zc = ( 0  Enc  * ) ,
C                                  ( 0   0   * )
C
C             where Ec is NR-by-NR and Enc is (N1-NR)-by-(N1-NR).
C             Both Ec and Enc are upper triangular and Enc is
C             nonsingular.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,L).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading L-by-M part of this array must
C             contain the L-by-M input matrix B in the partitioned
C             form
C                      ( B1 )
C                  B = (    ) ,
C                      ( 0  )
C
C             where B1 is N1-by-M.
C             On exit, the leading L-by-M part of this array contains
C             the transformed input matrix
C
C                               ( Bc )
C                       Qc'*B = (    ) ,
C                               ( 0  )
C
C             where Bc is NR-by-M.
C             The matrix ( Bc Ac ) is in the controlability
C             staircase form (1).
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,L).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed matrix C*Zc.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,L)
C             If COMPQ = 'N': Q is not referenced.
C             If COMPQ = 'I': on entry, Q need not be set;
C                             on exit, the leading L-by-L part of this
C                             array contains the orthogonal matrix Q,
C                             where Q' is the product of transformations
C                             which are applied to A, E, and B on
C                             the left.
C             If COMPQ = 'U': on entry, the leading L-by-L part of this
C                             array must contain an orthogonal matrix
C                             Qc;
C                             on exit, the leading L-by-L part of this
C                             array contains the orthogonal matrix
C                             Qc*Q.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= 1,        if COMPQ = 'N';
C             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If COMPZ = 'N': Z is not referenced.
C             If COMPZ = 'I': on entry, Z need not be set;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix Z,
C                             which is the product of transformations
C                             applied to A, E, and C on the right.
C             If COMPZ = 'U': on entry, the leading N-by-N part of this
C                             array must contain an orthogonal matrix
C                             Zc;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix
C                             Zc*Z.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.
C             LDZ >= 1,        if COMPZ = 'N';
C             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'.
C
C     NR      (output) INTEGER
C             The order of the reduced matrices Ac and Ec, and the
C             number of rows of the reduced matrix Bc; also the order of
C             the controllable part of the pair (B, A-lambda*E).
C
C     NRBLCK  (output) INTEGER                      _
C             The number k, of full row rank blocks Ai,i in the
C             staircase form of the pencil (Bc Ac-lambda*Ec) (see (1)
C             and (2)).
C
C     RTAU    (output) INTEGER array, dimension (N1)
C             RTAU(i), for i = 1, ..., NRBLCK, is the row dimension of
C                                     _
C             the full row rank block Ai,i-1 in the staircase form (1).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determinations when
C             transforming (A-lambda*E, B). If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             reciprocal condition numbers in rank determinations; a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = L*N*EPS,  is used instead, where
C             EPS is the machine precision (see LAPACK Library routine
C             DLAMCH).  TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M)
C
C     DWORK   DOUBLE PRECISION array, dimension MAX(N,L,2*M)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The subroutine is based on the reduction algorithm of [1].
C
C     REFERENCES
C
C     [1] A. Varga
C         Computation of Irreducible Generalized State-Space
C         Realizations.
C         Kybernetika, vol. 26, pp. 89-106, 1990.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable and requires
C     0( N*N1**2 )  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the RASP routine RPDS05.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 1999,
C     May 2003, Nov. 2003.
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, Nov. 2003.
C
C     KEYWORDS
C
C     Controllability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            IMAX, IMIN
      PARAMETER          ( IMAX = 1, IMIN = 2 )
      DOUBLE PRECISION   ONE, P05, ZERO
      PARAMETER          ( ONE = 1.0D0, P05 = 0.05D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            INFO, L, LBE, LDA, LDB, LDC, LDE, LDQ, LDZ, M,
     $                   N, N1, NR, NRBLCK, P
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      INTEGER            IWORK( * ), RTAU( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK( * ), E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      LOGICAL            ILQ, ILZ, WITHC
      INTEGER            I, IC, ICOL, ICOMPQ, ICOMPZ, IROW, ISMAX,
     $                   ISMIN, J, K, MN, NF, NR1, RANK, TAUIM1
      DOUBLE PRECISION   CO, C1, C2, RCOND, SMAX, SMAXPR, SMIN, SMINPR,
     $                   SVLMAX, S1, S2, SI, T, TT
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLANGE, DLAPY2, DNRM2, IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL           DLACPY, DLARF, DLARFG, DLARTG, DLASET, DROT,
     $                   DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
C
C     .. Executable Statements ..
C
C     Decode COMPQ.
C
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'U' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
C
C     Decode COMPZ.
C
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'U' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
C
C     Test the input scalar parameters.
C
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( N1.LT.0 .OR. N1.GT.MIN( L, N ) ) THEN
         INFO = -7
      ELSE IF( LBE.LT.0 .OR. LBE.GT.MAX( 0, N1-1 ) ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -10
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, L ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( ( ILQ .AND. LDQ.LT.L ) .OR. LDQ.LT.1 ) THEN
         INFO = -18
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -20
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -24
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01HX', -INFO )
         RETURN
      END IF
C
C     Initialize Q and Z if necessary.
C
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', L, L, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Initialize output variables.
C
      NR = 0
      NRBLCK = 0
C
C     Quick return if possible.
C
      IF( M.EQ.0 .OR. N1.EQ.0 ) THEN
         RETURN
      END IF
C
      WITHC  = P.GT.0
      SVLMAX = DLAPY2( DLANGE( 'F', L, M, B, LDB, DWORK ),
     $                 DLANGE( 'F', L, N, A, LDA, DWORK ) )
      RCOND = TOL
      IF ( RCOND.LE.ZERO ) THEN
C
C        Use the default tolerance in controllability determination.
C
         RCOND = DBLE( L*N )*DLAMCH( 'EPSILON' )
      END IF
C
      IF ( SVLMAX.LT.RCOND )
     $   SVLMAX = ONE
C
C     Reduce E to upper triangular form if necessary.
C
      IF( LBE.GT.0 ) THEN
         DO 10 I = 1, N1-1
C
C           Generate elementary reflector H(i) to annihilate
C           E(i+1:i+lbe,i).
C
            K = MIN( LBE, N1-I ) + 1
            CALL DLARFG( K, E(I,I), E(I+1,I), 1, TT )
            T = E(I,I)
            E(I,I) = ONE
C
C           Apply H(i) to E(i:n1,i+1:n) from the left.
C
            CALL DLARF( 'Left', K, N-I, E(I,I), 1, TT,
     $                  E(I,I+1), LDE, DWORK )
C
C           Apply H(i) to A(i:n1,1:n) from the left.
C
            CALL DLARF( 'Left', K, N, E(I,I), 1, TT,
     $                  A(I,1), LDA, DWORK )
C
C           Apply H(i) to B(i:n1,1:m) from the left.
C
            CALL DLARF( 'Left', K, M, E(I,I), 1, TT,
     $                  B(I,1), LDB, DWORK )
            IF( ILQ ) THEN
C
C              Apply H(i) to Q(1:l,i:n1) from the right.
C
               CALL DLARF( 'Right', L, K, E(I,I), 1, TT,
     $                     Q(1,I), LDQ, DWORK )
            END IF
            E(I,I) = T
   10    CONTINUE
         IF( N1.GT.1 )
     $      CALL DLASET( 'Lower', N1-1, N1-1, ZERO, ZERO, E(2,1), LDE )
      END IF
C
      ISMIN = 1
      ISMAX = ISMIN + M
      IC = -M
      TAUIM1 = M
      NF = N1
C
   20 CONTINUE
      NRBLCK = NRBLCK + 1
      RANK = 0
      IF( NF.GT.0 ) THEN
C
C        IROW will point to the current pivot line in B,
C        ICOL+1 will point to the first active columns of A.
C
         ICOL = IC + TAUIM1
         IROW = NR
         NR1 = NR + 1
         IF( NR.GT.0 )
     $      CALL DLACPY( 'Full', NF, TAUIM1, A(NR1,IC+1), LDA,
     $                   B(NR1,1), LDB )
C
C        Perform QR-decomposition with column pivoting on the current B
C        while keeping E upper triangular.
C        The current B is at first iteration B and for subsequent
C        iterations the NF-by-TAUIM1 matrix delimited by rows
C        NR + 1 to N1 and columns IC + 1 to IC + TAUIM1 of A.
C        The rank of current B is computed in RANK.
C
         IF( TAUIM1.GT.1 ) THEN
C
C           Compute column norms.
C
            DO 30 J = 1, TAUIM1
               DWORK(J) = DNRM2( NF, B(NR1,J), 1 )
               DWORK(M+J) = DWORK(J)
               IWORK(J) = J
   30       CONTINUE
         END IF
C
         MN = MIN( NF, TAUIM1 )
C
   40    CONTINUE
         IF( RANK.LT.MN ) THEN
            J = RANK + 1
            IROW = IROW + 1
C
C           Pivot if necessary.
C
            IF( J.NE.TAUIM1 ) THEN
               K = ( J - 1 ) + IDAMAX( TAUIM1-J+1, DWORK(J), 1 )
               IF( K.NE.J ) THEN
                  CALL DSWAP( NF, B(NR1,J), 1, B(NR1,K), 1 )
                  I = IWORK(K)
                  IWORK(K) = IWORK(J)
                  IWORK(J) = I
                  DWORK(K) = DWORK(J)
                  DWORK(M+K) = DWORK(M+J)
               END IF
            END IF
C
C           Zero elements below the current diagonal element of B.
C
            DO 50 I = N1-1, IROW, -1
C
C              Rotate rows I and I+1 to zero B(I+1,J).
C
               T = B(I,J)
               CALL DLARTG( T, B(I+1,J), CO, SI, B(I,J) )
               B(I+1,J) = ZERO
               CALL DROT( N-I+1, E(I,I), LDE, E(I+1,I), LDE, CO, SI )
               IF( J.LT.TAUIM1 )
     $             CALL DROT( TAUIM1-J, B(I,J+1), LDB,
     $                        B(I+1,J+1), LDB, CO, SI )
               CALL DROT( N-ICOL, A(I,ICOL+1), LDA,
     $                    A(I+1,ICOL+1), LDA, CO, SI )
               IF( ILQ ) CALL DROT( L, Q(1,I), 1, Q(1,I+1), 1, CO, SI )
C
C              Rotate columns I, I+1 to zero E(I+1,I).
C
               T = E(I+1,I+1)
               CALL DLARTG( T, E(I+1,I), CO, SI, E(I+1,I+1) )
               E(I+1,I) = ZERO
               CALL DROT( I, E(1,I+1), 1, E(1,I), 1, CO, SI )
               CALL DROT( N1, A(1,I+1), 1, A(1,I), 1, CO, SI )
               IF( ILZ ) CALL DROT( N, Z(1,I+1), 1, Z(1,I), 1, CO, SI )
               IF( WITHC )
     $            CALL DROT( P, C(1,I+1), 1, C(1,I), 1, CO, SI )
   50       CONTINUE
C
            IF( RANK.EQ.0 ) THEN
C
C              Initialize; exit if matrix is zero (RANK = 0).
C
               SMAX = ABS( B(NR1,1) )
               IF ( SMAX.EQ.ZERO ) GO TO 80
               SMIN = SMAX
               SMAXPR = SMAX
               SMINPR = SMIN
               C1 = ONE
               C2 = ONE
            ELSE
C
C              One step of incremental condition estimation.
C
               CALL DLAIC1( IMIN, RANK, DWORK(ISMIN), SMIN,
     $                      B(NR1,J), B(IROW,J), SMINPR, S1, C1 )
               CALL DLAIC1( IMAX, RANK, DWORK(ISMAX), SMAX,
     $                      B(NR1,J), B(IROW,J), SMAXPR, S2, C2 )
            END IF
C
C           Check the rank; finish the loop if rank loss occurs.
C
            IF( SVLMAX*RCOND.LE.SMAXPR ) THEN
               IF( SVLMAX*RCOND.LE.SMINPR ) THEN
                  IF( SMAXPR*RCOND.LE.SMINPR ) THEN
C
C                    Finish the loop if last row.
C
                     IF( IROW.EQ.N1 ) THEN
                        RANK = RANK + 1
                        GO TO 80
                     END IF
C
C                    Update partial column norms.
C
                     DO 60 I = J + 1, TAUIM1
                        IF( DWORK(I).NE.ZERO ) THEN
                           T = ONE - ( ABS( B(IROW,I) )/DWORK(I) )**2
                           T = MAX( T, ZERO )
                           TT = ONE + P05*T*( DWORK(I)/DWORK(M+I) )**2
                           IF( TT.NE.ONE ) THEN
                              DWORK(I) = DWORK(I)*SQRT( T )
                           ELSE
                              DWORK(I) = DNRM2( NF-J, B(IROW+1,I), 1 )
                              DWORK(M+I) = DWORK(I)
                           END IF
                        END IF
   60                CONTINUE
C
                     DO 70 I = 1, RANK
                        DWORK(ISMIN+I-1) = S1*DWORK(ISMIN+I-1)
                        DWORK(ISMAX+I-1) = S2*DWORK(ISMAX+I-1)
   70                CONTINUE
C
                     DWORK(ISMIN+RANK) = C1
                     DWORK(ISMAX+RANK) = C2
                     SMIN = SMINPR
                     SMAX = SMAXPR
                     RANK = RANK + 1
                     GO TO 40
                  END IF
               END IF
            END IF
            IF( NR.GT.0 ) THEN
               CALL DLASET( 'Full', N1-IROW+1, TAUIM1-J+1, ZERO, ZERO,
     $                      B(IROW,J), LDB )
            END IF
            GO TO 80
         END IF
      END IF
C
   80 IF( RANK.GT.0 ) THEN
         RTAU(NRBLCK) = RANK
C
C        Back permute interchanged columns.
C
         IF( TAUIM1.GT.1 ) THEN
            DO 100 J = 1, TAUIM1
               IF( IWORK(J).GT.0 ) THEN
                  K = IWORK(J)
                  IWORK(J) = -K
   90             CONTINUE
                  IF( K.NE.J ) THEN
                     CALL DSWAP( RANK, B(NR1,J), 1, B(NR1,K), 1 )
                     IWORK(K) = -IWORK(K)
                     K = -IWORK(K)
                     GO TO 90
                  END IF
               END IF
  100       CONTINUE
         END IF
      END IF
      IF( NR.GT.0 )
     $   CALL DLACPY( 'Full', NF, TAUIM1, B(NR1,1), LDB,
     $                A(NR1,IC+1), LDA )
      IF( RANK.GT.0 ) THEN
         NR = NR + RANK
         NF = NF - RANK
         IC = IC + TAUIM1
         TAUIM1 = RANK
         GO TO 20
      ELSE
         NRBLCK = NRBLCK - 1
      END IF
C
      IF( NRBLCK.GT.0 ) RANK = RTAU(1)
      IF( RANK.LT.N1 )
     $   CALL DLASET( 'Full', N1-RANK, M, ZERO, ZERO, B(RANK+1,1), LDB )
C
      RETURN
C *** Last line of TG01HX ***
      END
