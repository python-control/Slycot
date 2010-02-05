      SUBROUTINE TF01MY( N, M, P, NY, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   U, LDU, X, Y, LDY, DWORK, LDWORK, INFO )
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
C     To compute the output sequence of a linear time-invariant
C     open-loop system given by its discrete-time state-space model
C     (A,B,C,D), where A is an N-by-N general matrix.
C
C     The initial state vector x(1) must be supplied by the user.
C
C     This routine differs from SLICOT Library routine TF01MD in the
C     way the input and output trajectories are stored.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NY      (input) INTEGER
C             The number of output vectors y(k) to be computed.
C             NY >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input matrix B of the system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             output matrix C of the system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             direct link matrix D of the system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     U       (input) DOUBLE PRECISION array, dimension (LDU,M)
C             The leading NY-by-M part of this array must contain the
C             input vector sequence u(k), for k = 1,2,...,NY.
C             Specifically, the k-th row of U must contain u(k)'.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,NY).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the initial state vector
C             x(1) which consists of the N initial states of the system.
C             On exit, this array contains the final state vector
C             x(NY+1) of the N states of the system at instant NY+1.
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY,P)
C             The leading NY-by-P part of this array contains the output
C             vector sequence y(1),y(2),...,y(NY) such that the k-th
C             row of Y contains y(k)' (the outputs at instant k),
C             for k = 1,2,...,NY.
C
C     LDY     INTEGER
C             The leading dimension of array Y.  LDY >= MAX(1,NY).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= N.
C             For better performance, LDWORK should be larger.
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
C     Given an initial state vector x(1), the output vector sequence
C     y(1), y(2),..., y(NY) is obtained via the formulae
C
C        x(k+1) = A x(k) + B u(k)
C        y(k)   = C x(k) + D u(k),
C
C     where each element y(k) is a vector of length P containing the
C     outputs at instant k and k = 1,2,...,NY.
C
C     REFERENCES
C
C     [1] Luenberger, D.G.
C         Introduction to Dynamic Systems: Theory, Models and
C         Applications.
C         John Wiley & Sons, New York, 1979.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately (N + M) x (N + P) x NY
C     multiplications and additions.
C
C     FURTHER COMMENTS
C
C     The implementation exploits data locality and uses BLAS 3
C     operations as much as possible, given the workspace length.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Discrete-time system, multivariable system, state-space model,
C     state-space representation, time response.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDU, LDWORK, LDY, M,
     $                  N, NY, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), U(LDU,*), X(*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           IK, IREM, IS, IYL, MAXN, NB, NS
      DOUBLE PRECISION  UPD
C     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      MAXN = MAX( 1, N )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( NY.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAXN ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAXN ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.MAX( 1, NY ) ) THEN
         INFO = -14
      ELSE IF( LDY.LT.MAX( 1, NY ) ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.N ) THEN
         INFO = -19
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TF01MY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( NY, P ).EQ.0 ) THEN
         RETURN
      ELSE IF ( N.EQ.0 ) THEN
C
C        Non-dynamic system: compute the output vectors.
C
         IF ( M.EQ.0 ) THEN
            CALL DLASET( 'Full', NY, P, ZERO, ZERO, Y, LDY )
         ELSE
            CALL DGEMM( 'No transpose', 'Transpose', NY, P, M, ONE,
     $                  U, LDU, D, LDD, ZERO, Y, LDY )
         END IF
         RETURN
      END IF
C
C     Determine the block size (taken as for LAPACK routine DGETRF).
C
      NB = ILAENV( 1, 'DGETRF', ' ', NY, MAX( M, P ), -1, -1 )
C
C     Find the number of state vectors that can be accommodated in
C     the provided workspace and initialize.
C
      NS = MIN( LDWORK/N, NB*NB/N, NY )
C
      IF ( NS.LE.1 .OR. NY*MAX( M, P ).LE.NB*NB ) THEN
C
C        LDWORK < 2*N or small problem:
C                     only BLAS 2 calculations are used in the loop
C                     for computing the output corresponding to D = 0.
C        One row of the array Y is computed for each loop index value.
C
         DO 10 IK = 1, NY
            CALL DGEMV( 'No transpose', P, N, ONE, C, LDC, X, 1, ZERO,
     $                  Y(IK,1), LDY )
C
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA, X, 1, ZERO,
     $                  DWORK, 1 )
            CALL DGEMV( 'No transpose', N, M, ONE, B, LDB, U(IK,1), LDU,
     $                  ONE, DWORK, 1 )
C
            CALL DCOPY( N, DWORK, 1, X, 1 )
   10    CONTINUE
C
      ELSE
C
C        LDWORK >= 2*N and large problem:
C        some BLAS 3 calculations can also be used.
C
         IYL = ( NY/NS )*NS
         IF ( M.EQ.0 ) THEN
            UPD = ZERO
         ELSE
            UPD = ONE
         END IF
C
         CALL DCOPY( N, X, 1, DWORK, 1 )
C
         DO 30 IK = 1, IYL, NS
C
C           Compute the current NS-1 state vectors in the workspace.
C
            CALL DGEMM( 'No transpose', 'Transpose', N, NS-1, M, ONE,
     $                  B, LDB, U(IK,1), LDU, ZERO, DWORK(N+1), MAXN )
C
            DO 20 IS = 1, NS - 1
               CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                     DWORK((IS-1)*N+1), 1, UPD, DWORK(IS*N+1), 1 )
   20       CONTINUE
C
C           Initialize the current NS output vectors.
C
            CALL DGEMM( 'Transpose', 'Transpose', NS, P, N, ONE, DWORK,
     $                  MAXN, C, LDC, ZERO, Y(IK,1), LDY )
C
C           Prepare the next iteration.
C
            CALL DGEMV( 'No transpose', N, M, ONE, B, LDB,
     $                  U(IK+NS-1,1), LDU, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                  DWORK((NS-1)*N+1), 1, UPD, DWORK, 1 )
   30    CONTINUE
C
         IREM = NY - IYL
C
         IF ( IREM.GT.1 ) THEN
C
C           Compute the last IREM output vectors.
C           First, compute the current IREM-1 state vectors.
C
            IK = IYL + 1
            CALL DGEMM( 'No transpose', 'Transpose', N, IREM-1, M, ONE,
     $                  B, LDB, U(IK,1), LDU, ZERO, DWORK(N+1), MAXN )
C
            DO 40 IS = 1, IREM - 1
               CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                     DWORK((IS-1)*N+1), 1, UPD, DWORK(IS*N+1), 1 )
   40       CONTINUE
C
C           Initialize the last IREM output vectors.
C
            CALL DGEMM( 'Transpose', 'Transpose', IREM, P, N, ONE,
     $                  DWORK, MAXN, C, LDC, ZERO, Y(IK,1), LDY )
C
C           Prepare the final state vector.
C
            CALL DGEMV( 'No transpose', N, M, ONE, B, LDB,
     $                  U(IK+IREM-1,1), LDU, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                  DWORK((IREM-1)*N+1), 1, UPD, DWORK, 1 )
C
         ELSE IF ( IREM.EQ.1 ) THEN
C
C           Compute the last 1 output vectors.
C
            CALL DGEMV( 'No transpose', P, N, ONE, C, LDC, DWORK, 1,
     $                   ZERO, Y(IK,1), LDY )
C
C           Prepare the final state vector.
C
            CALL DCOPY( N, DWORK, 1, DWORK(N+1), 1 )
            CALL DGEMV( 'No transpose', N, M, ONE, B, LDB,
     $                  U(IK,1), LDU, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                  DWORK(N+1), 1, UPD, DWORK, 1 )
         END IF
C
C        Set the final state vector.
C
         CALL DCOPY( N, DWORK, 1, X, 1 )
C
      END IF
C
C     Add the direct contribution of the input to the output vectors.
C
      CALL DGEMM( 'No transpose', 'Transpose', NY, P, M, ONE, U, LDU,
     $            D, LDD, ONE, Y, LDY )
C
      RETURN
C *** Last line of TF01MY ***
      END
