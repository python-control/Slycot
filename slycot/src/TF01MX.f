      SUBROUTINE TF01MX( N, M, P, NY, S, LDS, U, LDU, X, Y, LDY,
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
C     To compute the output sequence of a linear time-invariant
C     open-loop system given by its discrete-time state-space model
C     with an (N+P)-by-(N+M) general system matrix S,
C
C            ( A  B )
C        S = (      ) .
C            ( C  D )
C
C     The initial state vector x(1) must be supplied by the user.
C
C     The input and output trajectories are stored as in the SLICOT
C     Library routine TF01MY.
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
C     S       (input) DOUBLE PRECISION array, dimension (LDS,N+M)
C             The leading (N+P)-by-(N+M) part of this array must contain
C             the system matrix S.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,N+P).
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
C             The length of the array DWORK.
C             LDWORK >= 0,        if MIN(N,P,NY) = 0;  otherwise,
C             LDWORK >= N+P,      if M = 0;
C             LDWORK >= 2*N+M+P,  if M > 0.
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
C        ( x(k+1) )     ( x(k) )
C        (        ) = S (      ) ,
C        (  y(k)  )     ( u(k) )
C
C     where each element y(k) is a vector of length P containing the
C     outputs at instant k, and k = 1,2,...,NY.
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
C     The implementation exploits data locality as much as possible,
C     given the workspace length.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 2002.
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
      INTEGER           INFO, LDS, LDU, LDWORK, LDY, M, N, NY, P
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), S(LDS,*), U(LDU,*), X(*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           I, IC, IU, IW, IY, J, JW, K, N2M, N2P, NB, NF,
     $                  NM, NP, NS
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
      NP = N  + P
      NM = N  + M
      IW = NM + NP
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( NY.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDS.LT.MAX( 1, NP ) ) THEN
         INFO = -6
      ELSE IF( LDU.LT.MAX( 1, NY ) ) THEN
         INFO = -8
      ELSE IF( LDY.LT.MAX( 1, NY ) ) THEN
         INFO = -11
      ELSE
         IF( MIN( N, P, NY ).EQ.0 ) THEN
            JW = 0
         ELSE IF( M.EQ.0 ) THEN
            JW = NP
         ELSE
            JW = IW
         END IF
         IF( LDWORK.LT.JW )
     $      INFO = -13
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TF01MX', -INFO )
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
     $                  U, LDU, S, LDS, ZERO, Y, LDY )
         END IF
         RETURN
      END IF
C
C     Determine the block size (taken as for LAPACK routine DGETRF).
C
      NB = ILAENV( 1, 'DGETRF', ' ', NY, MAX( M, P ), -1, -1 )
C
C     Find the number of state vectors, extended with inputs (if M > 0)
C     and outputs, that can be accommodated in the provided workspace.
C
      NS  = MIN( LDWORK/JW, NB*NB/JW, NY )
      N2P = N + NP
C
      IF ( M.EQ.0 ) THEN
C
C        System with no inputs.
C        Workspace: need   N + P;
C                   prefer larger.
C
         IF( NS.LE.1 .OR. NY*P.LE.NB*NB ) THEN
            IY = N + 1
C
C           LDWORK < 2*(N+P), or small problem.
C           One row of array Y is computed for each loop index value.
C
            DO 10 I = 1, NY
C
C              Compute
C
C              /x(i+1)\    /A\
C              |      | =  | | * x(i).
C              \ y(i) /    \C/
C
               CALL DGEMV( 'NoTranspose', NP, N, ONE, S, LDS, X, 1,
     $                     ZERO, DWORK, 1 )
               CALL DCOPY( N, DWORK, 1, X, 1 )
               CALL DCOPY( P, DWORK(IY), 1, Y(I,1), LDY )
   10       CONTINUE
C
         ELSE
C
C           LDWORK >= 2*(N+P), and large problem.
C           NS rows of array Y are computed before being saved.
C
            NF = ( NY/NS )*NS
            CALL DCOPY( N, X, 1, DWORK, 1 )
C
            DO 40 I = 1, NF, NS
C
C              Compute the current NS extended state vectors in the
C              workspace:
C
C              /x(i+1)\    /A\
C              |      | =  | | * x(i),  i = 1 : ns - 1.
C              \ y(i) /    \C/
C
               DO 20 IC = 1, ( NS - 1 )*NP, NP
                  CALL DGEMV( 'No transpose', NP, N, ONE, S, LDS,
     $                        DWORK(IC), 1, ZERO, DWORK(IC+NP), 1 )
   20          CONTINUE
C
C              Prepare the next iteration.
C
               CALL DGEMV( 'No transpose', NP, N, ONE, S, LDS,
     $                     DWORK((NS-1)*NP+1), 1, ZERO, DWORK, 1 )
C
C              Transpose the NS output vectors in the corresponding part
C              of Y (column-wise).
C
               DO 30 J = 1, P
                  CALL DCOPY( NS-1, DWORK(N2P+J), NP, Y(I,J), 1 )
                  Y(I+NS-1,J) = DWORK(N+J)
   30          CONTINUE
C
   40       CONTINUE
C
            NS = NY - NF
C
            IF ( NS.GT.1 ) THEN
C
C              Compute similarly the last NS output vectors.
C
               DO 50 IC = 1, ( NS - 1 )*NP, NP
                  CALL DGEMV( 'No transpose', NP, N, ONE, S, LDS,
     $                        DWORK(IC), 1, ZERO, DWORK(IC+NP), 1 )
   50          CONTINUE
C
               CALL DGEMV( 'No transpose', NP, N, ONE, S, LDS,
     $                     DWORK((NS-1)*NP+1), 1, ZERO, DWORK, 1 )
C
               DO 60 J = 1, P
                  CALL DCOPY( NS-1, DWORK(N2P+J), NP, Y(NF+1,J), 1 )
                  Y(NF+NS,J) = DWORK(N+J)
   60          CONTINUE
C
            ELSE IF ( NS.EQ.1 ) THEN
C
C              Compute similarly the last NS = 1 output vectors.
C
               CALL DCOPY( N, DWORK, 1, DWORK(NP+1), 1 )
               CALL DGEMV( 'No transpose', NP, N, ONE, S, LDS,
     $                     DWORK(NP+1), 1, ZERO, DWORK, 1 )
               CALL DCOPY( P, DWORK(N+1), 1, Y(NF+1,1), LDY )
C
            END IF
C
C           Set the final state vector.
C
            CALL DCOPY( N, DWORK, 1, X, 1 )
C
         END IF
C
      ELSE
C
C        General case.
C        Workspace: need   2*N + M + P;
C                   prefer larger.
C
         CALL DCOPY( N, X, 1, DWORK, 1 )
C
         IF( NS.LE.1 .OR. NY*( M + P ).LE.NB*NB ) THEN
            IU = N  + 1
            JW = IU + M
            IY = JW + N
C
C           LDWORK < 2*(2*N+M+P), or small problem.
C           One row of array Y is computed for each loop index value.
C
            DO 70 I = 1, NY
C
C              Compute
C
C              /x(i+1)\    /A, B\   /x(i)\
C              |      | =  |    | * |    | .
C              \ y(i) /    \C, D/   \u(i)/
C
               CALL DCOPY( M, U(I,1), LDU, DWORK(IU), 1 )
               CALL DGEMV( 'NoTranspose', NP, NM, ONE, S, LDS, DWORK, 1,
     $                     ZERO, DWORK(JW), 1 )
               CALL DCOPY( N, DWORK(JW), 1, DWORK, 1 )
               CALL DCOPY( P, DWORK(IY), 1, Y(I,1), LDY )
   70       CONTINUE
C
         ELSE
C
C           LDWORK >= 2*(2*N+M+P), and large problem.
C           NS rows of array Y are computed before being saved.
C
            NF  = ( NY/NS )*NS
            N2M = N + NM
C
            DO 110 I = 1, NF, NS
               JW = 1
C
C              Compute the current NS extended state vectors in the
C              workspace:
C
C              /x(i+1)\    /A, B\   /x(i)\
C              |      | =  |    | * |    | ,  i = 1 : ns - 1.
C              \ y(i) /    \C, D/   \u(i)/
C
               DO 80 J = 1, M
                  CALL DCOPY( NS, U(I,J), 1, DWORK(N+J), IW )
   80          CONTINUE
C
               DO 90 K = 1, NS - 1
                  CALL DGEMV( 'No transpose', NP, NM, ONE, S, LDS,
     $                        DWORK(JW), 1, ZERO, DWORK(JW+NM), 1 )
                  JW = JW + NM
                  CALL DCOPY( N, DWORK(JW), 1, DWORK(JW+NP), 1 )
                  JW = JW + NP
   90          CONTINUE
C
C              Prepare the next iteration.
C
               CALL DGEMV( 'No transpose', NP, NM, ONE, S, LDS,
     $                     DWORK(JW), 1, ZERO, DWORK(JW+NM), 1 )
               CALL DCOPY( N, DWORK(JW+NM), 1, DWORK, 1 )
C
C              Transpose the NS output vectors in the corresponding part
C              of Y (column-wise).
C
               DO 100 J = 1, P
                  CALL DCOPY( NS, DWORK(N2M+J), IW, Y(I,J), 1 )
  100          CONTINUE
C
  110       CONTINUE
C
            NS = NY - NF
C
            IF ( NS.GT.1 ) THEN
               JW = 1
C
C              Compute similarly the last NS output vectors.
C
               DO 120 J = 1, M
                  CALL DCOPY( NS, U(NF+1,J), 1, DWORK(N+J), IW )
  120          CONTINUE
C
               DO 130 K = 1, NS - 1
                  CALL DGEMV( 'No transpose', NP, NM, ONE, S, LDS,
     $                        DWORK(JW), 1, ZERO, DWORK(JW+NM), 1 )
                  JW = JW + NM
                  CALL DCOPY( N, DWORK(JW), 1, DWORK(JW+NP), 1 )
                  JW = JW + NP
  130          CONTINUE
C
               CALL DGEMV( 'No transpose', NP, NM, ONE, S, LDS,
     $                     DWORK(JW), 1, ZERO, DWORK(JW+NM), 1 )
               CALL DCOPY( N, DWORK(JW+NM), 1, DWORK, 1 )
C
               DO 140 J = 1, P
                  CALL DCOPY( NS, DWORK(N2M+J), IW, Y(NF+1,J), 1 )
  140          CONTINUE
C
            ELSE IF ( NS.EQ.1 ) THEN
C
C              Compute similarly the last NS = 1 output vectors.
C
               CALL DCOPY( N, DWORK, 1, DWORK(NP+1), 1 )
               CALL DCOPY( M, U(NF+1,1), LDU, DWORK(N2P+1), 1 )
               CALL DGEMV( 'No transpose', NP, NM, ONE, S, LDS,
     $                     DWORK(NP+1), 1, ZERO, DWORK, 1 )
               CALL DCOPY( P, DWORK(N+1), 1, Y(NF+1,1), LDY )
C
            END IF
C
         END IF
C
C        Set the final state vector.
C
         CALL DCOPY( N, DWORK, 1, X, 1 )
C
      END IF
C
      RETURN
C *** Last line of TF01MX ***
      END
