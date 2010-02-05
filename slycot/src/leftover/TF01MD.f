      SUBROUTINE TF01MD( N, M, P, NY, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   U, LDU, X, Y, LDY, DWORK, INFO )
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
C     U       (input) DOUBLE PRECISION array, dimension (LDU,NY)
C             The leading M-by-NY part of this array must contain the
C             input vector sequence u(k), for k = 1,2,...,NY.
C             Specifically, the k-th column of U must contain u(k).
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,M).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the initial state vector
C             x(1) which consists of the N initial states of the system.
C             On exit, this array contains the final state vector
C             x(NY+1) of the N states of the system at instant NY.
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY,NY)
C             The leading P-by-NY part of this array contains the output
C             vector sequence y(1),y(2),...,y(NY) such that the k-th
C             column of Y contains y(k) (the outputs at instant k),
C             for k = 1,2,...,NY.
C
C     LDY     INTEGER
C             The leading dimension of array Y.  LDY >= MAX(1,P).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
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
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TF01AD by S. Van Huffel, Katholieke
C     Univ. Leuven, Belgium.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
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
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDU, LDY, M, N, NY, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), U(LDU,*), X(*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           IK
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
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( NY.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.MAX( 1, M ) ) THEN
         INFO = -14
      ELSE IF( LDY.LT.MAX( 1, P ) ) THEN
         INFO = -17
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TF01MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( P, NY ).EQ.0 ) THEN
         RETURN
      ELSE IF ( N.EQ.0 ) THEN
C
C        Non-dynamic system: compute the output vectors.
C
         IF ( M.EQ.0 ) THEN
            CALL DLASET( 'Full', P, NY, ZERO, ZERO, Y, LDY )
         ELSE
            CALL DGEMM( 'No transpose', 'No transpose', P, NY, M, ONE,
     $                  D, LDD, U, LDU, ZERO, Y, LDY )
         END IF
         RETURN
      END IF
C
      DO 10 IK = 1, NY
         CALL DGEMV( 'No transpose', P, N, ONE, C, LDC, X, 1, ZERO,
     $               Y(1,IK), 1 )
C
         CALL DGEMV( 'No transpose', N, N, ONE, A, LDA, X, 1, ZERO,
     $               DWORK, 1 )
         CALL DGEMV( 'No transpose', N, M, ONE, B, LDB, U(1,IK), 1, ONE,
     $               DWORK, 1 )
C
         CALL DCOPY( N, DWORK, 1, X, 1 )
   10 CONTINUE
C
      CALL DGEMM( 'No transpose', 'No transpose', P, NY, M, ONE, D, LDD,
     $            U, LDU, ONE, Y, LDY )
C
      RETURN
C *** Last line of TF01MD ***
      END
