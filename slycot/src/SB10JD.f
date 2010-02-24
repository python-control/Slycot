      SUBROUTINE SB10JD( N, M, NP, A, LDA, B, LDB, C, LDC, D, LDD, E,
     $                   LDE, NSYS, DWORK, LDWORK, INFO )
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
C     To convert the descriptor state-space system
C
C     E*dx/dt = A*x + B*u
C           y = C*x + D*u
C
C     into regular state-space form
C
C      dx/dt = Ad*x + Bd*u
C          y = Cd*x + Dd*u .
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the descriptor system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     NP      (input) INTEGER
C             The row size of the matrix C.  NP >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state matrix A of the descriptor system.
C             On exit, the leading NSYS-by-NSYS part of this array
C             contains the state matrix Ad of the converted system.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the descriptor system.
C             On exit, the leading NSYS-by-M part of this array
C             contains the input matrix Bd of the converted system.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading NP-by-N part of this array must
C             contain the output matrix C of the descriptor system.
C             On exit, the leading NP-by-NSYS part of this array
C             contains the output matrix Cd of the converted system.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading NP-by-M part of this array must
C             contain the matrix D of the descriptor system.
C             On exit, the leading NP-by-M part of this array contains
C             the matrix Dd of the converted system.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix E of the descriptor system.
C             On exit, this array contains no useful information.
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= max(1,N).
C
C     NSYS    (output) INTEGER
C             The order of the converted state-space system.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= max( 1, 2*N*N + 2*N + N*MAX( 5, N + M + NP ) ).
C             For good performance, LDWORK must generally be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the iteration for computing singular value
C                   decomposition did not converge.
C
C     METHOD
C
C     The routine performs the transformations described in [1].
C
C     REFERENCES
C
C     [1] Chiang, R.Y. and Safonov, M.G.
C         Robust Control Toolbox User's Guide.
C         The MathWorks Inc., Natick, Mass., 1992.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000,
C     Feb. 2001.
C
C     KEYWORDS
C
C     Descriptor systems, state-space models.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDE, LDWORK, M, N,
     $                   NP, NSYS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), DWORK( * ),  E( LDE, * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IA12, IA21, IB2, IC2, INFO2, IS, ISA, IU,
     $                   IV, IWRK, J, K, LWA, LWAMAX, MINWRK, NS1
      DOUBLE PRECISION   EPS, SCALE, TOL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DGESVD, DLACPY, DLASET, DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( NP.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, NP ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -11
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -13
      END IF
C
C     Compute workspace.
C
      MINWRK = MAX( 1, 2*N*( N + 1 ) + N*MAX( 5, N + M + NP ) )
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10JD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         NSYS = 0
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Set tol.
C
      EPS = DLAMCH( 'Epsilon' )
      TOL = SQRT( EPS )
C
C     Workspace usage.
C
      IS = 0
      IU = IS + N
      IV = IU + N*N
C
      IWRK = IV + N*N
C
C     Compute the SVD of E.
C     Additional workspace:  need   5*N; prefer larger.
C
      CALL DGESVD( 'S', 'S', N, N, E, LDE, DWORK( IS+1 ), DWORK( IU+1 ),
     $             N, DWORK( IV+1 ), N, DWORK( IWRK+1 ), LDWORK-IWRK,
     $             INFO2 )
      IF( INFO2.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      LWAMAX = MAX( MINWRK, INT( DWORK( IWRK+1 ) + IWRK ) )
C
C     Determine the rank of E.
C
      NS1 = 0
      DO 10 I = 1, N
         IF( DWORK( IS+I ).GT.TOL ) NS1 = NS1 + 1
   10 CONTINUE
      IF( NS1.GT.0 ) THEN
C
C        Transform A.
C        Additional workspace:  need   N*max(N,M,NP).
C
         CALL DGEMM( 'T', 'N', N, N, N, ONE, DWORK( IU+1 ), N, A, LDA,
     $               ZERO, DWORK( IWRK+1 ), N )
         CALL DGEMM( 'N', 'T', N, N, N, ONE, DWORK( IWRK+1 ), N,
     $               DWORK( IV+1 ), N, ZERO, A, LDA )
C
C        Transform B.
C
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK( IWRK+1 ), N )
         CALL DGEMM( 'T', 'N', N, M, N, ONE, DWORK( IU+1 ), N,
     $               DWORK( IWRK+1 ), N, ZERO, B, LDB )
C
C        Transform C.
C
         CALL DLACPY( 'Full', NP, N, C, LDC, DWORK( IWRK+1 ), NP )
         CALL DGEMM( 'N', 'T', NP, N, N, ONE, DWORK( IWRK+1 ), NP,
     $               DWORK( IV+1 ), N, ZERO, C, LDC )
C
         K = N - NS1
         IF( K.GT.0 ) THEN
            ISA  = IU  + K*K
            IV   = ISA + K
            IWRK = IV  + K*MAX( K, NS1 )
C
C           Compute the SVD of A22.
C           Additional workspace:  need   5*K; prefer larger.
C
            CALL DGESVD( 'S', 'S', K, K, A( NS1+1, NS1+1 ), LDA,
     $                   DWORK( ISA+1 ), DWORK( IU+1 ), K,
     $                   DWORK( IV+1 ), K, DWORK( IWRK+1 ), LDWORK-IWRK,
     $                   INFO2 )
            IF( INFO2.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            IA12 = IWRK
            IB2  = IA12 + NS1*K
            IC2  = IB2  + K*M
C
            LWA = INT( DWORK( IWRK+1 ) ) + IWRK
            LWAMAX = MAX( LWA, LWAMAX, IC2 + K*NP )
C
C           Compute the transformed A12.
C
            CALL DGEMM( 'N', 'T', NS1, K, K, ONE, A( 1, NS1+1 ), LDA,
     $                  DWORK( IV+1 ), K, ZERO, DWORK( IA12+1 ), NS1 )
C
C           Compute CC2.
C
            CALL DGEMM( 'N', 'T', NP, K, K, ONE, C( 1, NS1+1 ), LDC,
     $                  DWORK( IV+1 ), K, ZERO, DWORK( IC2+1 ), NP )
C
C           Compute the transformed A21.
C
            IA21 = IV
            CALL DGEMM( 'T', 'N', K, NS1, K, ONE, DWORK( IU+1 ), K,
     $                  A( NS1+1, 1 ), LDA, ZERO, DWORK( IA21+1 ), K )
C
C           Compute BB2.
C
            CALL DGEMM( 'T', 'N', K, M, K, ONE, DWORK( IU+1 ), K,
     $                  B( NS1+1, 1 ), LDB, ZERO, DWORK( IB2+1 ), K )
C
C           Compute A12*pinv(A22) and CC2*pinv(A22).
C
            DO 20 J = 1, K
               SCALE = ZERO
               IF( DWORK( ISA+J ).GT.TOL ) SCALE = ONE/DWORK( ISA+J )
               CALL DSCAL( NS1, SCALE, DWORK( IA12+(J-1)*NS1+1 ), 1 )
               CALL DSCAL( NP,  SCALE, DWORK( IC2+(J-1)*NP+1 ), 1 )
  20        CONTINUE
C
C           Compute Ad.
C
            CALL DGEMM( 'N', 'N', NS1, NS1, K, -ONE, DWORK( IA12+1 ),
     $                  NS1, DWORK( IA21+1 ), K, ONE, A, LDA )
C
C           Compute Bd.
C
            CALL DGEMM( 'N', 'N', NS1, M, K, -ONE, DWORK( IA12+1 ), NS1,
     $                  DWORK( IB2+1 ), K, ONE, B, LDB )
C
C           Compute Cd.
C
            CALL DGEMM( 'N', 'N', NP, NS1, K, -ONE, DWORK( IC2+1 ), NP,
     $                  DWORK( IA21+1 ), K, ONE, C, LDC )
C
C           Compute Dd.
C
            CALL DGEMM( 'N', 'N', NP, M, K, -ONE, DWORK( IC2+1 ), NP,
     $                  DWORK( IB2+1 ), K, ONE, D, LDD )
         END IF
         DO 30 I = 1, NS1
            SCALE = ONE/SQRT( DWORK( IS+I ) )
            CALL DSCAL( NS1, SCALE, A( I, 1 ), LDA )
            CALL DSCAL( M,   SCALE, B( I, 1 ), LDB )
  30     CONTINUE
         DO 40 J = 1, NS1
            SCALE = ONE/SQRT( DWORK( IS+J ) )
            CALL DSCAL( NS1, SCALE, A( 1, J ), 1 )
            CALL DSCAL( NP,  SCALE, C( 1, J ), 1 )
  40     CONTINUE
         NSYS = NS1
      ELSE
         CALL DLASET( 'F', N,  N, ZERO, -ONE/EPS, A, LDA )
         CALL DLASET( 'F', N,  M, ZERO, ZERO, B, LDB )
         CALL DLASET( 'F', NP, N, ZERO, ZERO, C, LDC )
         NSYS = N
      END IF
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10JD ***
      END
