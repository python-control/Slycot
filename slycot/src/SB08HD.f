      SUBROUTINE SB08HD( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, CR,
     $                   LDCR, DR, LDDR, IWORK, DWORK, INFO )
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
C     To construct the state-space representation for the system
C     G = (A,B,C,D) from the factors Q = (AQR,BQR,CQ,DQ) and
C     R = (AQR,BQR,CR,DR) of its right coprime factorization
C                       -1
C              G = Q * R  ,
C
C     where G, Q and R are the corresponding transfer-function matrices.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. Also the number of rows of the
C             matrix B and the number of columns of the matrices C and
C             CR. N represents the order of the systems Q and R.
C             N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector. Also the number of columns
C             of the matrices B, D and DR and the number of rows of the
C             matrices CR and DR.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of output vector. Also the number of rows
C             of the matrices C and D.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix AQR of the systems
C             Q and R.
C             On exit, the leading N-by-N part of this array contains
C             the state dynamics matrix of the system G.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix BQR of the systems Q and R.
C             On exit, the leading N-by-M part of this array contains
C             the input/state matrix of the system G.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix CQ of the system Q.
C             On exit, the leading P-by-N part of this array contains
C             the state/output matrix of the system G.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the input/output matrix DQ of the system Q.
C             On exit, the leading P-by-M part of this array contains
C             the input/output matrix of the system G.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     CR      (input) DOUBLE PRECISION array, dimension (LDCR,N)
C             The leading M-by-N part of this array must contain the
C             state/output matrix CR of the system R.
C
C     LDCR    INTEGER
C             The leading dimension of array CR.  LDCR >= MAX(1,M).
C
C     DR      (input/output) DOUBLE PRECISION array, dimension (LDDR,M)
C             On entry, the leading M-by-M part of this array must
C             contain the input/output matrix DR of the system R.
C             On exit, the leading M-by-M part of this array contains
C             the LU factorization of the matrix DR, as computed by
C             LAPACK Library routine DGETRF.
C
C     LDDR    INTEGER
C             The leading dimension of array DR.  LDDR >= MAX(1,M).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M)
C
C     DWORK   DOUBLE PRECISION array, dimension (MAX(1,4*M))
C             On exit, DWORK(1) contains an estimate of the reciprocal
C             condition number of the matrix DR.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the matrix DR is singular;
C             = 2:  the matrix DR is numerically singular (warning);
C                   the calculations continued.
C
C     METHOD
C
C     The subroutine computes the matrices of the state-space
C     representation G = (A,B,C,D) by using the formulas:
C
C                       -1                   -1
C     A = AQR - BQR * DR  * CR,  B = BQR * DR  ,
C                      -1                   -1
C     C = CQ  - DQ * DR  * CR,   D = DQ * DR  .
C
C     REFERENCES
C
C     [1] Varga A.
C         Coprime factors model reduction method based on
C         square-root balancing-free techniques.
C         System Analysis, Modelling and Simulation,
C         vol. 11, pp. 303-311, 1993.
C
C     CONTRIBUTOR
C
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine RCFI.
C     V. Sima, Research Institute for Informatics, Bucharest, Nov. 1998,
C     full BLAS 3 version.
C
C     REVISIONS
C
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     Mar. 2000, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     Coprime factorization, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDCR, LDD, LDDR, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), CR(LDCR,*),
     $                  D(LDD,*), DR(LDDR,*), DWORK(*)
      INTEGER           IWORK(*)
C     .. Local Scalars
      DOUBLE PRECISION  DRNORM, RCOND
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE
C     .. External Subroutines ..
      EXTERNAL          DGECON, DGEMM, DGETRF, DTRSM, MA02GD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO = 0
C
C     Check the scalar input parameters.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDCR.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDDR.LT.MAX( 1, M ) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 )THEN
C
C        Error return.
C
         CALL XERBLA( 'SB08HD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( M.EQ.0 )THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Factor the matrix  DR.  First, compute the 1-norm.
C
      DRNORM = DLANGE( '1-norm', M, M, DR, LDDR, DWORK )
      CALL DGETRF( M, M, DR, LDDR, IWORK, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = 1
         DWORK(1) = ZERO
         RETURN
      END IF
C                         -1
C     Compute B = BQR * DR  , using the factorization P*DR = L*U.
C
      CALL DTRSM( 'Right', 'Upper', 'NoTranspose', 'NonUnit', N, M, ONE,
     $            DR, LDDR, B, LDB )
      CALL DTRSM( 'Right', 'Lower', 'NoTranspose', 'Unit', N, M, ONE,
     $            DR, LDDR, B, LDB )
      CALL MA02GD( N, B, LDB, 1, M, IWORK, -1 )
C                               -1
C     Compute A = AQR - BQR * DR  * CR.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, M, -ONE, B, LDB,
     $            CR, LDCR, ONE, A, LDA )
C                        -1
C     Compute D = DQ * DR  .
C
      CALL DTRSM( 'Right', 'Upper', 'NoTranspose', 'NonUnit', P, M, ONE,
     $            DR, LDDR, D, LDD )
      CALL DTRSM( 'Right', 'Lower', 'NoTranspose', 'Unit', P, M, ONE,
     $            DR, LDDR, D, LDD )
      CALL MA02GD( P, D, LDD, 1, M, IWORK, -1 )
C                             -1
C     Compute C = CQ - DQ * DR  * CR.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, N, M, -ONE, D, LDD,
     $            CR, LDCR, ONE, C, LDC )
C
C     Estimate the reciprocal condition number of DR.
C     Workspace  4*M.
C
      CALL DGECON( '1-norm', M, DR, LDDR, DRNORM, RCOND, DWORK, IWORK,
     $             INFO )
      IF( RCOND.LE.DLAMCH( 'Epsilon' ) )
     $   INFO = 2
C
      DWORK(1) = RCOND
C
      RETURN
C *** Last line of SB08HD ***
      END
