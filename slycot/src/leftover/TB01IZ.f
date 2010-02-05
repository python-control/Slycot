      SUBROUTINE TB01IZ( JOB, N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                   SCALE, INFO )
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
C     To reduce the 1-norm of a system matrix
C
C             S =  ( A  B )
C                  ( C  0 )
C
C     corresponding to the triple (A,B,C), by balancing. This involves
C     a diagonal similarity transformation inv(D)*A*D applied
C     iteratively to A to make the rows and columns of
C                           -1
C                  diag(D,I)  * S * diag(D,I)
C
C     as close in norm as possible.
C
C     The balancing can be performed optionally on the following
C     particular system matrices
C
C              S = A,    S = ( A  B )    or    S = ( A )
C                                                  ( C )
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates which matrices are involved in balancing, as
C             follows:
C             = 'A':  All matrices are involved in balancing;
C             = 'B':  B and A matrices are involved in balancing;
C             = 'C':  C and A matrices are involved in balancing;
C             = 'N':  B and C matrices are not involved in balancing.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A, the number of rows of matrix B
C             and the number of columns of matrix C.
C             N represents the dimension of the state vector.  N >= 0.
C
C     M       (input) INTEGER.
C             The number of columns of matrix B.
C             M represents the dimension of input vector.  M >= 0.
C
C     P       (input) INTEGER.
C             The number of rows of matrix C.
C             P represents the dimension of output vector.  P >= 0.
C
C     MAXRED  (input/output) DOUBLE PRECISION
C             On entry, the maximum allowed reduction in the 1-norm of
C             S (in an iteration) if zero rows or columns are
C             encountered.
C             If MAXRED > 0.0, MAXRED must be larger than one (to enable
C             the norm reduction).
C             If MAXRED <= 0.0, then the value 10.0 for MAXRED is
C             used.
C             On exit, if the 1-norm of the given matrix S is non-zero,
C             the ratio between the 1-norm of the given matrix and the
C             1-norm of the balanced matrix.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the system state matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the balanced matrix inv(D)*A*D.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input/output) COMPLEX*16 array, dimension (LDB,M)
C             On entry, if M > 0, the leading N-by-M part of this array
C             must contain the system input matrix B.
C             On exit, if M > 0, the leading N-by-M part of this array
C             contains the balanced matrix inv(D)*B.
C             The array B is not referenced if M = 0.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,N) if M > 0.
C             LDB >= 1        if M = 0.
C
C     C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C             On entry, if P > 0, the leading P-by-N part of this array
C             must contain the system output matrix C.
C             On exit, if P > 0, the leading P-by-N part of this array
C             contains the balanced matrix C*D.
C             The array C is not referenced if P = 0.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,P).
C
C     SCALE   (output) DOUBLE PRECISION array, dimension (N)
C             The scaling factors applied to S.  If D(j) is the scaling
C             factor applied to row and column j, then SCALE(j) = D(j),
C             for j = 1,...,N.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit.
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Balancing consists of applying a diagonal similarity
C     transformation
C                           -1
C                  diag(D,I)  * S * diag(D,I)
C
C     to make the 1-norms of each row of the first N rows of S and its
C     corresponding column nearly equal.
C
C     Information about the diagonal matrix D is returned in the vector
C     SCALE.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Balancing, eigenvalue, matrix algebra, matrix operations,
C     similarity transformation.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 1.0D+1 )
      DOUBLE PRECISION   FACTOR, MAXR
      PARAMETER          ( FACTOR = 0.95D+0, MAXR = 10.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, LDA, LDB, LDC, M, N, P
      DOUBLE PRECISION   MAXRED
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
      DOUBLE PRECISION   SCALE( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOCONV, WITHB, WITHC
      INTEGER            I, ICA, IRA, J
      DOUBLE PRECISION   CA, CO, F, G, MAXNRM, RA, RO, S, SFMAX1,
     $                   SFMAX2, SFMIN1, SFMIN2, SNORM, SRED
      COMPLEX*16         CDUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           DLAMCH, DZASUM, IZAMAX, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
C     ..
C     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input arguments.
C
      INFO  = 0
      WITHB = LSAME( JOB, 'A' ) .OR. LSAME( JOB, 'B' )
      WITHC = LSAME( JOB, 'A' ) .OR. LSAME( JOB, 'C' )
C
      IF( .NOT.WITHB .AND. .NOT.WITHC .AND. .NOT.LSAME( JOB, 'N' ) )
     $   THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( MAXRED.GT.ZERO .AND. MAXRED.LT.ONE ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( ( M.GT.0 .AND. LDB.LT.MAX( 1, N ) ) .OR.
     $         ( M.EQ.0 .AND. LDB.LT.1 ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TB01IZ', -INFO )
         RETURN
      END IF
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Compute the 1-norm of the required part of matrix S and exit if
C     it is zero.
C
      SNORM = ZERO
C
      DO 10 J = 1, N
         SCALE( J ) = ONE
         CO = DZASUM( N, A( 1, J ), 1 )
         IF( WITHC .AND. P.GT.0 )
     $      CO = CO + DZASUM( P, C( 1, J ), 1 )
         SNORM = MAX( SNORM, CO )
   10 CONTINUE
C
      IF( WITHB ) THEN
C
         DO 20 J = 1, M
            SNORM = MAX( SNORM, DZASUM( N, B( 1, J ), 1 ) )
   20    CONTINUE
C
      END IF
C
      IF( SNORM.EQ.ZERO )
     $   RETURN
C
C     Set some machine parameters and the maximum reduction in the
C     1-norm of S if zero rows or columns are encountered.
C
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
C
      SRED = MAXRED
      IF( SRED.LE.ZERO ) SRED = MAXR
C
      MAXNRM = MAX( SNORM/SRED, SFMIN1 )
C
C     Balance the matrix.
C
C     Iterative loop for norm reduction.
C
   30 CONTINUE
      NOCONV = .FALSE.
C
      DO 90 I = 1, N
         CO = ZERO
         RO = ZERO
C
         DO 40 J = 1, N
            IF( J.EQ.I )
     $         GO TO 40
            CO = CO + CABS1( A( J, I ) )
            RO = RO + CABS1( A( I, J ) )
   40    CONTINUE
C
         ICA = IZAMAX( N, A( 1, I ), 1 )
         CA  = ABS( A( ICA, I ) )
         IRA = IZAMAX( N, A( I, 1 ), LDA )
         RA  = ABS( A( I, IRA ) )
C
         IF( WITHC .AND. P.GT.0 ) THEN
            CO  = CO + DZASUM( P, C( 1, I ), 1 )
            ICA = IZAMAX( P, C( 1, I ), 1 )
            CA  = MAX( CA, ABS( C( ICA, I ) ) )
         END IF
C
         IF( WITHB .AND. M.GT.0 ) THEN
            RO  = RO + DZASUM( M, B( I, 1 ), LDB )
            IRA = IZAMAX( M, B( I, 1 ), LDB )
            RA  = MAX( RA, ABS( B( I, IRA ) ) )
         END IF
C
C        Special case of zero CO and/or RO.
C
         IF( CO.EQ.ZERO .AND. RO.EQ.ZERO )
     $      GO TO 90
         IF( CO.EQ.ZERO ) THEN
            IF( RO.LE.MAXNRM )
     $         GO TO 90
            CO = MAXNRM
         END IF
         IF( RO.EQ.ZERO ) THEN
            IF( CO.LE.MAXNRM )
     $         GO TO 90
            RO = MAXNRM
         END IF
C
C        Guard against zero CO or RO due to underflow.
C
         G = RO / SCLFAC
         F = ONE
         S = CO + RO
   50    CONTINUE
         IF( CO.GE.G .OR. MAX( F, CO, CA ).GE.SFMAX2 .OR.
     $       MIN( RO, G, RA ).LE.SFMIN2 )GO TO 60
         F  =  F*SCLFAC
         CO = CO*SCLFAC
         CA = CA*SCLFAC
         G  =  G / SCLFAC
         RO = RO / SCLFAC
         RA = RA / SCLFAC
         GO TO 50
C
   60    CONTINUE
         G = CO / SCLFAC
   70    CONTINUE
         IF( G.LT.RO .OR. MAX( RO, RA ).GE.SFMAX2 .OR.
     $       MIN( F, CO, G, CA ).LE.SFMIN2 )GO TO 80
         F  =  F / SCLFAC
         CO = CO / SCLFAC
         CA = CA / SCLFAC
         G  =  G / SCLFAC
         RO = RO*SCLFAC
         RA = RA*SCLFAC
         GO TO 70
C
C        Now balance.
C
   80    CONTINUE
         IF( ( CO+RO ).GE.FACTOR*S )
     $      GO TO 90
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 )
     $         GO TO 90
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F )
     $         GO TO 90
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
C
         CALL ZDSCAL( N, G, A( I, 1 ), LDA )
         CALL ZDSCAL( N, F, A( 1, I ), 1 )
         IF( M.GT.0 ) CALL ZDSCAL( M, G, B( I, 1 ), LDB )
         IF( P.GT.0 ) CALL ZDSCAL( P, F, C( 1, I ), 1 )
C
   90 CONTINUE
C
      IF( NOCONV )
     $   GO TO 30
C
C     Set the norm reduction parameter.
C
      MAXRED = SNORM
      SNORM  = ZERO
C
      DO 100 J = 1, N
         CO = DZASUM( N, A( 1, J ), 1 )
         IF( WITHC .AND. P.GT.0 )
     $      CO = CO + DZASUM( P, C( 1, J ), 1 )
         SNORM = MAX( SNORM, CO )
  100 CONTINUE
C
      IF( WITHB ) THEN
C
         DO 110 J = 1, M
            SNORM = MAX( SNORM, DZASUM( N, B( 1, J ), 1 ) )
  110    CONTINUE
C
      END IF
      MAXRED = MAXRED/SNORM
      RETURN
C *** Last line of TB01IZ ***
      END
