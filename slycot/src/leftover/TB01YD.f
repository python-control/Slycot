      SUBROUTINE TB01YD( N, M, P, A, LDA, B, LDB, C, LDC, INFO )
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
C     To apply a special similarity transformation to a system given as
C     a triple (A,B,C),
C
C        A <-- P * A * P,  B <-- P * B,  C <-- C * P,
C
C     where P is a matrix with 1 on the secondary diagonal, and with 0
C     in the other entries.
C
C     ARGUMENTS
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
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the system state matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the transformed matrix P*A*P.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the system input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed matrix P*B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,N) if M > 0.
C             LDB >= 1        if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the system output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed matrix C*P.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,P).
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
C     The rows and/or columns of the matrices of the triplet (A,B,C)
C     are swapped in a special way.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004.
C
C     KEYWORDS
C
C     Matrix algebra, matrix operations, similarity transformation.
C
C  *********************************************************************
C
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, M, N, P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
C     ..
C     .. Local Scalars ..
      INTEGER            J, NBY2
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSWAP, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input arguments.
C
      INFO  = 0
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -9
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01YD', -INFO )
         RETURN
      END IF
C
      IF( N.LE.1 )
     $   RETURN
C
C     Transform the matrix A.
C
      NBY2 = N/2
C
      DO 10 J = 1, NBY2
         CALL DSWAP( N, A( 1, J ), -1, A( 1, N-J+1 ), 1 )
   10 CONTINUE
C
      IF( MOD( N, 2 ).NE.0 .AND. N.GT.2 )
     $   CALL DSWAP( NBY2, A( NBY2+2, NBY2+1 ), -1, A( 1, NBY2+1 ), 1 )
C
      IF( M.GT.0 ) THEN
C
C        Transform the matrix B.
C
         DO 20 J = 1, NBY2
            CALL DSWAP( M, B( J, 1 ), LDB, B( N-J+1, 1 ), LDB )
   20    CONTINUE
C
      END IF
C
      IF( P.GT.0 ) THEN
C
C        Transform the matrix C.
C
         DO 30 J = 1, NBY2
            CALL DSWAP( P, C( 1, J ), 1, C( 1, N-J+1 ), 1 )
   30    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of TB01YD ***
      END
