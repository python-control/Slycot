      SUBROUTINE MA02BD( SIDE, M, N, A, LDA )
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
C     To reverse the order of rows and/or columns of a given matrix A
C     by pre-multiplying and/or post-multiplying it, respectively, with
C     a permutation matrix P, where P is a square matrix of appropriate
C     order, with ones down the secondary diagonal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies the operation to be performed, as follows:
C             = 'L': the order of rows of A is to be reversed by
C                    pre-multiplying A with P;
C             = 'R': the order of columns of A is to be reversed by
C                    post-multiplying A with P;
C             = 'B': both the order of rows and the order of columns
C                    of A is to be reversed by pre-multiplying and
C                    post-multiplying A with P.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the given matrix whose rows and/or columns are to
C             be permuted.
C             On exit, the leading M-by-N part of this array contains
C             the matrix P*A if SIDE = 'L', or A*P if SIDE = 'R', or
C             P*A*P if SIDE = 'B'.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine PAP.
C
C     REVISIONS
C
C     -
C
C    ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*)
C     .. Local Scalars ..
      LOGICAL            BSIDES
      INTEGER            I, J, K, M2, N2
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DSWAP
C     .. Executable Statements ..
C
      BSIDES  = LSAME( SIDE, 'B' )
C
      IF( ( LSAME( SIDE, 'L' ) .OR. BSIDES ) .AND. M.GT.1 ) THEN
C
C        Compute P*A.
C
         M2 = M/2
         K = M - M2 + 1
         DO 10 J = 1, N
            CALL DSWAP( M2, A(1,J), -1, A(K,J), 1 )
   10    CONTINUE
      END IF
      IF( ( LSAME( SIDE, 'R' ) .OR. BSIDES ) .AND. N.GT.1 ) THEN
C
C        Compute A*P.
C
         N2 = N/2
         K = N - N2 + 1
         DO 20 I = 1, M
            CALL DSWAP( N2, A(I,1), -LDA, A(I,K), LDA )
   20    CONTINUE
      END IF
C
      RETURN
C *** Last line of MA02BD ***
      END
