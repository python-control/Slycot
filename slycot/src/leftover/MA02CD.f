      SUBROUTINE MA02CD( N, KL, KU, A, LDA )
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
C     To compute the pertranspose of a central band of a square matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the square matrix A.  N >= 0.
C
C     KL      (input) INTEGER
C             The number of subdiagonals of A to be pertransposed.
C             0 <= KL <= N-1.
C
C     KU      (input) INTEGER
C             The number of superdiagonals of A to be pertransposed.
C             0 <= KU <= N-1.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain a square matrix whose central band formed from
C             the KL subdiagonals, the main diagonal and the KU
C             superdiagonals will be pertransposed.
C             On exit, the leading N-by-N part of this array contains
C             the matrix A with its central band (the KL subdiagonals,
C             the main diagonal and the KU superdiagonals) pertransposed
C             (that is the elements of each antidiagonal appear in
C             reversed order). This is equivalent to forming P*B'*P,
C             where B is the matrix formed from the central band of A
C             and P is a permutation matrix with ones down the secondary
C             diagonal.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine DMPTR.
C
C     REVISIONS
C
C     A. Varga, December 2001.
C     V. Sima, March 2004.
C
C    ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER          KL, KU, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)
C     .. Local Scalars ..
      INTEGER          I, I1, LDA1
C     .. External Subroutines ..
      EXTERNAL         DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC        MIN
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF( N.LE.1 )
     $   RETURN
C
      LDA1 = LDA + 1
C
C     Pertranspose the KL subdiagonals.
C
      DO 10 I = 1, MIN( KL, N-2 )
         I1 = (N-I) / 2
         IF( I1.GT.0 )
     $      CALL DSWAP( I1, A(I+1,1), LDA1, A(N-I1+1,N-I1+1-I), -LDA1 )
   10 CONTINUE
C
C     Pertranspose the KU superdiagonals.
C
      DO 20 I = 1, MIN( KU, N-2 )
         I1 = (N-I) / 2
         IF( I1.GT.0 )
     $      CALL DSWAP( I1, A(1,I+1), LDA1, A(N-I1+1-I,N-I1+1), -LDA1 )
   20 CONTINUE
C
C     Pertranspose the diagonal.
C
      I1 = N / 2
      IF( I1.GT.0 )
     $   CALL DSWAP( I1, A(1,1), LDA1, A(N-I1+1,N-I1+1), -LDA1 )
C
      RETURN
C *** Last line of MA02CD ***
      END
