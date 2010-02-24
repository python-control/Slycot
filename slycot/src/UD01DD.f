      SUBROUTINE UD01DD( M, N, NIN, A, LDA, INFO )
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
C     To read the elements of a sparse matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     NIN     (input) INTEGER
C             The input channel from which the elements of A are read.
C             NIN >= 0.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading M-by-N part of this array contains the sparse
C             matrix A. The not assigned elements are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1 : if a row index i is read with i < 1 or i > M or
C                   a column index j is read with j < 1 or j > N.
C                   This is a warning.
C
C     METHOD
C
C     First, the elements A(i,j) with 1 <= i <= M and 1 <= j <= N are
C     set to zero. Next the nonzero elements are read from the input
C     file NIN. Each line of NIN must contain consecutively the values
C     i, j, A(i,j). The routine terminates after the last line has been
C     read.
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998.
C     Based on routine RDSPAR by A.J. Geurts, Eindhoven University of
C     Technology, Holland.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N, NIN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           I, J
      DOUBLE PRECISION  AIJ
C     .. External Subroutines ..
      EXTERNAL          DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     .. Executable statements ..
C
      INFO = 0
C
C     Check the input scalar arguments.
C
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NIN.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'UD01DD', -INFO )
         RETURN
      END IF
C
      CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
C
C     Read (i, j, A(i,j)) of the nonzero elements one by one.
C
   10 READ( NIN, FMT = *, END = 20 ) I, J, AIJ
      IF ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         INFO = 1
      ELSE
         A(I,J) = AIJ
      END IF
      GO TO 10
   20 CONTINUE
C
      RETURN
C *** Last line of UD01DD ***
      END
