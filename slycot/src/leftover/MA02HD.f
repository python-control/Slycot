      LOGICAL FUNCTION MA02HD( JOB, M, N, DIAG, A, LDA )
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
C     To check if A = DIAG*I, where I is an M-by-N matrix with ones on
C     the diagonal and zeros elsewhere.
C
C     FUNCTION VALUE
C
C     MA02HD  LOGICAL
C             The function value is set to .TRUE. if A = DIAG*I, and to
C             .FALSE., otherwise.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the part of the matrix A to be checked out,
C             as follows:
C             = 'U': Upper triangular/trapezoidal part;
C             = 'L': Lower triangular/trapezoidal part.
C             Otherwise:  All of the matrix A.
C
C     Input/Output Parameters
C
C     M      (input) INTEGER
C            The number of rows of the matrix A.  M >= 0.
C
C     N      (input) INTEGER
C            The number of columns of the matrix A.  N >= 0.
C
C     DIAG   (input) DOUBLE PRECISION
C            The scalar DIAG.
C
C     A      (input) DOUBLE PRECISION array, dimension (LDA,N)
C            The leading M-by-N part of this array must contain the
C            matrix A.  If JOB = 'U', only the upper triangle or
C            trapezoid is accessed; if JOB = 'L', only the lower
C            triangle or trapezoid is accessed.
C
C     LDA    INTEGER
C            The leading dimension of the array A.  LDA >= max(1,M).
C
C     METHOD
C
C     The routine returns immediately after detecting a diagonal element
C     which differs from DIAG, or a nonzero off-diagonal element in the
C     searched part of A.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2001.
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2003.
C
C     KEYWORDS
C
C     Elementary operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            LDA, M, N
      DOUBLE PRECISION   DIAG
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*)
C     .. Local Scalars ..
      INTEGER            I, J
C     .. External Functions
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C
C     .. Executable Statements ..
C
C     Do not check parameters, for efficiency.
C
      IF( LSAME( JOB, 'U' ) ) THEN
C
         DO 20 J = 1, N
C
            DO 10 I = 1, MIN( J-1, M )
               IF( A(I,J).NE.ZERO ) THEN
                  MA02HD = .FALSE.
                  RETURN
               END IF
   10       CONTINUE
C
            IF( J.LE.M ) THEN
               IF( A(J,J).NE.DIAG ) THEN
                  MA02HD = .FALSE.
                  RETURN
               END IF
            END IF
   20    CONTINUE
C
      ELSE IF( LSAME( JOB, 'L' ) ) THEN
C
         DO 40 J = 1, MIN( M, N )
            IF( A(J,J).NE.DIAG ) THEN
               MA02HD = .FALSE.
               RETURN
            END IF
C
            IF ( J.NE.M ) THEN
C
               DO 30 I = MIN( J+1, M ), M
                  IF( A(I,J).NE.ZERO ) THEN
                     MA02HD = .FALSE.
                     RETURN
                  END IF
   30          CONTINUE
C
            END IF
   40    CONTINUE
C
      ELSE
C
         DO 70 J = 1, N
C
            DO 50 I = 1, MIN( J-1, M )
               IF( A(I,J).NE.ZERO ) THEN
                  MA02HD = .FALSE.
                  RETURN
               END IF
   50       CONTINUE
C
            IF( J.LE.M ) THEN
               IF( A(J,J).NE.DIAG ) THEN
                  MA02HD = .FALSE.
                  RETURN
               END IF
            END IF
C
            IF ( J.LT.M ) THEN
C
               DO 60 I = MIN( J+1, M ), M
                  IF( A(I,J).NE.ZERO ) THEN
                     MA02HD = .FALSE.
                     RETURN
                  END IF
   60          CONTINUE
C
            END IF
   70    CONTINUE
C
      END IF
C
      MA02HD = .TRUE.
C
      RETURN
C *** Last line of MA02HD ***
      END
