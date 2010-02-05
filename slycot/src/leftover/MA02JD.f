      DOUBLE PRECISION FUNCTION MA02JD( LTRAN1, LTRAN2, N, Q1, LDQ1, Q2,
     $                                  LDQ2, RES, LDRES )
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
C     To compute || Q^T Q - I ||_F for a matrix of the form
C
C                       [  op( Q1 )  op( Q2 ) ]
C                  Q =  [                     ],
C                       [ -op( Q2 )  op( Q1 ) ]
C
C     where Q1 and Q2 are N-by-N matrices. This residual can be used to
C     test wether Q is numerically an orthogonal symplectic matrix.
C
C     FUNCTION VALUE
C
C     MA02JD  DOUBLE PRECISION
C             The computed residual.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LTRAN1  LOGICAL
C             Specifies the form of op( Q1 ) as follows:
C             = .FALSE.:  op( Q1 ) = Q1;
C             = .TRUE. :  op( Q1 ) = Q1'.
C
C     LTRAN2  LOGICAL
C             Specifies the form of op( Q2 ) as follows:
C             = .FALSE.:  op( Q2 ) = Q2;
C             = .TRUE. :  op( Q2 ) = Q2'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices Q1 and Q2.  N >= 0.
C
C     Q1      (input) DOUBLE PRECISION array, dimension (LDQ1,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix op( Q1 ).
C
C     LDQ1    INTEGER
C             The leading dimension of the array Q1.  LDQ1 >= MAX(1,N).
C
C     Q2      (input) DOUBLE PRECISION array, dimension (LDQ2,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix op( Q2 ).
C
C     LDQ2    INTEGER
C             The leading dimension of the array Q2.  LDQ2 >= MAX(1,N).
C
C     Workspace
C
C     RES     DOUBLE PRECISION array, dimension (LDRES,N)
C
C     LDRES   INTEGER
C             The leading dimension of the array RES.  LDRES >= MAX(1,N).
C
C     METHOD
C
C     The routine computes the residual by simple elementary operations.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAORS).
C
C     KEYWORDS
C
C     Elementary operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           LTRAN1, LTRAN2
      INTEGER           LDQ1, LDQ2, LDRES, N
C     .. Array Arguments ..
      DOUBLE PRECISION  Q1(LDQ1,*), Q2(LDQ2,*), RES(LDRES,*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  TEMP
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Subroutines ..
      EXTERNAL          DGEMM
C     .. External Functions ..
      DOUBLE PRECISION  DLANGE, DLAPY2
      EXTERNAL          DLANGE, DLAPY2
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C
C     .. Executable Statements ..
C
      IF ( LTRAN1 ) THEN
         CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q1, LDQ1, ZERO, RES, LDRES )
      ELSE
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q1, LDQ1, ZERO, RES, LDRES )
      END IF
      IF ( LTRAN2 ) THEN
         CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q2, LDQ2, ONE, RES, LDRES )
      ELSE
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q2, LDQ2, ONE, RES, LDRES )
      END IF
      DO 10 I = 1, N
         RES(I,I) = RES(I,I) - ONE
   10 CONTINUE
      TEMP = DLANGE( 'Frobenius', N, N, RES, LDRES, DUMMY )
      IF ( LTRAN1 .AND. LTRAN2 ) THEN
         CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q1, LDQ1, ZERO, RES, LDRES )
         CALL DGEMM( 'No Transpose', 'Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q2, LDQ2, -ONE, RES, LDRES )
      ELSE IF ( LTRAN1 ) THEN
         CALL DGEMM( 'Transpose', 'Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q1, LDQ1, ZERO, RES, LDRES )
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q2, LDQ2, -ONE, RES, LDRES )
      ELSE IF ( LTRAN2 ) THEN
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q1, LDQ1, ZERO, RES, LDRES )
         CALL DGEMM( 'Transpose', 'Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q2, LDQ2, -ONE, RES, LDRES )
      ELSE
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE, Q2,
     $               LDQ2, Q1, LDQ1, ZERO, RES, LDRES )
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE, Q1,
     $               LDQ1, Q2, LDQ2, -ONE, RES, LDRES )
      END IF
      TEMP = DLAPY2( TEMP, DLANGE( 'Frobenius', N, N, RES, LDRES,
     $                             DUMMY ) )
      MA02JD = SQRT( TWO )*TEMP
      RETURN
C *** Last line of MA02JD ***
      END
