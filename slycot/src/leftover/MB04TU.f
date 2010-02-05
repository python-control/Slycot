      SUBROUTINE MB04TU( N, X, INCX, Y, INCY, C, S )
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
C     To perform the Givens transformation, defined by C (cos) and S
C     (sin), and interchange the vectors involved, i.e.
C
C        |X(i)|    | 0   1 |   | C   S |   |X(i)|
C        |    | := |       | x |       | x |    |, i = 1,...N.
C        |Y(i)|    | 1   0 |   |-S   C |   |Y(i)|
C
C     REMARK. This routine is a modification of DROT from BLAS.
C             This routine is called only by the SLICOT routines MB04TX
C             and MB04VX.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB04FU by Th.G.J. Beelen,
C     Philips Glass Eindhoven, Holland.
C
C     REVISIONS
C
C     January 26, 1998.
C
C     KEYWORDS
C
C     Othogonal transformation.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INCX, INCY, N
      DOUBLE PRECISION  C, S
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DTEMP
      INTEGER           I, IX, IY
C     .. Executable Statements ..
C
      IF ( N.LE.0 ) RETURN
      IF ( ( INCX.NE.1 ) .OR. ( INCY.NE.1 ) ) THEN
C
C        Code for unequal increments or equal increments not equal to 1.
C
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = (-N+1)*INCX + 1
         IF ( INCY.LT.0 ) IY = (-N+1)*INCY + 1
C
         DO 20 I = 1, N
            DTEMP = C*Y(IY) - S*X(IX)
            Y(IY) = C*X(IX) + S*Y(IY)
            X(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
   20    CONTINUE
C
      ELSE
C
C        Code for both increments equal to 1.
C
         DO 40 I = 1, N
            DTEMP = C*Y(I) - S*X(I)
            Y(I) = C*X(I) + S*Y(I)
            X(I) = DTEMP
   40    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MB04TU ***
      END
