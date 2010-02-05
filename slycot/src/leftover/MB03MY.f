      DOUBLE PRECISION FUNCTION MB03MY( NX, X, INCX )
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
C     To compute the absolute minimal value of NX elements in an array.
C     The function returns the value zero if NX < 1.
C
C     ARGUMENTS
C
C     NX      (input) INTEGER
C             The number of elements in X to be examined.
C
C     X       (input) DOUBLE PRECISION array, dimension (NX * INCX)
C             The one-dimensional array of which the absolute minimal
C             value of the elements is to be computed.
C             This array is not referenced if NX < 1.
C
C     INCX    (input) INTEGER
C             The increment to be taken in the array X, defining the
C             distance between two consecutive elements.  INCX >= 1.
C             INCX = 1, if all elements are contiguous in memory.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MB03AZ by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     June 16, 1997.
C
C     KEYWORDS
C
C     None.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INCX, NX
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  DX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF ( NX.LE.0 ) THEN
         MB03MY = ZERO
         RETURN
      END IF
C
      MB03MY = ABS( X(1) )
C
      DO 20 I = 1+INCX, NX*INCX, INCX
         DX = ABS( X(I) )
         IF ( DX.LT.MB03MY ) MB03MY = DX
   20 CONTINUE
C
      RETURN
C *** Last line of MB03MY ***
      END
