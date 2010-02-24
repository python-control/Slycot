      INTEGER FUNCTION MC01SX( LB, UB, E, MANT )
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
C     To compute the variation V of the exponents of a series of
C     non-zero floating-point numbers: a(j) = MANT(j) * beta**(E(j)),
C     where beta is the base of the machine representation of
C     floating-point numbers, i.e.,
C     V = max(E(j)) - min(E(j)), j = LB,...,UB and MANT(j) non-zero.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC01GX by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION        ZERO
      PARAMETER               ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER                 LB, UB
C     .. Array Arguments ..
      INTEGER                 E(*)
      DOUBLE PRECISION        MANT(*)
C     .. Local Scalars ..
      INTEGER                 J, MAXE, MINE
C     .. Intrinsic Functions ..
      INTRINSIC               MAX, MIN
C     .. Executable Statements ..
C
      MAXE = E(LB)
      MINE = MAXE
C
      DO 20 J = LB + 1, UB
         IF ( MANT(J).NE.ZERO ) THEN
            MAXE = MAX( MAXE, E(J) )
            MINE = MIN( MINE, E(J) )
         END IF
   20 CONTINUE
C
      MC01SX = MAXE - MINE
C
      RETURN
C *** Last line of MC01SX ***
      END
