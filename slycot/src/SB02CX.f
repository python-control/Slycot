      LOGICAL FUNCTION SB02CX( REIG, IEIG )
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
C     To select the purely imaginary eigenvalues in computing the
C     H-infinity norm of a system.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     REIG    (input) DOUBLE PRECISION
C             The real part of the current eigenvalue considered.
C
C     IEIG    (input) DOUBLE PRECISION
C             The imaginary part of the current eigenvalue considered.
C
C     METHOD
C
C     The function value SB02CX is set to .TRUE. for a purely imaginary
C     eigenvalue and to .FALSE., otherwise.
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
C     P. Hr. Petkov, Technical University of Sofia, May, 1999.
C
C     REVISIONS
C
C     P. Hr. Petkov, Technical University of Sofia, Oct. 2000.
C
C     KEYWORDS
C
C     H-infinity norm, robust control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  HUNDRD
      PARAMETER         ( HUNDRD = 100.0D+0 )
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION  IEIG, REIG
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, TOL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     ..
C     .. Executable Statements ..
C
C     Get the machine precision.
C
      EPS = DLAMCH( 'Epsilon' )
C
C     Set the tolerance in the determination of the purely
C     imaginary eigenvalues.
C
      TOL = HUNDRD*EPS
      SB02CX = ABS( REIG ).LT.TOL
C
      RETURN
C *** Last line of SB02CX ***
      END
