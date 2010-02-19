      LOGICAL FUNCTION SB02OW( ALPHAR, ALPHAI, BETA )
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
C     To select the stable generalized eigenvalues for solving the
C     continuous-time algebraic Riccati equation.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     ALPHAR  (input) DOUBLE PRECISION
C             The real part of the numerator of the current eigenvalue
C             considered.
C
C     ALPHAI  (input) DOUBLE PRECISION
C             The imaginary part of the numerator of the current
C             eigenvalue considered.
C
C     BETA    (input) DOUBLE PRECISION
C             The (real) denominator of the current eigenvalue
C             considered. It is assumed that BETA <> 0 (regular case).
C
C     METHOD
C
C     The function value SB02OW is set to .TRUE. for a stable eigenvalue
C     and to .FALSE., otherwise.
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
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997.
C     Supersedes Release 2.0 routine SB02CW by P. Van Dooren, Philips
C     Research Laboratory, Brussels, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, optimal regulator, Schur form.
C
C     ******************************************************************
C
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHAR, ALPHAI, BETA
C     .. Executable Statements ..
C
      SB02OW = ( ALPHAR.LT.ZERO .AND. BETA.GT.ZERO ) .OR.
     $         ( ALPHAR.GT.ZERO .AND. BETA.LT.ZERO )
C
      RETURN
C *** Last line of SB02OW ***
      END
