      SUBROUTINE DG01NY( INDI, N, XR, XI )
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
C     For efficiency, no tests of the input scalar parameters are
C     performed.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, EIGHT
      PARAMETER         ( ZERO=0.0D0, HALF=0.5D0, ONE = 1.0D0,
     $                    TWO=2.0D0, EIGHT=8.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         INDI
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  XI(*), XR(*)
C     .. Local Scalars ..
      LOGICAL           LINDI
      INTEGER           I, J, N2
      DOUBLE PRECISION  AI, AR, BI, BR, HELPI, HELPR, PI2, WHELP, WI,
     $                  WR, WSTPI, WSTPR
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, DBLE, SIN
C     .. Executable Statements ..
C
      LINDI = LSAME( INDI, 'D' )
C
C     Initialisation.
C
      PI2 = EIGHT*ATAN( ONE )
      IF ( LINDI ) PI2 = -PI2
C
      WHELP = PI2/DBLE( 2*N )
      WSTPI = SIN( WHELP )
      WHELP = SIN( HALF*WHELP )
      WSTPR = -TWO*WHELP*WHELP
      WI = ZERO
C
      IF ( LINDI ) THEN
         WR = ONE
         XR(N+1) = XR(1)
         XI(N+1) = XI(1)
      ELSE
         WR = -ONE
      END IF
C
C     Recursion.
C
      N2 = N/2 + 1
      DO 10 I = 1, N2
         J = N + 2 - I
         AR = XR(I) + XR(J)
         AI = XI(I) - XI(J)
         BR = XI(I) + XI(J)
         BI = XR(J) - XR(I)
         IF ( LINDI ) THEN
            AR = HALF*AR
            AI = HALF*AI
            BR = HALF*BR
            BI = HALF*BI
         END IF
         HELPR = WR*BR - WI*BI
         HELPI = WR*BI + WI*BR
         XR(I) = AR + HELPR
         XI(I) = AI + HELPI
         XR(J) = AR - HELPR
         XI(J) = HELPI - AI
         WHELP = WR
         WR = WR + WR*WSTPR - WI*WSTPI
         WI = WI + WI*WSTPR + WHELP*WSTPI
   10 CONTINUE
C
      RETURN
C *** Last line of DG01NY ***
      END
