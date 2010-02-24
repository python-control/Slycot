      SUBROUTINE SB01BX( REIG, N, XR, XI, WR, WI, S, P )
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
C     To choose a real eigenvalue or a pair of complex conjugate
C     eigenvalues at "minimal" distance to a given real or complex
C     value.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     REIG    LOGICAL
C             Specifies the type of eigenvalues as follows:
C             = .TRUE.,  a real eigenvalue is to be selected;
C             = .FALSE., a pair of complex eigenvalues is to be
C                        selected.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of eigenvalues contained in the arrays WR
C             and WI.  N >= 1.
C
C     XR,XI   (input) DOUBLE PRECISION
C             If REIG = .TRUE., XR must contain the real value and XI
C             is assumed zero and therefore not referenced.
C             If REIG = .FALSE., XR must contain the real part and XI
C             the imaginary part, respectively, of the complex value.
C
C     WR,WI   (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, if REIG = .TRUE., WR must contain the real
C             eigenvalues from which an eigenvalue at minimal distance
C             to XR is to be selected. In this case, WI is considered
C             zero and therefore not referenced.
C             On entry, if REIG = .FALSE., WR and WI must contain the
C             real and imaginary parts, respectively, of the eigenvalues
C             from which a pair of complex conjugate eigenvalues at
C             minimal "distance" to XR + jXI is to be selected.
C             The eigenvalues of each pair of complex conjugate
C             eigenvalues must appear consecutively.
C             On exit, the elements of these arrays are reordered such
C             that the selected eigenvalue(s) is (are) found in the
C             last element(s) of these arrays.
C
C     S,P     (output) DOUBLE PRECISION
C             If REIG = .TRUE., S (and also P) contains the value of
C             the selected real eigenvalue.
C             If REIG = .FALSE., S and P contain the sum and product,
C             respectively, of the selected complex conjugate pair of
C             eigenvalues.
C
C     FURTHER COMMENTS
C
C     For efficiency reasons, |x| + |y| is used for a complex number
C     x + jy, instead of its modulus.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     February 1999. Based on the RASP routine PMDIST.
C
C     REVISIONS
C
C     March 30, 1999, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     Feb. 15, 2004, V. Sima, Research Institute for Informatics,
C     Bucharest.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      LOGICAL          REIG
      INTEGER          N
      DOUBLE PRECISION P, S, XI ,XR
C     .. Array Arguments ..
      DOUBLE PRECISION WI(*), WR(*)
C     .. Local Scalars ..
      INTEGER          I, J, K
      DOUBLE PRECISION X, Y
C     .. Intrinsic Functions ..
      INTRINSIC        ABS
C     .. Executable Statements ..
C
      J = 1
      IF( REIG ) THEN
         Y = ABS( WR(1)-XR )
         DO 10 I = 2, N
            X = ABS( WR(I)-XR )
            IF( X .LT. Y ) THEN
               Y = X
               J = I
            END IF
   10    CONTINUE
         S = WR(J)
         K = N - J
         IF( K .GT. 0 ) THEN
            DO 20 I = J, J + K - 1
               WR(I) = WR(I+1)
   20       CONTINUE
            WR(N) = S
         END IF
         P = S
      ELSE
         Y = ABS( WR(1)-XR ) + ABS( WI(1)-XI )
         DO 30 I = 3, N, 2
            X = ABS( WR(I)-XR ) + ABS( WI(I)-XI )
            IF( X .LT. Y ) THEN
               Y = X
               J = I
            END IF
   30    CONTINUE
         X = WR(J)
         Y = WI(J)
         K = N - J - 1
         IF( K .GT. 0 ) THEN
            DO 40 I = J, J + K - 1
               WR(I) = WR(I+2)
               WI(I) = WI(I+2)
   40       CONTINUE
            WR(N-1) = X
            WI(N-1) = Y
            WR(N) = X
            WI(N) = -Y
         END IF
         S = X + X
         P = X * X + Y * Y
      END IF
C
      RETURN
C *** End of SB01BX ***
      END
