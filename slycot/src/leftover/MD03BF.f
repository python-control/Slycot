      SUBROUTINE MD03BF( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, DPAR2,
     $                   LDPAR2, X, NFEVL, E, J, LDJ, DWORK, LDWORK,
     $                   INFO )
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
C     This is the FCN routine for solving a standard nonlinear least
C     squares problem using SLICOT Library routine MD03BD. See the
C     parameter FCN in the routine MD03BD for the description of
C     parameters.
C
C     The example programmed in this routine is adapted from that
C     accompanying the MINPACK routine LMDER.
C
C     ******************************************************************
C
C     .. Parameters ..
C     .. NOUT is the unit number for printing intermediate results ..
      INTEGER           NOUT
      PARAMETER         ( NOUT = 6 )
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           IFLAG, INFO, LDJ, LDPAR1, LDPAR2, LDWORK, LIPAR,
     $                  M, N, NFEVL
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DPAR1(*), DPAR2(*), DWORK(*), E(*), J(LDJ,*),
     $                  X(*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  ERR, TMP1, TMP2, TMP3, TMP4
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. DATA Statements ..
      DOUBLE PRECISION  Y(15)
      DATA              Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), Y(7), Y(8),
     $                  Y(9), Y(10), Y(11), Y(12), Y(13), Y(14), Y(15)
     $                  / 1.4D-1, 1.8D-1, 2.2D-1, 2.5D-1, 2.9D-1,
     $                    3.2D-1, 3.5D-1, 3.9D-1, 3.7D-1, 5.8D-1,
     $                    7.3D-1, 9.6D-1, 1.34D0, 2.1D0,  4.39D0 /
C
C     .. Executable Statements ..
C
      INFO = 0
      IF ( IFLAG.EQ.1 ) THEN
C
C        Compute the error function values.
C
         DO 10 I = 1, 15
            TMP1 = I
            TMP2 = 16 - I
            IF ( I.GT.8 ) THEN
               TMP3 = TMP2
            ELSE
               TMP3 = TMP1
            END IF
            E(I) = Y(I) - ( X(1) + TMP1/( X(2)*TMP2 + X(3)*TMP3 ) )
   10    CONTINUE
C
      ELSE IF ( IFLAG.EQ.2 ) THEN
C
C        Compute the Jacobian.
C
         DO 30 I = 1, 15
            TMP1 = I
            TMP2 = 16 - I
            IF ( I.GT.8 ) THEN
               TMP3 = TMP2
            ELSE
               TMP3 = TMP1
            END IF
            TMP4 = ( X(2)*TMP2 + X(3)*TMP3 )**2
            J(I,1) = -ONE
            J(I,2) = TMP1*TMP2/TMP4
            J(I,3) = TMP1*TMP3/TMP4
   30    CONTINUE
C
         NFEVL = 0
C
      ELSE IF ( IFLAG.EQ.3 ) THEN
C
C        Set the parameter LDJ, the length of the array J, and the sizes
C        of the workspace for FCN (IFLAG = 1 or 2), MD03BA and MD03BB.
C
         LDJ = M
         IPAR(1) = M*N
         IPAR(2) = 0
         IPAR(3) = 0
         IPAR(4) = 4*N + 1
         IPAR(5) = 4*N
C
      ELSE IF ( IFLAG.EQ.0 ) THEN
C
C        Special call for printing intermediate results.
C
         ERR = DNRM2( M, E, 1 )
         WRITE( 1, '('' Norm of current error = '', D15.6)') ERR
C
      END IF
C
      RETURN
C
C *** Last line of MD03BF ***
      END
