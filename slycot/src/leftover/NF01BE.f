      SUBROUTINE NF01BE( IFLAG, NSMP, N, IPAR, LIPAR, Z, LDZ, Y, LDY, X,
     $                   NFEVL, E, J, LDJ, DWORK, LDWORK, INFO )
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
C     This is the FCN routine for optimizing the parameters of the
C     nonlinear part of a Wiener system (initialization phase), using
C     SLICOT Library routine MD03BD. See the argument FCN in the
C     routine MD03BD for the description of parameters. Note that
C     NF01BE is called for each output of the Wiener system.
C
C     ******************************************************************
C
C     .. Parameters ..
C     .. CJTE is initialized to avoid the calculation of J'*e ..
C     .. NOUT is the unit number for printing intermediate results ..
      CHARACTER         CJTE
      PARAMETER         ( CJTE = 'N' )
      INTEGER           NOUT
      PARAMETER         ( NOUT = 6 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           IFLAG, INFO, LDJ, LDWORK, LDY, LDZ, LIPAR, N,
     $                  NFEVL, NSMP
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DWORK(*), E(*), J(LDJ,*), X(*), Y(LDY,*),
     $                  Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ERR
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, NF01AY, NF01BY
C
C     .. Executable Statements ..
C
      INFO = 0
      IF ( IFLAG.EQ.1 ) THEN
C
C        Call NF01AY to compute the output y of the Wiener system (in E)
C        and then the error functions (also in E). The array Z must
C        contain the output of the linear part of the Wiener system, and
C        Y must contain the original output Y of the Wiener system.
C        IPAR(2) must contain the number of outputs.
C        Workspace: need:    2*NN, NN = IPAR(3) (number of neurons);
C                   prefer:  larger.
C
         CALL NF01AY( NSMP, IPAR(2), 1, IPAR(3), LIPAR-2, X, N, Z, LDZ,
     $                E, NSMP, DWORK, LDWORK, INFO )
         CALL DAXPY( NSMP, -ONE, Y, 1, E, 1 )
         DWORK(1) = 2*IPAR(3)
C
      ELSE IF ( IFLAG.EQ.2 ) THEN
C
C        Call NF01BY to compute the Jacobian in a compressed form.
C        IPAR(2), IPAR(3) must have the same content as for IFLAG = 1.
C        Workspace: need:    0.
C
         CALL NF01BY( CJTE, NSMP, IPAR(2), 1, IPAR(3), LIPAR-2, X, N, Z,
     $                LDZ, E, J, LDJ, DWORK, DWORK, LDWORK, INFO )
         NFEVL = 0
         DWORK(1) = ZERO
C
      ELSE IF ( IFLAG.EQ.3 ) THEN
C
C        Set the parameter LDJ, the length of the array J, and the sizes
C        of the workspace for FCN (IFLAG = 1 or 2), QRFACT and LMPARM.
C
         LDJ     = NSMP
         IPAR(1) = NSMP*N
         IPAR(2) = 2*IPAR(3)
         IPAR(3) = 0
         IPAR(4) = 4*N + 1
         IPAR(5) = 4*N
C
      ELSE IF ( IFLAG.EQ.0 ) THEN
C
C        Special call for printing intermediate results.
C
         ERR = DNRM2( NSMP, E, 1 )
         WRITE( NOUT, '('' Norm of current error = '', D15.6)') ERR
      END IF
      RETURN
C
C *** Last line of NF01BE ***
      END
