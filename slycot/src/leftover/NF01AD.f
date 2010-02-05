      SUBROUTINE NF01AD( NSMP, M, L, IPAR, LIPAR, X, LX, U, LDU, Y, LDY,
     $                   DWORK, LDWORK, INFO )
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
C     To calculate the output y of the Wiener system
C
C        x(t+1) = A*x(t) + B*u(t)
C        z(t)   = C*x(t) + D*u(t),
C
C        y(t)   = f(z(t),wb(1:L)),
C
C     where t = 1, 2, ..., NSMP, and f is a nonlinear function,
C     evaluated by the SLICOT Library routine NF01AY. The parameter
C     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ),
C     where wb(i), i = 1:L, correspond to the nonlinear part, theta
C     corresponds to the linear part, and the notation is fully
C     described below.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NSMP    (input) INTEGER
C             The number of training samples.  NSMP >= 0.
C
C     M       (input) INTEGER
C             The length of each input sample.  M >= 0.
C
C     L       (input) INTEGER
C             The length of each output sample.  L >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters needed.
C             IPAR(1)  must contain the order of the linear part,
C                      referred to as N below.  N >= 0.
C             IPAR(2)  must contain the number of neurons for the
C                      nonlinear part, referred to as NN below.
C                      NN >= 0.
C
C     LIPAR   (input) INTEGER
C             The length of IPAR.  LIPAR >= 2.
C
C     X       (input) DOUBLE PRECISION array, dimension (LX)
C             The parameter vector, partitioned as
C             X = (wb(1), ..., wb(L), theta), where the vectors
C             wb(i), of length NN*(L+2)+1, are parameters for the
C             static nonlinearity, which is simulated by the
C             SLICOT Library routine NF01AY. See the documentation of
C             NF01AY for further details. The vector theta, of length
C             N*(M + L + 1) + L*M, represents the matrices A, B, C,
C             D and x(1), and it can be retrieved from these matrices
C             by SLICOT Library routine TB01VD and retranslated by
C             TB01VY.
C
C     LX      (input) INTEGER
C             The length of the array X.
C             LX >= ( NN*(L+2)+1 )*L + N*(M + L + 1) + L*M.
C
C     U       (input) DOUBLE PRECISION array, dimension (LDU, M)
C             The leading NSMP-by-M part of this array must contain the
C             set of input samples,
C             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ).
C
C     LDU     INTEGER
C             The leading dimension of the array U.  LDU >= MAX(1,NSMP).
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY, L)
C             The leading NSMP-by-L part of this array contains the
C             simulated output.
C
C     LDY     INTEGER
C             The leading dimension of the array Y.  LDY >= MAX(1,NSMP).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N +
C                                     MAX( N*(N + L), N + M + L ) )
C                                                              if M > 0;
C             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*N + 2*N +
C                                     MAX( N*(N + L), L ) ),   if M = 0.
C             A larger value of LDWORK could improve the efficiency.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C     METHOD
C
C     BLAS routines are used for the matrix-vector multiplications and
C     the routine NF01AY is called for the calculation of the nonlinear
C     function.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Mar. 2001, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Dec. 2001.
C
C     KEYWORDS
C
C     Nonlinear system, output normal form, simulation, state-space
C     representation, Wiener system.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDU, LDWORK, LDY, LX, LIPAR, M, NSMP
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DWORK(*), U(LDU,*), X(*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           AC, BD, IX, JW, LDAC, LTHS, N, NN, NTHS, Z
C     .. External Subroutines ..
      EXTERNAL          NF01AY, TB01VY, TF01MX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      IF ( NSMP.LT.0 ) THEN
         INFO = -1
      ELSEIF ( M.LT.0 ) THEN
         INFO = -2
      ELSEIF ( L.LT.0 ) THEN
         INFO = -3
      ELSEIF ( LIPAR.LT.2 ) THEN
         INFO = -5
      ELSE
C
         N    = IPAR(1)
         NN   = IPAR(2)
         LDAC = N + L
         NTHS = ( NN*( L + 2 ) + 1 )*L
         LTHS = N*( M + L + 1 ) + L*M
C
         IF ( N.LT.0 .OR. NN.LT.0 ) THEN
            INFO = -4
         ELSEIF ( LX.LT.NTHS + LTHS ) THEN
            INFO = -7
         ELSEIF ( LDU.LT.MAX( 1, NSMP ) ) THEN
            INFO = -9
         ELSEIF ( LDY.LT.MAX( 1, NSMP ) ) THEN
            INFO = -11
         ELSE
            IF ( M.GT.0 ) THEN
               JW = MAX( N*LDAC, N + M + L )
            ELSE
               JW = MAX( N*LDAC, L )
            END IF
            IF ( LDWORK.LT.NSMP*L + MAX( 2*NN, LDAC*( N + M ) + 2*N +
     $            JW ) )
     $         INFO = -13
         ENDIF
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01AD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MIN( NSMP, L ).EQ.0 )
     $   RETURN
C
C     Compute the output of the linear part.
C     Workspace: need   NSMP*L + (N + L)*(N + M) + N + N*(N + L + 1).
C     (NSMP*L locations are reserved for the output of the linear part.)
C
      Z  = 1
      AC = Z  + NSMP*L
      BD = AC + LDAC*N
      IX = BD + LDAC*M
      JW = IX + N
C
      CALL TB01VY( 'Apply', N, M, L, X(NTHS+1), LTHS, DWORK(AC), LDAC,
     $             DWORK(BD), LDAC, DWORK(AC+N), LDAC, DWORK(BD+N),
     $             LDAC, DWORK(IX), DWORK(JW), LDWORK-JW+1, INFO )
C
C     Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, if M>0;
C                       NSMP*L + (N + L)*N + 2*N + L,           if M=0;
C                prefer larger.
C
      CALL TF01MX( N, M, L, NSMP, DWORK(AC), LDAC, U, LDU, DWORK(IX),
     $             DWORK(Z), NSMP, DWORK(JW), LDWORK-JW+1, INFO )
C
C     Simulate the static nonlinearity.
C     Workspace: need   NSMP*L + 2*NN;
C                prefer larger.
C
      JW = AC
      CALL NF01AY( NSMP, L, L, IPAR(2), LIPAR-1, X, NTHS, DWORK(Z),
     $             NSMP, Y, LDY, DWORK(JW), LDWORK-JW+1, INFO )
C
      RETURN
C
C *** Last line of NF01AD ***
      END
