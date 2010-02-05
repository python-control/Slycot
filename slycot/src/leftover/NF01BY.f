      SUBROUTINE NF01BY( CJTE, NSMP, NZ, L, IPAR, LIPAR, WB, LWB, Z,
     $                   LDZ, E, J, LDJ, JTE, DWORK, LDWORK, INFO )
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
C     To compute the Jacobian of the error function for a neural network
C     of the structure
C
C             - tanh(w1*z+b1) -
C           /      :            \
C         z ---    :          --- sum(ws(i)*...)+ b(n+1)  --- y,
C           \      :            /
C             - tanh(wn*z+bn) -
C
C     for the single-output case. The Jacobian has the form
C
C                d e(1)  / d WB(1)   ...    d e(1)  / d WB(NWB)
C         J =            :                          :           ,
C              d e(NSMP) / d WB(1)   ...  d e(NSMP) / d WB(NWB)
C
C     where e(z) is the error function, WB is the set of weights and
C     biases of the network (for the considered output), and NWB is
C     the number of elements of this set, NWB = IPAR(1)*(NZ+2)+1
C     (see below).
C
C     In the multi-output case, this routine should be called for each
C     output.
C
C     NOTE: this routine must have the same arguments as SLICOT Library
C     routine NF01BD.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     CJTE    CHARACTER*1
C             Specifies whether the matrix-vector product J'*e should be
C             computed or not, as follows:
C             = 'C' :  compute J'*e;
C             = 'N' :  do not compute J'*e.
C
C     Input/Output Parameters
C
C     NSMP    (input) INTEGER
C             The number of training samples.  NSMP >= 0.
C
C     NZ      (input) INTEGER
C             The length of each input sample.  NZ >= 0.
C
C     L       (input) INTEGER
C             The length of each output sample.
C             Currently, L must be 1.
C
C     IPAR    (input/output) INTEGER array, dimension (LIPAR)
C             The integer parameters needed.
C             On entry, the first element of this array must contain
C             a value related to the number of neurons, n; specifically,
C             n = abs(IPAR(1)), since setting IPAR(1) < 0 has a special
C             meaning (see below).
C             On exit, if IPAR(1) < 0 on entry, then no computations are
C             performed, except the needed tests on input parameters,
C             but the following values are returned:
C             IPAR(1) contains the length of the array J, LJ;
C             LDJ     contains the leading dimension of array J.
C             Otherwise, IPAR(1) and LDJ are unchanged on exit.
C
C     LIPAR   (input) INTEGER
C             The length of the vector IPAR.  LIPAR >= 1.
C
C     WB      (input) DOUBLE PRECISION array, dimension (LWB)
C             The leading NWB = IPAR(1)*(NZ+2)+1 part of this array
C             must contain the weights and biases of the network,
C             WB = ( w(1,1), ..., w(1,NZ), ..., w(n,1), ...,  w(n,NZ),
C                    ws(1), ..., ws(n), b(1), ..., b(n+1) ),
C             where w(i,j) are the weights of the hidden layer,
C             ws(i) are the weights of the linear output layer and
C             b(i) are the biases.
C
C     LWB     (input) INTEGER
C             The length of array WB.  LWB >= NWB.
C
C     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ)
C             The leading NSMP-by-NZ part of this array must contain the
C             set of input samples,
C             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ).
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= MAX(1,NSMP).
C
C     E       (input) DOUBLE PRECISION array, dimension (NSMP)
C             If CJTE = 'C', this array must contain the error vector e.
C             If CJTE = 'N', this array is not referenced.
C
C     J       (output) DOUBLE PRECISION array, dimension (LDJ, NWB)
C             The leading NSMP-by-NWB part of this array contains the
C             Jacobian of the error function.
C
C     LDJ     INTEGER
C             The leading dimension of array J.  LDJ >= MAX(1,NSMP).
C             Note that LDJ is an input parameter, except for
C             IPAR(1) < 0 on entry, when it is an output parameter.
C
C     JTE     (output) DOUBLE PRECISION array, dimension (NWB)
C             If CJTE = 'C', this array contains the matrix-vector
C             product J'*e.
C             If CJTE = 'N', this array is not referenced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             This argument is included for combatibility with SLICOT
C             Library routine NF01BD.
C
C     LDWORK  INTEGER
C             Normally, the length of the array DWORK.  LDWORK >= 0.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The Jacobian is computed analytically.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Input output description, neural network, nonlinear system,
C     optimization, system response.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         CJTE
      INTEGER           INFO, L, LDJ, LDWORK, LDZ, LIPAR, LWB, NSMP, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), E(*), J(LDJ,*), JTE(*), WB(*),
     $                  Z(LDZ,*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      LOGICAL           WJTE
      INTEGER           BP1, DI, I, IB, K, M, NN, NWB, WS
      DOUBLE PRECISION  BIGNUM, SMLNUM, TMP
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLABAD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      WJTE = LSAME( CJTE, 'C' )
      INFO = 0
      NN   = IPAR(1)
      NWB  = NN*( NZ + 2 ) + 1
      IF( .NOT.( WJTE .OR. LSAME( CJTE, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( NSMP.LT.0 ) THEN
         INFO = -2
      ELSEIF ( NZ.LT.0 ) THEN
         INFO = -3
      ELSEIF ( L.NE.1 ) THEN
         INFO = -4
      ELSEIF ( LIPAR.LT.1 ) THEN
         INFO = -6
      ELSEIF ( IPAR(1).LT.0 ) THEN
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'NF01BY', -INFO )
         ELSE
            IPAR(1) = NSMP*( ABS( NN )*( NZ + 2 ) + 1 )
            LDJ     = NSMP
         ENDIF
         RETURN
      ELSEIF ( LWB.LT.NWB ) THEN
         INFO = -8
      ELSEIF ( LDZ.LT.MAX( 1, NSMP ) ) THEN
         INFO = -10
      ELSEIF ( LDJ.LT.MAX( 1, NSMP ) ) THEN
         INFO = -13
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BY', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MIN( NSMP, NZ ).EQ.0 )
     $   RETURN
C
C     Set parameters to avoid overflows and increase accuracy for
C     extreme values.
C
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = LOG( SMLNUM )
      BIGNUM = LOG( BIGNUM )
C
      WS  = NZ*NN + 1
      IB  = WS + NN
      BP1 = IB + NN
C
      J(1, BP1) = ONE
      CALL DCOPY( NSMP, J(1, BP1), 0, J(1, BP1), 1 )
C
      DO 10 I = 0, NN - 1
         CALL DCOPY( NSMP, WB(IB+I), 0, J(1, WS+I), 1 )
   10 CONTINUE
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NSMP, NN, NZ, -TWO, Z,
     $            LDZ, WB, NZ, -TWO, J(1, WS), LDJ )
      DI = 1
C
      DO 50 I = 0, NN - 1
C
         DO 20 K = 1, NSMP
            TMP = J(K, WS+I)
            IF ( ABS( TMP ).GE.BIGNUM ) THEN
               IF ( TMP.GT.ZERO ) THEN
                  J(K, WS+I) = -ONE
               ELSE
                  J(K, WS+I) = ONE
               END IF
            ELSE IF ( ABS( TMP ).LE.SMLNUM ) THEN
               J(K, WS+I) = ZERO
            ELSE
               J(K, WS+I) = TWO/( ONE + EXP( TMP ) ) - ONE
            END IF
            J(K, IB+I) = WB(WS+I)*( ONE - J(K, WS+I)**2 )
   20    CONTINUE
C
         DO 40 K = 0, NZ - 1
C
             DO 30 M = 1, NSMP
                J(M, DI+K) = J(M, IB+I)*Z(M, K+1)
   30       CONTINUE
C
   40    CONTINUE
C
         DI = DI + NZ
   50 CONTINUE
C
      IF ( WJTE ) THEN
C
C        Compute J'e.
C
         CALL DGEMV( 'Transpose', NSMP, NWB, ONE, J, LDJ, E, 1, ZERO,
     $               JTE, 1 )
      END IF
C
      RETURN
C
C *** Last line of NF01BY ***
      END
