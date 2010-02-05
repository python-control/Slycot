      SUBROUTINE NF01AY( NSMP, NZ, L, IPAR, LIPAR, WB, LWB, Z, LDZ,
     $                   Y, LDY, DWORK, LDWORK, INFO )
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
C     To calculate the output of a set of neural networks with the
C     structure
C
C             - tanh(w1'*z+b1) -
C           /      :             \
C         z ---    :           --- sum(ws(i)*...)+ b(n+1)  --- y,
C           \      :             /
C             - tanh(wn'*z+bn) -
C
C     given the input z and the parameter vectors wi, ws, and b,
C     where z, w1, ..., wn are vectors of length NZ, ws is a vector
C     of length n, b(1), ..., b(n+1) are scalars, and n is called the
C     number of neurons in the hidden layer, or just number of neurons.
C     Such a network is used for each L output variables.
C
C     ARGUMENTS
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
C             The length of each output sample.  L >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters needed.
C             IPAR(1) must contain the number of neurons, n, per output
C             variable, denoted NN in the sequel.  NN >= 0.
C
C     LIPAR   (input) INTEGER
C             The length of the vector IPAR.  LIPAR >= 1.
C
C     WB      (input) DOUBLE PRECISION array, dimension (LWB)
C             The leading (NN*(NZ+2)+1)*L part of this array must
C             contain the weights and biases of the network. This vector
C             is partitioned into L vectors of length NN*(NZ+2)+1,
C             WB = [ wb(1), ..., wb(L) ]. Each wb(k), k = 1, ..., L,
C             corresponds to one output variable, and has the structure
C             wb(k) = [ w1(1), ..., w1(NZ), ..., wn(1), ..., wn(NZ),
C                       ws(1), ..., ws(n), b(1), ..., b(n+1) ],
C             where wi(j) are the weights of the hidden layer,
C             ws(i) are the weights of the linear output layer, and
C             b(i) are the biases, as in the scheme above.
C
C     LWB     (input) INTEGER
C             The length of the array WB.
C             LWB >= ( NN*(NZ + 2) + 1 )*L.
C
C     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ)
C             The leading NSMP-by-NZ part of this array must contain the
C             set of input samples,
C             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ).
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= MAX(1,NSMP).
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY, L)
C             The leading NSMP-by-L part of this array contains the set
C             of output samples,
C             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ).
C
C     LDY     INTEGER
C             The leading dimension of the array Y.  LDY >= MAX(1,NSMP).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 2*NN.
C             For better performance, LDWORK should be larger.
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
C     BLAS routines are used to compute the matrix-vector products.
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
C     simulation, system response.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDWORK, LDY, LDZ, LIPAR, LWB, NSMP, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), WB(*), Y(LDY,*), Z(LDZ,*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      LOGICAL           LAST
      INTEGER           I, IB, J, K, LDWB, LJ, LK, M, MF, NN, NV, WS
      DOUBLE PRECISION  BIGNUM, DF, SMLNUM, TMP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH
      EXTERNAL          DDOT, DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLABAD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, MIN, MOD
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      NN   = IPAR(1)
      LDWB = NN*( NZ + 2 ) + 1
      IF ( NSMP.LT.0 ) THEN
         INFO = -1
      ELSEIF ( NZ.LT.0 ) THEN
         INFO = -2
      ELSEIF ( L.LT.0 ) THEN
         INFO = -3
      ELSEIF ( NN.LT.0 ) THEN
         INFO = -4
      ELSEIF ( LIPAR.LT.1 ) THEN
         INFO = -5
      ELSEIF ( LWB.LT.LDWB*L ) THEN
         INFO = -7
      ELSEIF ( LDZ.LT.MAX( 1, NSMP ) ) THEN
         INFO = -9
      ELSEIF ( LDY.LT.MAX( 1, NSMP ) ) THEN
         INFO = -11
      ELSEIF ( LDWORK.LT.2*NN ) THEN
          INFO = -13
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01AY', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MIN( NSMP, L ).EQ.0 )
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
      WS = NZ*NN + 1
      IB = WS + NN - 1
      LK = 0
      IF ( MIN( NZ, NN ).EQ.0 ) THEN
         NV = 2
      ELSE
         NV = ( LDWORK - NN )/NN
      END IF
C
      IF ( NV.GT.2 ) THEN
         MF   = ( NSMP/NV )*NV
         LAST = MOD( NSMP, NV ).NE.0
C
C        Some BLAS 3 calculations can be used.
C
         DO 70 K = 0, L - 1
            TMP = WB(IB+NN+1+LK)
C
            DO 10 J = 1, NN
               DWORK(J) = TWO*WB(IB+J+LK)
   10       CONTINUE
C
            DO 40 I = 1, MF, NV
C
C              Compute -2*[w1 w2 ... wn]'*Z', where
C              Z = [z(i)';...; z(i+NV-1)'].
C
               CALL DGEMM( 'Transpose', 'Transpose', NN, NV, NZ, -TWO,
     $                     WB(1+LK), NZ, Z(I,1), LDZ, ZERO, DWORK(NN+1),
     $                     NN )
               LJ = NN
C
               DO 30 M = 1, NV
                  DO 20 J = 1, NN
C
C                    Compute tanh(wj'*z(i) + bj), j = 1:n.
C
                     LJ = LJ + 1
                     DF = DWORK(LJ) - DWORK(J)
                     IF ( ABS( DF ).GE.BIGNUM ) THEN
                        IF ( DF.GT.ZERO ) THEN
                           DWORK(LJ) = -ONE
                        ELSE
                           DWORK(LJ) = ONE
                        END IF
                     ELSE IF ( ABS( DF ).LE.SMLNUM ) THEN
                        DWORK(LJ) = ZERO
                     ELSE
                        DWORK(LJ) = TWO/( ONE + EXP( DF ) ) - ONE
                     END IF
   20             CONTINUE
C
   30          CONTINUE
C
               Y(I, K+1) = TMP
               CALL DCOPY( NV-1, Y(I, K+1), 0, Y(I+1, K+1), 1 )
               CALL DGEMV( 'Transpose', NN, NV, ONE, DWORK(NN+1), NN,
     $                     WB(WS+LK), 1, ONE, Y(I, K+1), 1 )
   40       CONTINUE
C
            IF ( LAST ) THEN
C
C              Process the last samples.
C
               NV = NSMP - MF
               I  = MF + 1
C
C              Compute -2*[w1 w2 ... wn]'*Z', where
C              Z = [z(i)';...; z(NSMP)'].
C
               CALL DGEMM( 'Transpose', 'Transpose', NN, NV, NZ, -TWO,
     $                     WB(1+LK), NZ, Z(I,1), LDZ, ZERO, DWORK(NN+1),
     $                     NN )
               LJ = NN
C
               DO 60 M = 1, NV
                  DO 50 J = 1, NN
C
C                    Compute tanh(wj'*z(i) + bj), j = 1:n.
C
                     LJ = LJ + 1
                     DF = DWORK(LJ) - DWORK(J)
                     IF ( ABS( DF ).GE.BIGNUM ) THEN
                        IF ( DF.GT.ZERO ) THEN
                           DWORK(LJ) = -ONE
                        ELSE
                           DWORK(LJ) = ONE
                        END IF
                     ELSE IF ( ABS( DF ).LE.SMLNUM ) THEN
                        DWORK(LJ) = ZERO
                     ELSE
                        DWORK(LJ) = TWO/( ONE + EXP( DF ) ) - ONE
                     END IF
   50             CONTINUE
C
   60          CONTINUE
C
               Y(I, K+1) = TMP
               IF ( NV.GT.1 )
     $            CALL DCOPY( NV-1, Y(I, K+1), 0, Y(I+1, K+1), 1 )
               CALL DGEMV( 'Transpose', NN, NV, ONE, DWORK(NN+1), NN,
     $                     WB(WS+LK), 1, ONE, Y(I, K+1), 1 )
            END IF
C
            LK = LK + LDWB
   70    CONTINUE
C
      ELSE
C
C        BLAS 2 calculations only can be used.
C
         DO 110 K = 0, L - 1
            TMP = WB(IB+NN+1+LK)
C
            DO 80 J = 1, NN
               DWORK(J) = TWO*WB(IB+J+LK)
   80       CONTINUE
C
            DO 100 I = 1, NSMP
C
C              Compute -2*[w1 w2 ... wn]'*z(i).
C
               IF ( NZ.EQ.0 ) THEN
                  DWORK(NN+1) = ZERO
                  CALL DCOPY( NN, DWORK(NN+1), 0, DWORK(NN+1), 1 )
               ELSE
                  CALL DGEMV( 'Transpose', NZ, NN, -TWO, WB(1+LK), NZ,
     $                        Z(I,1), LDZ, ZERO, DWORK(NN+1), 1 )
               END IF
C
               DO 90 J = NN + 1, 2*NN
C
C                 Compute tanh(wj'*z(i) + bj), j = 1:n.
C
                  DF = DWORK(J) - DWORK(J-NN)
                  IF ( ABS( DF ).GE.BIGNUM ) THEN
                     IF ( DF.GT.ZERO ) THEN
                        DWORK(J) = -ONE
                     ELSE
                        DWORK(J) = ONE
                     END IF
                  ELSE IF ( ABS( DF ).LE.SMLNUM ) THEN
                     DWORK(J) = ZERO
                  ELSE
                     DWORK(J) = TWO/( ONE + EXP( DF ) ) - ONE
                  END IF
   90          CONTINUE
C
               Y(I, K+1) = DDOT( NN, WB(WS+LK), 1, DWORK(NN+1), 1 ) +
     $                     TMP
  100       CONTINUE
C
            LK = LK + LDWB
  110    CONTINUE
C
      END IF
      RETURN
C
C *** Last line of NF01AY ***
      END
