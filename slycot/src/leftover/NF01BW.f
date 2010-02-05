      SUBROUTINE NF01BW( N, IPAR, LIPAR, DPAR, LDPAR, J, LDJ, X, INCX,
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
C     To compute the matrix-vector product x <-- (J'*J + c*I)*x, for the
C     Jacobian J as received from SLICOT Library routine NF01BD:
C
C          /  dy(1)/dwb(1)  |  dy(1)/dtheta  \
C     Jc = |       :        |       :        | .
C          \  dy(L)/dwb(L)  |  dy(L)/dtheta  /
C
C     This is a compressed representation of the actual structure
C
C         /   J_1    0    ..   0   |  L_1  \
C         |    0    J_2   ..   0   |  L_2  |
C     J = |    :     :    ..   :   |   :   | .
C         |    :     :    ..   :   |   :   |
C         \    0     0    ..  J_L  |  L_L  /
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the vector x.
C             N = BN*BSN + ST >= 0.  (See parameter description below.)
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters describing the structure of the
C             matrix J, as follows:
C             IPAR(1) must contain ST, the number of parameters
C                     corresponding to the linear part.  ST >= 0.
C             IPAR(2) must contain BN, the number of blocks, BN = L,
C                     for the parameters corresponding to the nonlinear
C                     part.  BN >= 0.
C             IPAR(3) must contain BSM, the number of rows of the blocks
C                     J_k = dy(k)/dwb(k), k = 1:BN, if BN > 0, or the
C                     number of rows of the matrix J, if BN <= 1.
C             IPAR(4) must contain BSN, the number of columns of the
C                     blocks J_k, k = 1:BN.  BSN >= 0.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 4.
C
C     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C             The real parameters needed for solving the problem.
C             The entry DPAR(1) must contain the real scalar c.
C
C     LDPAR   (input) INTEGER
C             The length of the array DPAR.  LDPAR >= 1.
C
C     J       (input) DOUBLE PRECISION array, dimension (LDJ, NC)
C             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1.
C             The leading NR-by-NC part of this array must contain
C             the (compressed) representation (Jc) of the Jacobian
C             matrix J, where NR = BSM if BN <= 1, and NR = BN*BSM,
C             if BN > 1.
C
C     LDJ     (input) INTEGER
C             The leading dimension of array J.  LDJ >= MAX(1,NR).
C
C     X       (input/output) DOUBLE PRECISION array, dimension
C             (1+(N-1)*INCX)
C             On entry, this incremented array must contain the
C             vector x.
C             On exit, this incremented array contains the value of the
C             matrix-vector product (J'*J + c*I)*x.
C
C     INCX    (input) INTEGER
C             The increment for the elements of X.  INCX >= 1.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= NR.
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
C     The associativity of matrix multiplications is used; the result
C     is obtained as:  x_out = J'*( J*x ) + c*x.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Mar. 2001, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001,
C     Mar. 2002.
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations,
C     Wiener system.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INCX, INFO, LDJ, LDPAR, LDWORK, LIPAR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DPAR(*), DWORK(*), J(LDJ,*), X(*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      INTEGER           BN, BSM, BSN, IBSM, IBSN, IX, JL, M, NTHS, ST,
     $                  XL
      DOUBLE PRECISION  C
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSEIF ( LIPAR.LT.4 ) THEN
         INFO = -3
      ELSEIF ( LDPAR.LT.1 ) THEN
         INFO = -5
      ELSEIF ( INCX.LT.1 ) THEN
         INFO = -9
      ELSE
         ST   = IPAR(1)
         BN   = IPAR(2)
         BSM  = IPAR(3)
         BSN  = IPAR(4)
         NTHS = BN*BSN
         IF ( BN.GT.1 ) THEN
            M = BN*BSM
         ELSE
            M = BSM
         END IF
         IF ( MIN( ST, BN, BSM, BSN ).LT.0 ) THEN
            INFO = -2
         ELSEIF ( N.NE.NTHS + ST ) THEN
            INFO = -1
         ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
            INFO = -7
         ELSEIF ( LDWORK.LT.M ) THEN
            INFO = -11
         END IF
      END IF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BW', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      C = DPAR(1)
C
      IF ( M.EQ.0 ) THEN
C
C        Special case, void Jacobian: x <-- c*x.
C
         CALL DSCAL( N, C, X, INCX )
         RETURN
      END IF
C
      IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
C
C        Special case, l <= 1 or BSN = 0: the Jacobian is represented
C        as a full matrix. Adapted code from NF01BX is included in-line.
C
         CALL DGEMV( 'NoTranspose', M, N, ONE, J, LDJ, X, INCX, ZERO,
     $               DWORK, 1 )
         CALL DGEMV( 'Transpose', M, N, ONE, J, LDJ, DWORK, 1, C, X,
     $               INCX )
         RETURN
      END IF
C
C     General case: l > 1, BSN > 0, BSM > 0.
C
      JL = BSN + 1
      IX = BSN*INCX
      XL = BN*IX + 1
C
      IF ( ST.GT.0 ) THEN
         CALL DGEMV( 'NoTranspose', M, ST, ONE, J(1,JL), LDJ, X(XL),
     $               INCX, ZERO, DWORK, 1 )
      ELSE
         DWORK(1) = ZERO
         CALL DCOPY( M, DWORK(1), 0, DWORK, 1 )
      END IF
      IBSN = 1
C
      DO 10 IBSM = 1, M, BSM
         CALL DGEMV( 'NoTranspose', BSM, BSN, ONE, J(IBSM,1), LDJ,
     $               X(IBSN), INCX, ONE, DWORK(IBSM), 1 )
         CALL DGEMV( 'Transpose', BSM, BSN, ONE, J(IBSM,1), LDJ,
     $               DWORK(IBSM), 1, C, X(IBSN), INCX )
         IBSN = IBSN + IX
   10 CONTINUE
C
      IF ( ST.GT.0 )
     $   CALL DGEMV( 'Transpose', M, ST, ONE, J(1,JL), LDJ, DWORK, 1, C,
     $               X(XL), INCX )
C
      RETURN
C
C *** Last line of NF01BW ***
      END
