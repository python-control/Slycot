      SUBROUTINE NF01BU( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, J,
     $                   LDJ, JTJ, LDJTJ, DWORK, LDWORK, INFO )
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
C     To compute the matrix J'*J + c*I, for the Jacobian J as received
C     from SLICOT Library routine NF01BD:
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
C     Mode Parameters
C
C     STOR    CHARACTER*1
C             Specifies the storage scheme for the symmetric
C             matrix J'*J + c*I, as follows:
C             = 'F' :  full storage is used;
C             = 'P' :  packed storage is used.
C
C     UPLO    CHARACTER*1
C             Specifies which part of the matrix J'*J + c*I is stored,
C             as follows:
C             = 'U' :  the upper triagular part is stored;
C             = 'L' :  the lower triagular part is stored.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix J'*J + c*I.
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
C     JTJ     (output) DOUBLE PRECISION array,
C                      dimension (LDJTJ,N),    if STOR = 'F',
C                      dimension (N*(N+1)/2),  if STOR = 'P'.
C             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if
C             STOR = 'P') part of this array contains the upper or
C             lower triangle of the matrix J'*J + c*I, depending on
C             UPLO = 'U', or UPLO = 'L', respectively, stored either as
C             a two-dimensional, or one-dimensional array, depending
C             on STOR.
C
C     LDJTJ   INTEGER
C             The leading dimension of the array JTJ.
C             LDJTJ >= MAX(1,N), if STOR = 'F'.
C             LDJTJ >= 1,        if STOR = 'P'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             Currently, this array is not used.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 0.
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
C     The matrix product is computed columnn-wise, exploiting the
C     symmetry. BLAS 3 routines DGEMM and DSYRK are used if STOR = 'F',
C     and BLAS 2 routine DGEMV is used if STOR = 'P'.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
C
C     REVISIONS
C
C     V. Sima, Dec. 2001, Mar. 2002.
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
      CHARACTER         STOR, UPLO
      INTEGER           INFO, LDJ, LDJTJ, LDPAR, LDWORK, LIPAR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DPAR(*), DWORK(*), J(LDJ,*), JTJ(*)
      INTEGER           IPAR(*)
C     .. Local Scalars ..
      LOGICAL           FULL, UPPER
      INTEGER           BN, BSM, BSN, I1, IBSM, IBSN, II, JL, K, M,
     $                  NBSN, NTHS, ST
      DOUBLE PRECISION  C
C     .. Local Arrays ..
      DOUBLE PRECISION  TMP(1)
      INTEGER           ITMP(1)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DLASET, DSYRK, NF01BV,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
C
      FULL  = LSAME( STOR, 'F' )
      UPPER = LSAME( UPLO, 'U' )
C
      IF( .NOT.( FULL .OR. LSAME( STOR, 'P' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSEIF ( N.LT.0 ) THEN
         INFO = -3
      ELSEIF ( LIPAR.LT.4 ) THEN
         INFO = -5
      ELSEIF ( LDPAR.LT.1 ) THEN
         INFO = -7
      ELSEIF ( LDJTJ.LT.1 .OR. ( FULL .AND. LDJTJ.LT.N ) ) THEN
         INFO = -11
      ELSEIF ( LDWORK.LT.0 ) THEN
         INFO = -13
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
            INFO = -4
         ELSEIF ( N.NE.NTHS + ST ) THEN
            INFO = -3
         ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
            INFO = -9
         END IF
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BU', -INFO )
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
      IF ( BN.LE.1 .OR. BSN.EQ.0 .OR. BSM.EQ.0 ) THEN
C
C        Special case, l <= 1 or BSN = 0 or BSM = 0: the Jacobian is
C        represented as a full matrix.
C
         ITMP(1) = M
         CALL NF01BV( STOR, UPLO, N, ITMP, 1, DPAR, 1, J, LDJ, JTJ,
     $                LDJTJ, DWORK, LDWORK, INFO )
         RETURN
      END IF
C
C     General case: l > 1, BSN > 0, BSM > 0.
C
      JL = BSN + 1
C
      IF ( FULL ) THEN
C
         NBSN = N*BSN
C
         IF ( UPPER ) THEN
C
C           Compute the leading upper triangular part (full storage).
C
            CALL DLASET( UPLO, BSN, BSN, ZERO, C, JTJ, LDJTJ )
            CALL DSYRK(  UPLO, 'Transpose', BSN, BSM, ONE, J, LDJ, ONE,
     $                   JTJ, LDJTJ )
            IBSN = BSN
            I1   = NBSN + 1
C
            DO 10 IBSM = BSM + 1, M, BSM
               II = I1 + IBSN
               CALL DLASET( 'Full', IBSN, BSN, ZERO, ZERO, JTJ(I1),
     $                      LDJTJ )
               I1 = I1 + NBSN
               CALL DLASET( UPLO, BSN, BSN, ZERO, C, JTJ(II), LDJTJ )
               CALL DSYRK(  UPLO, 'Transpose', BSN, BSM, ONE, J(IBSM,1),
     $                      LDJ, ONE, JTJ(II), LDJTJ )
               IBSN = IBSN + BSN
   10       CONTINUE
C
            IF ( ST.GT.0 ) THEN
C
C              Compute the last block column.
C
               DO 20 IBSM = 1, M, BSM
                  CALL DGEMM( 'Transpose', 'NoTranspose', BSN, ST, BSM,
     $                        ONE, J(IBSM,1), LDJ, J(IBSM,JL), LDJ,
     $                        ZERO, JTJ(I1), LDJTJ )
                  I1 = I1 + BSN
   20          CONTINUE
C
               CALL DLASET( UPLO, ST, ST, ZERO, C, JTJ(I1), LDJTJ )
               CALL DSYRK(  UPLO, 'Transpose', ST, M, ONE, J(1,JL),
     $                      LDJ, ONE, JTJ(I1), LDJTJ )
            END IF
C
         ELSE
C
C           Compute the leading lower triangular part (full storage).
C
            IBSN = NTHS
            II   = 1
C
            DO 30 IBSM = 1, M, BSM
               I1 = II + BSN
               CALL DLASET( UPLO, BSN, BSN, ZERO, C, JTJ(II), LDJTJ )
               CALL DSYRK(  UPLO, 'Transpose', BSN, BSM, ONE, J(IBSM,1),
     $                      LDJ, ONE, JTJ(II), LDJTJ )
               IBSN = IBSN - BSN
               CALL DLASET( 'Full', IBSN, BSN, ZERO, ZERO, JTJ(I1),
     $                      LDJTJ )
               II = I1 + NBSN
               IF ( ST.GT.0 )
     $            CALL DGEMM( 'Transpose', 'NoTranspose', ST, BSN, BSM,
     $                        ONE, J(IBSM,JL), LDJ, J(IBSM,1), LDJ,
     $                        ZERO, JTJ(I1+IBSN), LDJTJ )
   30       CONTINUE
C
            IF ( ST.GT.0 ) THEN
C
C              Compute the last diagonal block.
C
               CALL DLASET( UPLO, ST, ST, ZERO, C, JTJ(II), LDJTJ )
               CALL DSYRK(  UPLO, 'Transpose', ST, M, ONE, J(1,JL),
     $                      LDJ, ONE, JTJ(II), LDJTJ )
            END IF
C
         END IF
C
      ELSE
C
         TMP(1) = ZERO
C
         IF ( UPPER ) THEN
C
C           Compute the leading upper triangular part (packed storage).
C
            IBSN = 0
            I1   = 1
C
            DO 50 IBSM = 1, M, BSM
C
               DO 40 K = 1, BSN
                  II = I1 + IBSN
                  CALL DCOPY( IBSN, TMP, 0, JTJ(I1), 1 )
                  CALL DGEMV( 'Transpose', BSM, K, ONE, J(IBSM,1), LDJ,
     $                        J(IBSM,K), 1, ZERO, JTJ(II), 1 )
                  I1 = II + K
                  JTJ(I1-1) = JTJ(I1-1) + C
   40          CONTINUE
C
               IBSN = IBSN + BSN
   50       CONTINUE
C
C           Compute the last block column.
C
            DO 70 K = 1, ST
C
               DO 60 IBSM = 1, M, BSM
                  CALL DGEMV( 'Transpose', BSM, BSN, ONE, J(IBSM,1),
     $                        LDJ, J(IBSM,BSN+K), 1, ZERO, JTJ(I1), 1 )
                  I1 = I1 + BSN
   60          CONTINUE
C
               CALL DGEMV( 'Transpose', M, K, ONE, J(1,JL), LDJ,
     $                     J(1,BSN+K), 1, ZERO, JTJ(I1), 1 )
               I1 = I1 + K
               JTJ(I1-1) = JTJ(I1-1) + C
   70       CONTINUE
C
         ELSE
C
C           Compute the leading lower triangular part (packed storage).
C
            IBSN = NTHS
            II   = 1
C
            DO 90 IBSM = 1, M, BSM
               IBSN = IBSN - BSN
C
               DO 80 K = 1, BSN
                  I1 = II + BSN - K + 1
                  CALL DCOPY( IBSN, TMP, 0, JTJ(I1), 1 )
                  CALL DGEMV( 'Transpose', BSM, BSN-K+1, ONE, J(IBSM,K),
     $                        LDJ, J(IBSM,K), 1, ZERO, JTJ(II), 1 )
                  JTJ(II) = JTJ(II) + C
                  I1 = I1 + IBSN
                  II = I1 + ST
                  IF ( ST.GT.0 )
     $               CALL DGEMV( 'Transpose', BSM, ST, ONE, J(IBSM,JL),
     $                           LDJ, J(IBSM,K), 1, ZERO, JTJ(I1), 1 )
   80          CONTINUE
C
   90       CONTINUE
C
C           Compute the last diagonal block.
C
            DO 100 K = 1, ST
               CALL DGEMV( 'Transpose', M, ST-K+1, ONE, J(1,BSN+K), LDJ,
     $                     J(1,BSN+K), 1, ZERO, JTJ(II), 1 )
               JTJ(II) = JTJ(II) + C
               II = II + ST - K + 1
  100       CONTINUE
C
         END IF
C
      END IF
C
      RETURN
C
C *** Last line of NF01BU ***
      END
