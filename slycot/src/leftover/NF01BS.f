      SUBROUTINE NF01BS( N, IPAR, LIPAR, FNORM, J, LDJ, E, JNORMS,
     $                   GNORM, IPVT, DWORK, LDWORK, INFO )
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
C     To compute the QR factorization of the Jacobian matrix J, as
C     received in compressed form from SLICOT Library routine NF01BD,
C
C            /  dy(1)/dwb(1)  |  dy(1)/ dtheta  \
C       Jc = |       :        |       :         | ,
C            \  dy(L)/dwb(L)  |  dy(L)/ dtheta  /
C
C     and to apply the transformation Q on the error vector e (in-situ).
C     The factorization is J*P = Q*R, where Q is a matrix with
C     orthogonal columns, P a permutation matrix, and R an upper
C     trapezoidal matrix with diagonal elements of nonincreasing
C     magnitude for each block column (see below). The 1-norm of the
C     scaled gradient is also returned.
C
C     Actually, the Jacobian J has the block form
C
C       dy(1)/dwb(1)       0         .....       0        dy(1)/dtheta
C            0        dy(2)/dwb(2)   .....       0        dy(2)/dtheta
C          .....         .....       .....     .....         .....
C            0           .....         0    dy(L)/dwb(L)  dy(L)/dtheta
C
C     but the zero blocks are omitted. The diagonal blocks have the
C     same size and correspond to the nonlinear part. The last block
C     column corresponds to the linear part. It is assumed that the
C     Jacobian matrix has at least as many rows as columns. The linear
C     or nonlinear parts can be empty. If L <= 1, the Jacobian is
C     represented as a full matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of columns of the Jacobian matrix J.
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
C                     BN*BSM >= N, if BN > 0;
C                     BSM >= N,    if BN = 0.
C             IPAR(4) must contain BSN, the number of columns of the
C                     blocks J_k, k = 1:BN.  BSN >= 0.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 4.
C
C     FNORM   (input) DOUBLE PRECISION
C             The Euclidean norm of the vector e.  FNORM >= 0.
C
C     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, NC)
C             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1.
C             On entry, the leading NR-by-NC part of this array must
C             contain the (compressed) representation (Jc) of the
C             Jacobian matrix J, where NR = BSM if BN <= 1, and
C             NR = BN*BSM, if BN > 1.
C             On exit, the leading N-by-NC part of this array contains
C             a (compressed) representation of the upper triangular
C             factor R of the Jacobian matrix. The matrix R has the same
C             structure as the Jacobian matrix J, but with an additional
C             diagonal block. Note that for efficiency of the later
C             calculations, the matrix R is delivered with the leading
C             dimension MAX(1,N), possibly much smaller than the value
C             of LDJ on entry.
C
C     LDJ     (input/output) INTEGER
C             The leading dimension of array J.
C             On entry, LDJ >= MAX(1,NR).
C             On exit,  LDJ >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (NR)
C             On entry, this array contains the vector e,
C             e = vec( Y - y ), where Y is set of output samples, and
C             vec denotes the concatenation of the columns of a matrix.
C             On exit, this array contains the updated vector Z*Q'*e,
C             where Z is the block row permutation matrix used in the
C             QR factorization of J (see METHOD).
C
C     JNORMS  (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the Euclidean norms of the columns
C             of the Jacobian matrix, considered in the initial order.
C
C     GNORM   (output) DOUBLE PRECISION
C             If FNORM > 0, the 1-norm of the scaled vector J'*e/FNORM,
C             with each element i further divided by JNORMS(i) (if
C             JNORMS(i) is nonzero).
C             If FNORM = 0, the returned value of GNORM is 0.
C
C     IPVT    (output) INTEGER array, dimension (N)
C             This array defines the permutation matrix P such that
C             J*P = Q*R. Column j of P is column IPVT(j) of the identity
C             matrix.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 1,      if N = 0 or BN <= 1 and BSM = N = 1;
C                               otherwise,
C             LDWORK >= 4*N+1,  if BN <= 1 or  BSN = 0;
C             LDWORK >= JWORK,  if BN >  1 and BSN > 0, where JWORK is
C                               given by the following procedure:
C              JWORK  = BSN + MAX(3*BSN+1,ST);
C              JWORK  = MAX(JWORK,4*ST+1),         if BSM > BSN;
C              JWORK  = MAX(JWORK,(BSM-BSN)*(BN-1)),
C                                                  if BSN < BSM < 2*BSN.
C             For optimum performance LDWORK should be larger.
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
C     A QR factorization with column pivoting of the matrix J is
C     computed, J*P = Q*R.
C
C     If l = L > 1, the R factor of the QR factorization has the same
C     structure as the Jacobian, but with an additional diagonal block.
C     Denote
C
C         /   J_1    0    ..   0   |  L_1  \
C         |    0    J_2   ..   0   |  L_2  |
C     J = |    :     :    ..   :   |   :   | .
C         |    :     :    ..   :   |   :   |
C         \    0     0    ..  J_l  |  L_l  /
C
C     The algorithm consists in two phases. In the first phase, the
C     algorithm uses QR factorizations with column pivoting for each
C     block J_k, k = 1:l, and applies the orthogonal matrix Q'_k to the
C     corresponding part of the last block column and of e. After all
C     block rows have been processed, the block rows are interchanged
C     so that the zeroed submatrices in the first l block columns are
C     moved to the bottom part. The same block row permutation Z is
C     also applied to the vector e. At the end of the first phase,
C     the structure of the processed matrix J is
C
C         /   R_1    0    ..   0   |  L^1_1  \
C         |    0    R_2   ..   0   |  L^1_2  |
C         |    :     :    ..   :   |    :    | .
C         |    :     :    ..   :   |    :    |
C         |    0     0    ..  R_l  |  L^1_l  |
C         |    0     0    ..   0   |  L^2_1  |
C         |    :     :    ..   :   |    :    |
C         \    0     0    ..   0   |  L^2_l  /
C
C     In the second phase, the submatrix L^2_1:l is triangularized
C     using an additional QR factorization with pivoting. (The columns
C     of L^1_1:l are also permuted accordingly.) Therefore, the column
C     pivoting is restricted to each such local block column.
C
C     If l <= 1, the matrix J is triangularized in one phase, by one
C     QR factorization with pivoting. In this case, the column
C     pivoting is global.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     Feb. 22, 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, Jacobian matrix, matrix algebra,
C     matrix operations, Wiener system.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDJ, LDWORK, LIPAR, N
      DOUBLE PRECISION  FNORM, GNORM
C     .. Array Arguments ..
      INTEGER           IPAR(*), IPVT(*)
      DOUBLE PRECISION  DWORK(*), E(*), J(*), JNORMS(*)
C     .. Local Scalars ..
      INTEGER           BN, BSM, BSN, I, IBSM, IBSN, IBSNI, ITAU, JL,
     $                  JLM, JWORK, K, L, M, MMN, NTHS, ST, WRKOPT
      DOUBLE PRECISION  SUM
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEQP3, DLACPY, DLAPMT, DORMQR, DSWAP,
     $                  MD03BX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSEIF( LIPAR.LT.4 ) THEN
         INFO = -3
      ELSEIF ( FNORM.LT.ZERO ) THEN
         INFO = -4
      ELSEIF ( LDJ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE
         ST   = IPAR(1)
         BN   = IPAR(2)
         BSM  = IPAR(3)
         BSN  = IPAR(4)
         NTHS = BN*BSN
         MMN  = BSM - BSN
         IF ( BN.GT.0 ) THEN
            M = BN*BSM
         ELSE
            M = N
         END IF
         IF ( MIN( ST, BN, BSM, BSN ).LT.0 ) THEN
            INFO = -2
         ELSEIF ( N.NE.NTHS + ST ) THEN
            INFO = -1
         ELSEIF ( M.LT.N ) THEN
            INFO = -2
         ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
            INFO = -6
         ELSE
            IF ( N.EQ.0 ) THEN
               JWORK = 1
            ELSEIF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
               IF ( BN.LE.1 .AND. BSM.EQ.1 .AND. N.EQ.1 ) THEN
                  JWORK = 1
               ELSE
                  JWORK = 4*N + 1
               END IF
            ELSE
               JWORK = BSN + MAX( 3*BSN + 1, ST )
               IF ( BSM.GT.BSN ) THEN
                  JWORK = MAX( JWORK, 4*ST + 1 )
                  IF ( BSM.LT.2*BSN )
     $               JWORK = MAX( JWORK, MMN*( BN - 1 ) )
               END IF
            END IF
            IF ( LDWORK.LT.JWORK )
     $         INFO = -12
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'NF01BS', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      GNORM = ZERO
      IF ( N.EQ.0 ) THEN
         LDJ = 1
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
C
C        Special case, l <= 1 or BSN = 0: the Jacobian is represented
C        as a full matrix.
C        (Note: Comments in the code beginning "Workspace:" describe the
C        minimal amount of real workspace needed at that point in the
C        code, as well as the preferred amount for good performance.
C        NB refers to the optimal block size for the immediately
C        following subroutine, as returned by ILAENV.)
C
C        Workspace: need:    4*N + 1;
C                   prefer:  3*N + ( N+1 )*NB.
C
         CALL MD03BX( M, N, FNORM, J, LDJ, E, JNORMS, GNORM, IPVT,
     $                DWORK, LDWORK, INFO )
         RETURN
      END IF
C
C     General case: l > 1 and BSN > 0.
C     Initialize the column pivoting indices.
C
      DO 10 I = 1, N
         IPVT(I) = 0
   10 CONTINUE
C
C     Compute the QR factorization with pivoting of J.
C     Pivoting is done separately on each block column of J.
C
      WRKOPT = 1
      IBSN   = 1
      JL     = LDJ*BSN + 1
      JWORK  = BSN + 1
C
      DO 30 IBSM = 1, M, BSM
C
C        Compute the QR factorization with pivoting of J_k, and apply Q'
C        to the corresponding part of the last block-column and of e.
C        Workspace: need:    4*BSN + 1;
C                   prefer:  3*BSN + ( BSN+1 )*NB.
C
         CALL DGEQP3( BSM, BSN, J(IBSM), LDJ, IPVT(IBSN), DWORK,
     $                DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         IF ( IBSM.GT.1 ) THEN
C
C           Adjust the column pivoting indices.
C
            DO 20 I = IBSN, IBSN + BSN - 1
               IPVT(I) = IPVT(I) + IBSN - 1
   20       CONTINUE
C
         END IF
C
         IF ( ST.GT.0 ) THEN
C
C           Workspace: need:    BSN + ST;
C                      prefer:  BSN + ST*NB.
C
            CALL DORMQR( 'Left', 'Transpose', BSM, ST, BSN, J(IBSM),
     $                   LDJ, DWORK, J(JL), LDJ, DWORK(JWORK),
     $                   LDWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
         END IF
C
C        Workspace: need:    BSN + 1;
C                   prefer:  BSN + NB.
C
         CALL DORMQR( 'Left', 'Transpose', BSM, 1, BSN, J(IBSM), LDJ,
     $                DWORK, E(IBSM), BSM, DWORK(JWORK), LDWORK-JWORK+1,
     $                INFO )
         JL   = JL   + BSM
         IBSN = IBSN + BSN
   30 CONTINUE
C
      IF ( MMN.GT.0 ) THEN
C
C        Case BSM > BSN.
C        Compute the original column norms for the first block column
C        of Jc.
C        Permute the rows of the first block column to move the zeroed
C        submatrices to the bottom. In the same loops, reshape the
C        first block column of R to have the leading dimension N.
C
         L = IPVT(1)
         JNORMS(L) = ABS( J(1) )
         IBSM = BSM + 1
         IBSN = BSN + 1
C
         DO 40 K = 1, BN - 1
            J(IBSN) = J(IBSM)
            L = IPVT(IBSN)
            JNORMS(L) = ABS( J(IBSN) )
            IBSM = IBSM + BSM
            IBSN = IBSN + BSN
   40    CONTINUE
C
         IBSN = IBSN + ST
C
         DO 60 I = 2, BSN
            IBSM = ( I - 1 )*LDJ + 1
            JL   = I
C
            DO 50 K = 1, BN
C
               DO 45 L = 0, I - 1
                  J(IBSN+L) = J(IBSM+L)
   45          CONTINUE
C
               L = IPVT(JL)
               JNORMS(L) = DNRM2( I, J(IBSN), 1 )
               IBSM = IBSM + BSM
               IBSN = IBSN + BSN
               JL   = JL   + BSN
   50       CONTINUE
C
            IBSN = IBSN + ST
   60    CONTINUE
C
C        Permute the rows of the second block column of Jc and of
C        the vector e.
C
         JL = LDJ*BSN
         IF ( BSM.GE.2*BSN ) THEN
C
C           A swap operation can be used.
C
            DO 80 I = 1, ST
               IBSN = BSN + 1
C
               DO 70 IBSM = BSM + 1, M, BSM
                  CALL DSWAP( MMN, J(JL+IBSM), 1, J(JL+IBSN), 1 )
                  IBSN = IBSN + BSN
   70          CONTINUE
C
               JL = JL + LDJ
   80       CONTINUE
C
C           Permute the rows of e.
C
            IBSN = BSN + 1
C
            DO 90 IBSM = BSM + 1, M, BSM
               CALL DSWAP( MMN, E(IBSM), 1, E(IBSN), 1 )
               IBSN = IBSN + BSN
   90       CONTINUE
C
         ELSE
C
C           A swap operation cannot be used.
C           Workspace: need:    ( BSM-BSN )*( BN-1 ).
C
            DO 110 I = 1, ST
               IBSN  = BSN + 1
               JLM   = JL  + IBSN
               JWORK = 1
C
               DO 100 IBSM = BSM + 1, M, BSM
                  CALL DCOPY( MMN, J(JLM), 1, DWORK(JWORK), 1 )
C
                  DO 105 K = JL, JL + BSN - 1
                     J(IBSN+K) = J(IBSM+K)
  105             CONTINUE
C
                  JLM   = JLM   + BSM
                  IBSN  = IBSN  + BSN
                  JWORK = JWORK + MMN
  100          CONTINUE
C
               CALL DCOPY( MMN*( BN-1 ), DWORK, 1, J(JL+IBSN), 1 )
               JL = JL + LDJ
  110       CONTINUE
C
C           Permute the rows of e.
C
            IBSN  = BSN + 1
            JLM   = IBSN
            JWORK = 1
C
            DO 120 IBSM = BSM + 1, M, BSM
               CALL DCOPY( MMN, E(JLM),  1, DWORK(JWORK), 1 )
C
               DO 115 K = 0, BSN - 1
                  E(IBSN+K) = E(IBSM+K)
  115          CONTINUE
C
               JLM   = JLM   + BSM
               IBSN  = IBSN  + BSN
               JWORK = JWORK + MMN
  120       CONTINUE
C
            CALL DCOPY( MMN*( BN-1 ), DWORK, 1, E(IBSN), 1 )
         END IF
C
         IF ( ST.GT.0 ) THEN
C
C           Compute the QR factorization with pivoting of the submatrix
C           L^2_1:l, and apply Q' to the corresponding part of e.
C
C           Workspace: need:    4*ST + 1;
C                      prefer:  3*ST + ( ST+1 )*NB.
C
            JL    = ( LDJ + BN )*BSN + 1
            ITAU  = 1
            JWORK = ITAU + ST
            CALL DGEQP3( MMN*BN, ST, J(JL), LDJ, IPVT(NTHS+1),
     $                   DWORK(ITAU), DWORK(JWORK), LDWORK-JWORK+1,
     $                   INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C           Permute columns of the upper part of the second block
C           column of Jc.
C
            CALL DLAPMT( .TRUE., NTHS, ST, J(JL-NTHS), LDJ,
     $                   IPVT(NTHS+1) )
C
C           Adjust the column pivoting indices.
C
            DO 130 I = NTHS + 1, N
               IPVT(I) = IPVT(I) + NTHS
  130       CONTINUE
C
C           Workspace: need:    ST + 1;
C                      prefer:  ST + NB.
C
            CALL DORMQR( 'Left', 'Transpose', MMN*BN, 1, ST, J(JL), LDJ,
     $                   DWORK(ITAU), E(IBSN), LDJ, DWORK(JWORK),
     $                   LDWORK-JWORK+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C           Reshape the second block column of R to have the leading
C           dimension N.
C
            IBSN = N*BSN + 1
            CALL DLACPY( 'Full', N, ST, J(LDJ*BSN+1), LDJ, J(IBSN), N )
C
C           Compute the original column norms for the second block
C           column.
C
            DO 140 I = NTHS + 1, N
               L = IPVT(I)
               JNORMS(L) = DNRM2( I, J(IBSN), 1 )
               IBSN = IBSN + N
  140       CONTINUE
C
         END IF
C
      ELSE
C
C        Case BSM = BSN.
C        Compute the original column norms for the first block column
C        of Jc.
C
         IBSN = 1
C
         DO 160 I = 1, BSN
            JL = I
C
            DO 150 K = 1, BN
               L = IPVT(JL)
               JNORMS(L) = DNRM2( I, J(IBSN), 1 )
               IBSN = IBSN + BSN
               JL   = JL   + BSN
  150       CONTINUE
C
            IBSN = IBSN + ST
  160    CONTINUE
C
         DO 170 I = NTHS + 1, N
            IPVT(I) = I
  170    CONTINUE
C
      END IF
C
C     Compute the norm of the scaled gradient.
C
      IF ( FNORM.NE.ZERO ) THEN
C
         DO 190 IBSN = 1, NTHS, BSN
            IBSNI = IBSN
C
            DO 180 I = 1, BSN
               L = IPVT(IBSN+I-1)
               IF ( JNORMS(L).NE.ZERO ) THEN
                  SUM   = DDOT( I, J(IBSNI), 1, E(IBSN), 1 )/FNORM
                  GNORM = MAX( GNORM, ABS( SUM/JNORMS(L) ) )
               END IF
               IBSNI = IBSNI + N
  180       CONTINUE
C
  190    CONTINUE
C
         IBSNI = N*BSN + 1
C
         DO 200 I = NTHS + 1, N
            L = IPVT(I)
            IF ( JNORMS(L).NE.ZERO ) THEN
               SUM   = DDOT( I, J(IBSNI), 1, E, 1 )/FNORM
               GNORM = MAX( GNORM, ABS( SUM/JNORMS(L) ) )
            END IF
            IBSNI = IBSNI + N
  200    CONTINUE
C
      END IF
C
      LDJ = N
      DWORK(1) = WRKOPT
      RETURN
C
C *** Last line of NF01BS ***
      END
