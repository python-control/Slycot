      SUBROUTINE MB04PA( LHAM, N, K, NB, A, LDA, QG, LDQG, XA, LDXA,
     $                   XG, LDXG, XQ, LDXQ, YA, LDYA, CS, TAU, DWORK )
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
C     To reduce a Hamiltonian like matrix
C
C                   [  A   G  ]           T          T
C              H =  [       T ] ,    G = G ,    Q = Q,
C                   [  Q  -A  ]
C
C     or a skew-Hamiltonian like matrix
C
C                   [  A   G  ]            T          T
C              W =  [       T ] ,    G = -G ,   Q = -Q,
C                   [  Q   A  ]
C
C     so that elements below the (k+1)-th subdiagonal in the first nb
C     columns of the (k+n)-by-n matrix A, and offdiagonal elements
C     in the first nb columns and rows of the n-by-n matrix Q are zero.
C
C     The reduction is performed by an orthogonal symplectic
C     transformation UU'*H*UU and matrices U, XA, XG, XQ, and YA are
C     returned so that
C
C                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'+ XG*U' ]
C         UU'*H*UU = [                                           ].
C                    [ Qout + U*XQ'+ XQ*U'  -Aout'- XA*U'- U*YA' ]
C
C     Similarly,
C
C                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'- XG*U' ]
C         UU'*W*UU = [                                           ].
C                    [ Qout + U*XQ'- XQ*U'   Aout'+ XA*U'+ U*YA' ]
C
C     This is an auxiliary routine called by MB04PB.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LHAM    LOGICAL
C             Specifies the type of matrix to be reduced:
C             = .FALSE. :  skew-Hamiltonian like W;
C             = .TRUE.  :  Hamiltonian like H.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     K       (input) INTEGER
C             The offset of the reduction. Elements below the (K+1)-th
C             subdiagonal in the first NB columns of A are reduced
C             to zero.  K >= 0.
C
C     NB      (input) INTEGER
C             The number of columns/rows to be reduced.  N > NB >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading (K+N)-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading (K+N)-by-N part of this array
C             contains the matrix Aout and in the zero part
C             information about the elementary reflectors used to
C             compute the reduction.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,K+N).
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C                            (LDQG,N+1)
C             On entry, the leading N+K-by-N+1 part of this array must
C             contain in the bottom left part the lower triangular part
C             of the N-by-N matrix Q and in the remainder the upper
C             trapezoidal part of the last N columns of the N+K-by-N+K
C             matrix G.
C             On exit, the leading N+K-by-N+1 part of this array
C             contains parts of the matrices Q and G in the same fashion
C             as on entry only that the zero parts of Q contain
C             information about the elementary reflectors used to
C             compute the reduction. Note that if LHAM = .FALSE. then
C             the (K-1)-th and K-th subdiagonals are not referenced.
C
C     LDQG    INTEGER
C             The leading dimension of the array QG. LDQG >= MAX(1,N+K).
C
C     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB)
C             On exit, the leading N-by-(2*NB) part of this array
C             contains the matrix XA.
C
C     LDXA    INTEGER
C             The leading dimension of the array XA.  LDXA >= MAX(1,N).
C
C     XG      (output) DOUBLE PRECISION array, dimension (LDXG,2*NB)
C             On exit, the leading (K+N)-by-(2*NB) part of this array
C             contains the matrix XG.
C
C     LDXG    INTEGER
C             The leading dimension of the array XG. LDXG >= MAX(1,K+N).
C
C     XQ      (output) DOUBLE PRECISION array, dimension (LDXQ,2*NB)
C             On exit, the leading N-by-(2*NB) part of this array
C             contains the matrix XQ.
C
C     LDXQ    INTEGER
C             The leading dimension of the array XQ.  LDXQ >= MAX(1,N).
C
C     YA      (output) DOUBLE PRECISION array, dimension (LDYA,2*NB)
C             On exit, the leading (K+N)-by-(2*NB) part of this array
C             contains the matrix YA.
C
C     LDYA    INTEGER
C             The leading dimension of the array YA. LDYA >= MAX(1,K+N).
C
C     CS      (output) DOUBLE PRECISION array, dimension (2*NB)
C             On exit, the first 2*NB elements of this array contain the
C             cosines and sines of the symplectic Givens rotations used
C             to compute the reduction.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (NB)
C             On exit, the first NB elements of this array contain the
C             scalar factors of some of the elementary reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (3*NB)
C
C     METHOD
C
C     For details regarding the representation of the orthogonal
C     symplectic matrix UU within the arrays A, QG, CS, TAU see the
C     description of MB04PU.
C
C     The contents of A and QG on exit are illustrated by the following
C     example with n = 5, k = 2 and nb = 2:
C
C           ( a  r  r  a  a  )         ( g  g  r  r  g  g  )
C           ( a  r  r  a  a  )         ( g  g  r  r  g  g  )
C           ( a  r  r  a  a  )         ( q  g  r  r  g  g  )
C       A = ( r  r  r  r  r  ),   QG = ( t  r  r  r  r  r  ),
C           ( u2 r  r  r  r  )         ( u1 t  r  r  r  r  )
C           ( u2 u2 r  a  a  )         ( u1 u1 r  q  g  g  )
C           ( u2 u2 r  a  a  )         ( u1 u1 r  q  q  g  )
C
C     where a, g and q denote elements of the original matrices, r
C     denotes a modified element, t denotes a scalar factor of an
C     applied elementary reflector and ui denote elements of the
C     matrix U.
C
C     REFERENCES
C
C     [1] C. F. VAN LOAN:
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications, 61, pp. 233-251, 1984.
C
C     [2] D. KRESSNER:
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     REVISIONS
C
C     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DLAPVL).
C
C     KEYWORDS
C
C     Elementary matrix operations, Hamiltonian matrix,
C     skew-Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D+0 )
C     .. Scalar Arguments ..
      LOGICAL           LHAM
      INTEGER           K, LDA, LDQG, LDXA, LDXG, LDXQ, LDYA, N, NB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), CS(*), DWORK(*), QG(LDQG,*), TAU(*),
     $                  XA(LDXA,*), XG(LDXG,*), XQ(LDXQ,*), YA(LDYA,*)
C     .. Local Scalars ..
      INTEGER           I, J, NB1, NB2
      DOUBLE PRECISION  AKI, ALPHA, C, S, TAUQ, TEMP, TTEMP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMV, DLARFG, DLARTG, DROT, DSCAL,
     $                  DSYMV, MB01MD
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF ( N+K.LE.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      NB1 = NB + 1
      NB2 = NB + NB1
C
      IF ( LHAM ) THEN
         DO 50  I = 1, NB
C
C           Transform i-th columns of A and Q. See routine MB04PU.
C
            ALPHA = QG(K+I+1,I)
            CALL DLARFG( N-I, ALPHA, QG(K+MIN( I+2, N ),I), 1, TAUQ )
            QG(K+I+1,I) = ONE
            TEMP = -TAUQ*DDOT( N-I, QG(K+I+1,I), 1, A(K+I+1,I), 1 )
            CALL DAXPY( N-I, TEMP,  QG(K+I+1,I), 1, A(K+I+1,I), 1 )
            AKI = A(K+I+1,I)
            CALL DLARTG( AKI, ALPHA, C, S, A(K+I+1,I) )
            AKI = A(K+I+1,I)
            CALL DLARFG( N-I, AKI, A(K+MIN( I+2, N ),I), 1, TAU(I) )
            A(K+I+1,I) = ONE
C
C           Update XA with first Householder reflection.
C
C           xa = H(1:n,1:n)'*u1
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, A(K+I+1,I+1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, XA(I+1,I), 1 )
C           w1 = U1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, QG(K+I+1,1), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, DWORK, 1 )
C           xa = xa + XA1*w1
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
C           w2 = U2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB1), 1 )
C           xa = xa + XA2*w2
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB1), 1, ONE, XA(I+1,I), 1 )
C           temp = YA1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  QG(K+I+1,I), 1, ZERO, XA(1,I), 1 )
C           xa = xa + U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, XA(1,I), 1, ONE, XA(I+1,I), 1 )
C           temp = YA2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  QG(K+I+1,I), 1, ZERO, XA(1,I), 1 )
C           xa = xa + U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
C           xa = -tauq*xa
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update YA with first Householder reflection.
C
C           ya = H(1:n,1:n)*u1
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, YA(1,I), 1 )
C           temp = XA1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           ya = ya + U1*temp
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, DWORK(NB2), 1, ONE, YA(K+I+1,I), 1 )
C           temp = XA2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XA(I+1,NB1), LDXA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           ya = ya + U2*temp
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  DWORK(NB2), 1, ONE, YA(K+I+1,I), 1 )
C           ya = ya + YA1*w1
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
C           ya = ya + YA2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB1), 1, ONE, YA(1,I), 1 )
C           ya = -tauq*ya
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C           temp = -tauq*ya'*u1
            TEMP = -TAUQ*DDOT( N-I, QG(K+I+1,I), 1, YA(K+I+1,I), 1 )
C           ya = ya + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, YA(K+I+1,I), 1 )
C
C           Update (i+1)-th column of A.
C
C           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, QG(K+I+1,1), LDQG,
     $                  XA(I+1,1), LDXA, ONE, A(K+I+1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(K+I+1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', N+K, I, ONE, YA, LDYA,
     $                  QG(K+I+1,1), LDQG, ONE, A(1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', N+K, I-1, ONE, YA(1,NB1), LDYA,
     $                  A(K+I+1,1), LDA, ONE, A(1,I+1), 1 )
C
C           Update (i+1)-th row of A.
C
            IF ( N.GT.I+1 ) THEN
C              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XA(I+2,1),
     $                     LDXA, QG(K+I+1,1), LDQG, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, A(K+I+1,1), LDA, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, YA(K+I+1,1), LDYA, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, YA(K+I+1,NB1), LDYA, ONE, A(K+I+1,I+2),
     $                     LDA )
            END IF
C
C           Annihilate updated parts in YA.
C
            DO 10  J = 1, I
               YA(K+I+1,J) = ZERO
   10       CONTINUE
            DO 20  J = 1, I-1
               YA(K+I+1,NB+J) = ZERO
   20       CONTINUE
C
C           Update XQ with first Householder reflection.
C
C           xq = Q*u1
            CALL DSYMV( 'Lower', N-I, ONE, QG(K+I+1,I+1), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XQ(I+1,I), 1 )
C           xq = xq + XQ1*w1
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
C           xq = xq + XQ2*w2
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB1), 1, ONE, XQ(I+1,I), 1 )
C           temp = XQ1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, XQ(1,I), 1 )
C           xq = xq + U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
C           temp = XQ2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XQ(I+1,NB1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, XQ(1,I), 1 )
C           xq = xq + U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
C           xq = -tauq*xq
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C           temp = -tauq/2*xq'*u1
            TEMP = -HALF*TAUQ*DDOT( N-I, QG(K+I+1,I), 1, XQ(I+1,I), 1 )
C           xq = xq + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, XQ(I+1,I), 1 )
C
C           Update (i+1)-th column and row of Q.
C
C           Q(:,i+1) = Q(:,i+1) + U1 * XQ1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, QG(K+I+1,1), LDQG,
     $                  XQ(I+1,1), LDXQ, ONE, QG(K+I+1,I+1), 1 )
C           Q(:,i+1) = Q(:,i+1) + U2 * XQ2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XQ(I+1,NB1), LDXQ, ONE, QG(K+I+1,I+1), 1 )
C           Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  QG(K+I+1,1), LDQG, ONE, QG(K+I+1,I+1), 1 )
C           Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, A(K+I+1,1), LDA, ONE, QG(K+I+1,I+1), 1 )
C
C           Update XG with first Householder reflection.
C
C           xg = G*u1
            CALL DGEMV( 'No transpose', K+I, N-I, ONE, QG(1,I+2), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XG(1,I), 1 )
            CALL DSYMV( 'Upper', N-I, ONE, QG(K+I+1,I+2), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XG(K+I+1,I), 1 )
C           xg = xg + XG1*w1
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
C           xg = xg + XG2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB1), 1, ONE, XG(1,I), 1 )
C           temp = XG1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XG(K+I+1,1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           xg = xg + U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, DWORK(NB2), 1, ONE, XG(K+I+1,I), 1 )
C           temp = XG2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XG(K+I+1,NB1),
     $                  LDXQ, QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           xg = xg + U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  DWORK(NB2), 1, ONE, XG(K+I+1,I), 1 )
C           xg = -tauq*xg
            CALL DSCAL( N+K, -TAUQ, XG(1,I), 1 )
C           temp = -tauq/2*xq'*u1
            TEMP = -HALF*TAUQ*DDOT( N-I, QG(K+I+1,I), 1, XG(K+I+1,I),
     $                              1 )
C           xg = xg + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, XG(K+I+1,I), 1 )
C
C           Update (i+1)-th column and row of G.
C
C           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', K+I, I, ONE, XG, LDXG,
     $                  QG(K+I+1,1), LDQG, ONE, QG(1,I+2), 1 )
C           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', K+I, I-1, ONE, XG(1,NB1), LDXG,
     $                  A(K+I+1,1), LDA, ONE, QG(1,I+2), 1 )
C           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, XG(K+I+1,1), LDXG,
     $                  QG(K+I+1,1), LDQG, ONE, QG(K+I+1,I+2), LDQG )
C           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XG(K+I+1,NB1),
     $                  LDXG, A(K+I+1,1), LDA, ONE, QG(K+I+1,I+2),
     $                  LDQG )
C           G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, QG(K+I+1,1), LDQG,
     $                  XG(K+I+1,1), LDXG, ONE, QG(K+I+1,I+2), LDQG )
C           G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XG(K+I+1,NB1), LDXG, ONE, QG(K+I+1,I+2), LDQG )
C
C           Annihilate updated parts in XG.
C
            DO 30  J = 1, I
               XG(K+I+1,J) = ZERO
   30       CONTINUE
            DO 40  J = 1, I-1
               XG(K+I+1,NB+J) = ZERO
   40       CONTINUE
C
C           Apply orthogonal symplectic Givens rotation.
C
            CALL DROT( K+I, A(1,I+1), 1, QG(1,I+2), 1, C, S )
            IF ( N.GT.I+1 ) THEN
               CALL DROT( N-I-1, A(K+I+2,I+1), 1, QG(K+I+1,I+3), LDQG,
     $                    C, S )
               CALL DROT( N-I-1, A(K+I+1,I+2), LDA, QG(K+I+2,I+1), 1, C,
     $                    S )
            END IF
            TEMP  = A(K+I+1,I+1)
            TTEMP = QG(K+I+1,I+2)
            A(K+I+1,I+1)  = C*TEMP  + S*QG(K+I+1,I+1)
            QG(K+I+1,I+2) = C*TTEMP - S*TEMP
            QG(K+I+1,I+1) = -S*TEMP + C*QG(K+I+1,I+1)
            TTEMP = -S*TTEMP - C*TEMP
            TEMP  = A(K+I+1,I+1)
            QG(K+I+1,I+1) =  C*QG(K+I+1,I+1) + S*TTEMP
            A(K+I+1,I+1)  =  C*TEMP + S*QG(K+I+1,I+2)
            QG(K+I+1,I+2) = -S*TEMP + C*QG(K+I+1,I+2)
            CS(2*I-1) = C
            CS(2*I)   = S
            QG(K+I+1,I) = TAUQ
C
C           Update XA with second Householder reflection.
C
C           xa = H(1:n,1:n)'*u2
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, A(K+I+1,I+1), LDA,
     $                  A(K+I+1,I), 1, ZERO, XA(I+1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              w1 = U1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                    LDQG, A(K+I+2,I), 1, ZERO, DWORK, 1 )
C              xa = xa + XA1*w1
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XA(I+2,1),
     $                     LDXA, DWORK, 1, ONE, XA(I+2,NB+I), 1 )
C              w2 = U2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, A(K+I+2,I), 1, ZERO, DWORK(NB1), 1 )
C              xa = xa + XA2*w2
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, DWORK(NB1), 1, ONE, XA(I+2,NB+I), 1 )
C              temp = YA1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, YA(K+I+2,1),
     $                     LDYA, A(K+I+2,I), 1, ZERO, XA(1,NB+I), 1 )
C              xa = xa + U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, XA(1,NB+I), 1, ONE, XA(I+2,NB+I), 1 )
C              temp = YA2'*u1
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, YA(K+I+2,NB1),
     $                     LDYA, A(K+I+2,I), 1, ZERO, XA(1,NB+I), 1 )
C              xa = xa + U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, XA(1,NB+I), 1, ONE, XA(I+2,NB+I), 1 )
            END IF
C           xa = -tau*xa
            CALL DSCAL( N-I, -TAU(I), XA(I+1,NB+I), 1 )
C
C           Update YA with second Householder reflection.
C
C           ya = H(1:n,1:n)*u2
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  A(K+I+1,I), 1, ZERO, YA(1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              temp = XA1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                     A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              ya = ya + U1*temp
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, DWORK(NB2), 1, ONE, YA(K+I+2,NB+I), 1 )
C              temp = XA2'*u1
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              ya = ya + U2*temp
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, DWORK(NB2), 1, ONE, YA(K+I+2,NB+I), 1 )
            END IF
C           ya = ya + YA1*w1
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                     DWORK, 1, ONE, YA(1,NB+I), 1 )
C           ya = ya + YA2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB1), 1, ONE, YA(1,NB+I), 1 )
C           ya = -tau*ya
            CALL DSCAL( K+N, -TAU(I), YA(1,NB+I), 1 )
C           temp = -tau*ya'*u2
            TEMP = -TAU(I)*DDOT( N-I, A(K+I+1,I), 1, YA(K+I+1,NB+I), 1 )
C           ya = ya + temp*u2
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, YA(K+I+1,NB+I), 1 )
C
C           Update (i+1)-th column of A.
C
C           H(1:n,i+1) = H(1:n,i+1) + ya
            CALL DAXPY( K+N, ONE, YA(1,NB+I), 1, A(1,I+1), 1 )
C           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2
            CALL DAXPY( N-I, XA(I+1,NB+I), A(K+I+1,I), 1, A(K+I+1,I+1),
     $                  1 )
C
C           Update (i+1)-th row of A.
C
            IF ( N.GT.I+1 ) THEN
C              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)';
               CALL DAXPY( N-I-1, ONE, XA(I+2,NB+I), 1, A(K+I+1,I+2),
     $                     LDA )
C              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)'
               CALL DAXPY( N-I-1, YA(K+I+1,NB+I), A(K+I+2,I), 1,
     $                     A(K+I+1,I+2), LDA )
            END IF
C
C           Annihilate updated parts in YA.
C
            YA(K+I+1,NB+I) = ZERO
C
C           Update XQ with second Householder reflection.
C
C           xq = Q*u2
            CALL DSYMV( 'Lower', N-I, ONE, QG(K+I+1,I+1), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XQ(I+1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              xq = xq + XQ1*w1
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XQ(I+2,1),
     $                     LDXQ, DWORK, 1, ONE, XQ(I+2,NB+I), 1 )
C              xq = xq + XQ2*w2
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XQ(I+2,NB1),
     $                     LDXQ, DWORK(NB1), 1, ONE, XQ(I+2,NB+I), 1 )
C              temp = XQ1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                     A(K+I+2,I), 1, ZERO, XQ(1,NB+I), 1 )
C              xq = xq + U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, XQ(1,NB+I), 1, ONE, XQ(I+2,NB+I), 1 )
C              temp = XQ2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XQ(I+2,NB1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, XQ(1,NB+I), 1 )
C              xq = xq + U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, XQ(1,NB+I), 1, ONE, XQ(I+2,NB+I), 1 )
            END IF
C           xq = -tauq*xq
            CALL DSCAL( N-I, -TAU(I), XQ(I+1,NB+I), 1 )
C           temp = -tauq/2*xq'*u2
            TEMP = -HALF*TAU(I)*DDOT( N-I, A(K+I+1,I), 1, XQ(I+1,NB+I),
     $                              1 )
C           xq = xq + temp*u2
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, XQ(I+1,NB+I), 1 )
C
C           Update (i+1)-th column and row of Q.
C
            CALL DAXPY( N-I, ONE, XQ(I+1,NB+I), 1, QG(K+I+1,I+1), 1 )
C           H(1:n,n+i+1) = H(1:n,n+i+1) + U * XQ(i+1,:)';
            CALL DAXPY( N-I, XQ(I+1,NB+I), A(K+I+1,I), 1,
     $                  QG(K+I+1,I+1), 1 )
C
C           Update XG with second Householder reflection.
C
C           xg = G*u2
            CALL DGEMV( 'No transpose', K+I, N-I, ONE, QG(1,I+2), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XG(1,NB+I), 1 )
            CALL DSYMV( 'Upper', N-I, ONE, QG(K+I+1,I+2), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XG(K+I+1,NB+I), 1 )
C           xg = xg + XG1*w1
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,NB+I), 1 )
C           xg = xg + XG2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB1), 1, ONE, XG(1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              temp = XG1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              xg = xg + U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, DWORK(NB2), 1, ONE, XG(K+I+2,NB+I), 1 )
C              temp = XG2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XG(K+I+2,NB1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              xg = xg + U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, DWORK(NB2), 1, ONE, XG(K+I+2,NB+I), 1 )
            END IF
C           xg = -tauq*xg
            CALL DSCAL( N+K, -TAU(I), XG(1,NB+I), 1 )
C           temp = -tauq/2*xg'*u1
            TEMP = -HALF*TAU(I)*DDOT( N-I, A(K+I+1,I), 1,
     $                                XG(K+I+1,NB+I), 1 )
C           xg = xg + temp*u1
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, XG(K+I+1,NB+I), 1 )
C
C           Update (i+1)-th column and row of G.
C
            CALL DAXPY( K+I, ONE, XG(1,NB+I), 1, QG(1,I+2), 1 )
            CALL DAXPY( N-I, ONE, XG(K+I+1,NB+I), 1, QG(K+I+1,I+2),
     $                  LDQG )
            CALL DAXPY( N-I, XG(K+I+1,NB+I), A(K+I+1,I), 1,
     $                  QG(K+I+1,I+2), LDQG )
C
C           Annihilate updated parts in XG.
C
            XG(K+I+1,NB+I) = ZERO
C
            A(K+I+1,I) = AKI
   50    CONTINUE
      ELSE
         DO 100  I = 1, NB
C
C           Transform i-th columns of A and Q.
C
            ALPHA = QG(K+I+1,I)
            CALL DLARFG( N-I, ALPHA, QG(K+MIN( I+2, N ),I), 1, TAUQ )
            QG(K+I+1,I) = ONE
            TEMP = -TAUQ*DDOT( N-I, QG(K+I+1,I), 1, A(K+I+1,I), 1 )
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, A(K+I+1,I), 1 )
            AKI = A(K+I+1,I)
            CALL DLARTG( AKI, ALPHA, C, S, A(K+I+1,I) )
            AKI = A(K+I+1,I)
            CALL DLARFG( N-I, AKI, A(K+MIN( I+2, N ),I), 1, TAU(I) )
            A(K+I+1,I) = ONE
C
C           Update XA with first Householder reflection.
C
C           xa = H(1:n,1:n)'*u1
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, A(K+I+1,I+1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, XA(I+1,I), 1 )
C           w1 = U1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, QG(K+I+1,1), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, DWORK, 1 )
C           xa = xa + XA1*w1
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
C           w2 = U2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB1), 1 )
C           xa = xa + XA2*w2
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB1), 1, ONE, XA(I+1,I), 1 )
C           temp = YA1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  QG(K+I+1,I), 1, ZERO, XA(1,I), 1 )
C           xa = xa + U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, XA(1,I), 1, ONE, XA(I+1,I), 1 )
C           temp = YA2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  QG(K+I+1,I), 1, ZERO, XA(1,I), 1 )
C           xa = xa + U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
C           xa = -tauq*xa
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update YA with first Householder reflection.
C
C           ya = H(1:n,1:n)*u1
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  QG(K+I+1,I), 1, ZERO, YA(1,I), 1 )
C           temp = XA1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           ya = ya + U1*temp
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, QG(K+I+1,1),
     $                  LDQG, DWORK(NB2), 1, ONE, YA(K+I+1,I), 1 )
C           temp = XA2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XA(I+1,NB1), LDXA,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           ya = ya + U2*temp
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  DWORK(NB2), 1, ONE, YA(K+I+1,I), 1 )
C           ya = ya + YA1*w1
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
C           ya = ya + YA2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB1), 1, ONE, YA(1,I), 1 )
C           ya = -tauq*ya
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C           temp = -tauq*ya'*u1
            TEMP = -TAUQ*DDOT( N-I, QG(K+I+1,I), 1, YA(K+I+1,I), 1 )
C           ya = ya + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, YA(K+I+1,I), 1 )
C
C           Update (i+1)-th column of A.
C
C           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I, ONE, QG(K+I+1,1), LDQG,
     $                  XA(I+1,1), LDXA, ONE, A(K+I+1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)';
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(K+I+1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', N+K, I, ONE, YA, LDYA,
     $                  QG(K+I+1,1), LDQG, ONE, A(1,I+1), 1 )
C           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', N+K, I-1, ONE, YA(1,NB1), LDYA,
     $                  A(K+I+1,1), LDA, ONE, A(1,I+1), 1 )
C
C           Update (i+1)-th row of A.
C
            IF ( N.GT.I+1 ) THEN
C              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XA(I+2,1),
     $                     LDXA, QG(K+I+1,1), LDQG, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, A(K+I+1,1), LDA, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, YA(K+I+1,1), LDYA, ONE, A(K+I+1,I+2),
     $                     LDA )
C              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)'
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, YA(K+I+1,NB1), LDYA, ONE, A(K+I+1,I+2),
     $                     LDA )
            END IF
C
C           Annihilate updated parts in YA.
C
            DO 60  J = 1, I
               YA(K+I+1,J) = ZERO
   60       CONTINUE
            DO 70  J = 1, I-1
               YA(K+I+1,NB+J) = ZERO
   70       CONTINUE
C
C           Update XQ with first Householder reflection.
C
C           xq = Q*u1
            CALL MB01MD( 'Lower', N-I, ONE, QG(K+I+1,I+1), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XQ(I+1,I), 1 )
C           xq = xq + XQ1*w1
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
C           xq = xq + XQ2*w2
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB1), 1, ONE, XQ(I+1,I), 1 )
C           temp = XQ1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, XQ(1,I), 1 )
C           xq = xq - U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, -ONE, QG(K+I+1,1),
     $                  LDQG, XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
C           temp = XQ2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XQ(I+1,NB1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, XQ(1,I), 1 )
C           xq = xq - U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, -ONE, A(K+I+1,1), LDA,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
C           xq = -tauq*xq
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C           temp = -tauq/2*xq'*u1
            TEMP = -HALF*TAUQ*DDOT( N-I, QG(K+I+1,I), 1, XQ(I+1,I), 1 )
C           xq = xq + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, XQ(I+1,I), 1 )
C
C           Update (i+1)-th column and row of Q.
C
            IF ( N.GT.I+1 ) THEN
C              Q(:,i+1) = Q(:,i+1) - U1 * XQ1(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I, -ONE, QG(K+I+2,1),
     $                     LDQG, XQ(I+1,1), LDXQ, ONE, QG(K+I+2,I+1),
     $                     1 )
C              Q(:,i+1) = Q(:,i+1) - U2 * XQ2(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I-1, -ONE, A(K+I+2,1),
     $                     LDA, XQ(I+1,NB1), LDXQ, ONE, QG(K+I+2,I+1),
     $                     1 )
C              Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XQ(I+2,1),
     $                     LDXQ, QG(K+I+1,1), LDQG, ONE, QG(K+I+2,I+1),
     $                     1 )
C              Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XQ(I+2,NB1),
     $                     LDXQ, A(K+I+1,1), LDA, ONE, QG(K+I+2,I+1),
     $                     1 )
            END IF
C
C           Update XG with first Householder reflection.
C
C           xg = G*u1
            CALL DGEMV( 'No transpose', K+I, N-I, ONE, QG(1,I+2), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XG(1,I), 1 )
            CALL MB01MD( 'Upper', N-I, ONE, QG(K+I+1,I+2), LDQG,
     $                  QG(K+I+1,I), 1, ZERO, XG(K+I+1,I), 1 )
C           xg = xg + XG1*w1
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
C           xg = xg + XG2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB1), 1, ONE, XG(1,I), 1 )
C           temp = XG1'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XG(K+I+1,1), LDXQ,
     $                  QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           xg = xg - U1*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, -ONE, QG(K+I+1,1),
     $                  LDQG, DWORK(NB2), 1, ONE, XG(K+I+1,I), 1 )
C           temp = XG2'*u1
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, XG(K+I+1,NB1),
     $                  LDXQ, QG(K+I+1,I), 1, ZERO, DWORK(NB2), 1 )
C           xg = xg - U2*temp
            CALL DGEMV( 'No Transpose', N-I, I-1, -ONE, A(K+I+1,1), LDA,
     $                  DWORK(NB2), 1, ONE, XG(K+I+1,I), 1 )
C           xg = -tauq*xg
            CALL DSCAL( N+K, -TAUQ, XG(1,I), 1 )
C           temp = -tauq/2*xq'*u1
            TEMP = -HALF*TAUQ*DDOT( N-I, QG(K+I+1,I), 1, XG(K+I+1,I),
     $                              1 )
C           xg = xg + temp*u1
            CALL DAXPY( N-I, TEMP, QG(K+I+1,I), 1, XG(K+I+1,I), 1 )
C
C           Update (i+1)-th column and row of G.
C
C           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)';
            CALL DGEMV( 'No transpose', K+I, I, ONE, XG, LDXG,
     $                  QG(K+I+1,1), LDQG, ONE, QG(1,I+2), 1 )
C           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)';
            CALL DGEMV( 'No transpose', K+I, I-1, ONE, XG(1,NB1), LDXG,
     $                  A(K+I+1,1), LDA, ONE, QG(1,I+2), 1 )
            IF ( N.GT.I+1 ) THEN
C              G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I, -ONE, XG(K+I+2,1),
     $                     LDXG, QG(K+I+1,1), LDQG, ONE, QG(K+I+1,I+3),
     $                     LDQG )
C              G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I-1, -ONE,
     $                     XG(K+I+2,NB1), LDXG, A(K+I+1,1), LDA, ONE,
     $                     QG(K+I+1,I+3), LDQG )
C              G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, XG(K+I+1,1), LDXG, ONE, QG(K+I+1,I+3),
     $                     LDQG )
C              G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)';
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, XG(K+I+1,NB1), LDXG, ONE, QG(K+I+1,I+3),
     $                     LDQG )
            END IF
C
C           Annihilate updated parts in XG.
C
            DO 80  J = 1, I
               XG(K+I+1,J) = ZERO
   80       CONTINUE
            DO 90  J = 1, I-1
               XG(K+I+1,NB+J) = ZERO
   90       CONTINUE
C
C           Apply orthogonal symplectic Givens rotation.
C
            CALL DROT( K+I, A(1,I+1), 1, QG(1,I+2), 1, C, S )
            IF ( N.GT.I+1 ) THEN
               CALL DROT( N-I-1, A(K+I+2,I+1), 1, QG(K+I+1,I+3), LDQG,
     $                    C, -S )
               CALL DROT( N-I-1, A(K+I+1,I+2), LDA, QG(K+I+2,I+1), 1,
     $                    C, -S )
            END IF
            CS(2*I-1) = C
            CS(2*I)   = S
            QG(K+I+1,I) = TAUQ
C
C           Update XA with second Householder reflection.
C
C           xa = H(1:n,1:n)'*u2
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, A(K+I+1,I+1), LDA,
     $                  A(K+I+1,I), 1, ZERO, XA(I+1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              w1 = U1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                    LDQG, A(K+I+2,I), 1, ZERO, DWORK, 1 )
C              xa = xa + XA1*w1
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XA(I+2,1),
     $                     LDXA, DWORK, 1, ONE, XA(I+2,NB+I), 1 )
C              w2 = U2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, A(K+I+2,I), 1, ZERO, DWORK(NB1), 1 )
C              xa = xa + XA2*w2
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, DWORK(NB1), 1, ONE, XA(I+2,NB+I), 1 )
C              temp = YA1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, YA(K+I+2,1),
     $                     LDYA, A(K+I+2,I), 1, ZERO, XA(1,NB+I), 1 )
C              xa = xa + U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, XA(1,NB+I), 1, ONE, XA(I+2,NB+I), 1 )
C              temp = YA2'*u1
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, YA(K+I+2,NB1),
     $                     LDYA, A(K+I+2,I), 1, ZERO, XA(1,NB+I), 1 )
C              xa = xa + U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, XA(1,NB+I), 1, ONE, XA(I+2,NB+I), 1 )
            END IF
C           xa = -tau*xa
            CALL DSCAL( N-I, -TAU(I), XA(I+1,NB+I), 1 )
C
C           Update YA with second Householder reflection.
C
C           ya = H(1:n,1:n)*u2
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  A(K+I+1,I), 1, ZERO, YA(1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              temp = XA1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                     A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              ya = ya + U1*temp
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, QG(K+I+2,1),
     $                     LDQG, DWORK(NB2), 1, ONE, YA(K+I+2,NB+I), 1 )
C              temp = XA2'*u1
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XA(I+2,NB1),
     $                     LDXA, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              ya = ya + U2*temp
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, A(K+I+2,1),
     $                     LDA, DWORK(NB2), 1, ONE, YA(K+I+2,NB+I), 1 )
            END IF
C           ya = ya + YA1*w1
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                     DWORK, 1, ONE, YA(1,NB+I), 1 )
C           ya = ya + YA2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB1), 1, ONE, YA(1,NB+I), 1 )
C           ya = -tau*ya
            CALL DSCAL( K+N, -TAU(I), YA(1,NB+I), 1 )
C           temp = -tau*ya'*u2
            TEMP = -TAU(I)*DDOT( N-I, A(K+I+1,I), 1, YA(K+I+1,NB+I), 1 )
C           ya = ya + temp*u2
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, YA(K+I+1,NB+I), 1 )
C
C           Update (i+1)-th column of A.
C
C           H(1:n,i+1) = H(1:n,i+1) + ya
            CALL DAXPY( K+N, ONE, YA(1,NB+I), 1, A(1,I+1), 1 )
C           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2
            CALL DAXPY( N-I, XA(I+1,NB+I), A(K+I+1,I), 1, A(K+I+1,I+1),
     $                  1 )
C
C           Update (i+1)-th row of A.
C
            IF ( N.GT.I+1 ) THEN
C              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)';
               CALL DAXPY( N-I-1, ONE, XA(I+2,NB+I), 1, A(K+I+1,I+2),
     $                     LDA )
C              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)'
               CALL DAXPY( N-I-1, YA(K+I+1,NB+I), A(K+I+2,I), 1,
     $                     A(K+I+1,I+2), LDA )
            END IF
C
C           Annihilate updated parts in YA.
C
            YA(K+I+1,NB+I) = ZERO
C
C           Update XQ with second Householder reflection.
C
C           xq = Q*u2
            CALL MB01MD( 'Lower', N-I, ONE, QG(K+I+1,I+1), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XQ(I+1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              xq = xq + XQ1*w1
               CALL DGEMV( 'No transpose', N-I-1, I, ONE, XQ(I+2,1),
     $                     LDXQ, DWORK, 1, ONE, XQ(I+2,NB+I), 1 )
C              xq = xq + XQ2*w2
               CALL DGEMV( 'No transpose', N-I-1, I-1, ONE, XQ(I+2,NB1),
     $                     LDXQ, DWORK(NB1), 1, ONE, XQ(I+2,NB+I), 1 )
C              temp = XQ1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                     A(K+I+2,I), 1, ZERO, XQ(1,NB+I), 1 )
C              xq = xq - U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, -ONE, QG(K+I+2,1),
     $                     LDQG, XQ(1,NB+I), 1, ONE, XQ(I+2,NB+I), 1 )
C              temp = XQ2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XQ(I+2,NB1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, XQ(1,NB+I), 1 )
C              xq = xq - U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, -ONE, A(K+I+2,1),
     $                     LDA, XQ(1,NB+I), 1, ONE, XQ(I+2,NB+I), 1 )
            END IF
C           xq = -tauq*xq
            CALL DSCAL( N-I, -TAU(I), XQ(I+1,NB+I), 1 )
C           temp = -tauq/2*xq'*u2
            TEMP = -HALF*TAU(I)*DDOT( N-I, A(K+I+1,I), 1, XQ(I+1,NB+I),
     $                              1 )
C           xq = xq + temp*u2
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, XQ(I+1,NB+I), 1 )
C
C           Update (i+1)-th column and row of Q.
C
            IF ( N.GT.I+1 ) THEN
               CALL DAXPY( N-I-1, ONE, XQ(I+2,NB+I), 1, QG(K+I+2,I+1),
     $                     1 )
C              H(1:n,n+i+1) = H(1:n,n+i+1) - U * XQ(i+1,:)';
               CALL DAXPY( N-I-1, -XQ(I+1,NB+I), A(K+I+2,I), 1,
     $                     QG(K+I+2,I+1), 1 )
            END IF
C
C           Update XG with second Householder reflection.
C
C           xg = G*u2
            CALL DGEMV( 'No transpose', K+I, N-I, ONE, QG(1,I+2), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XG(1,NB+I), 1 )
            CALL MB01MD( 'Upper', N-I, ONE, QG(K+I+1,I+2), LDQG,
     $                  A(K+I+1,I), 1, ZERO, XG(K+I+1,NB+I), 1 )
C           xg = xg + XG1*w1
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,NB+I), 1 )
C           xg = xg + XG2*w2
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB1), 1, ONE, XG(1,NB+I), 1 )
            IF ( N.GT.I+1 ) THEN
C              temp = XG1'*u2
               CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              xg = xg - U1*temp
               CALL DGEMV( 'No Transpose', N-I-1, I, -ONE, QG(K+I+2,1),
     $                     LDQG, DWORK(NB2), 1, ONE, XG(K+I+2,NB+I), 1 )
C              temp = XG2'*u2
               CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, XG(K+I+2,NB1),
     $                     LDXQ, A(K+I+2,I), 1, ZERO, DWORK(NB2), 1 )
C              xg = xg - U2*temp
               CALL DGEMV( 'No Transpose', N-I-1, I-1, -ONE, A(K+I+2,1),
     $                     LDA, DWORK(NB2), 1, ONE, XG(K+I+2,NB+I), 1 )
            END IF
C           xg = -tauq*xg
            CALL DSCAL( N+K, -TAU(I), XG(1,NB+I), 1 )
C           temp = -tauq/2*xg'*u1
            TEMP = -HALF*TAU(I)*DDOT( N-I, A(K+I+1,I), 1,
     $                                XG(K+I+1,NB+I), 1 )
C           xg = xg + temp*u1
            CALL DAXPY( N-I, TEMP, A(K+I+1,I), 1, XG(K+I+1,NB+I), 1 )
C
C           Update (i+1)-th column and row of G.
C
            CALL DAXPY( K+I, ONE, XG(1,NB+I), 1, QG(1,I+2), 1 )
            IF ( N.GT.I+1 ) THEN
               CALL DAXPY( N-I-1, -ONE, XG(K+I+2,NB+I), 1,
     $                     QG(K+I+1,I+3), LDQG )
               CALL DAXPY( N-I-1, XG(K+I+1,NB+I), A(K+I+2,I), 1,
     $                     QG(K+I+1,I+3), LDQG )
            END IF
C
C           Annihilate updated parts in XG.
C
            XG(K+I+1,NB+I) = ZERO
C
            A(K+I+1,I) = AKI
  100    CONTINUE
      END IF
C
      RETURN
C *** Last line of MB04PA ***
      END
