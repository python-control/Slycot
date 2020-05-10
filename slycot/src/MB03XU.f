      SUBROUTINE MB03XU( LTRA, LTRB, N, K, NB, A, LDA, B, LDB, G, LDG,
     $                   Q, LDQ, XA, LDXA, XB, LDXB, XG, LDXG, XQ, LDXQ,
     $                   YA, LDYA, YB, LDYB, YG, LDYG, YQ, LDYQ, CSL,
     $                   CSR, TAUL, TAUR, DWORK )
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
C     To reduce 2*nb columns and rows of a real (k+2n)-by-(k+2n)
C     matrix H:
C
C             [ op(A)   G   ]
C         H = [             ],
C             [  Q    op(B) ]
C
C     so that elements in the first nb columns below the k-th
C     subdiagonal of the (k+n)-by-n matrix op(A), in the first nb
C     columns and rows of the n-by-n matrix Q and in the first nb rows
C     above the diagonal of the n-by-(k+n) matrix op(B) are zero.
C     The reduction is performed by orthogonal symplectic
C     transformations UU'*H*VV and matrices U, V, YA, YB, YG, YQ, XA,
C     XB, XG, and XQ are returned so that
C
C                    [ op(Aout)+U*YA'+XA*V'     G+U*YG'+XG*V'    ]
C         UU' H VV = [                                           ].
C                    [   Qout+U*YQ'+XQ*V'   op(Bout)+U*YB'+XB*V' ]
C
C     This is an auxiliary routine called by MB04TB.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LTRA    LOGICAL
C             Specifies the form of op( A ) as follows:
C             = .FALSE.:  op( A ) = A;
C             = .TRUE.:   op( A ) = A'.
C
C     LTRB    LOGICAL
C             Specifies the form of op( B ) as follows:
C             = .FALSE.:  op( B ) = B;
C             = .TRUE.:   op( B ) = B'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix Q. N >= 0.
C
C     K       (input) INTEGER
C             The offset of the reduction. Elements below the K-th
C             subdiagonal in the first NB columns of op(A) are
C             reduced to zero. K >= 0.
C
C     NB      (input) INTEGER
C             The number of columns/rows to be reduced. N > NB >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C                     (LDA,N)     if LTRA = .FALSE.
C                     (LDA,K+N)   if LTRA = .TRUE.
C             On entry with LTRA = .FALSE., the leading (K+N)-by-N part
C             of this array must contain the matrix A.
C             On entry with LTRA = .TRUE., the leading N-by-(K+N) part
C             of this array must contain the matrix A.
C             On exit with LTRA = .FALSE., the leading (K+N)-by-N part
C             of this array contains the matrix Aout and, in the zero
C             parts, information about the elementary reflectors used to
C             compute the reduction.
C             On exit with LTRA = .TRUE., the leading N-by-(K+N) part of
C             this array contains the matrix Aout and in the zero parts
C             information about the elementary reflectors.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,K+N),  if LTRA = .FALSE.;
C             LDA >= MAX(1,N),    if LTRA = .TRUE..
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C                     (LDB,K+N)   if LTRB = .FALSE.
C                     (LDB,N)     if LTRB = .TRUE.
C             On entry with LTRB = .FALSE., the leading N-by-(K+N) part
C             of this array must contain the matrix B.
C             On entry with LTRB = .TRUE., the leading (K+N)-by-N part
C             of this array must contain the matrix B.
C             On exit with LTRB = .FALSE., the leading N-by-(K+N) part
C             of this array contains the matrix Bout and, in the zero
C             parts, information about the elementary reflectors used to
C             compute the reduction.
C             On exit with LTRB = .TRUE., the leading (K+N)-by-N part of
C             this array contains the matrix Bout and in the zero parts
C             information about the elementary reflectors.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,N),    if LTRB = .FALSE.;
C             LDB >= MAX(1,K+N),  if LTRB = .TRUE..
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix G.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Gout.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix Q.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Qout and in the zero parts information about
C             the elementary reflectors used to compute the reduction.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1,N).
C
C     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB)
C             On exit, the leading N-by-(2*NB) part of this array
C             contains the matrix XA.
C
C     LDXA    INTEGER
C             The leading dimension of the array XA.  LDXA >= MAX(1,N).
C
C     XB      (output) DOUBLE PRECISION array, dimension (LDXB,2*NB)
C             On exit, the leading (K+N)-by-(2*NB) part of this array
C             contains the matrix XB.
C
C     LDXB    INTEGER
C             The leading dimension of the array XB. LDXB >= MAX(1,K+N).
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
C     YB      (output) DOUBLE PRECISION array, dimension (LDYB,2*NB)
C             On exit, the leading N-by-(2*NB) part of this array
C             contains the matrix YB.
C
C     LDYB    INTEGER
C             The leading dimension of the array YB.  LDYB >= MAX(1,N).
C
C     YG      (output) DOUBLE PRECISION array, dimension (LDYG,2*NB)
C             On exit, the leading (K+N)-by-(2*NB) part of this array
C             contains the matrix YG.
C
C     LDYG    INTEGER
C             The leading dimension of the array YG. LDYG >= MAX(1,K+N).
C
C     YQ      (output) DOUBLE PRECISION array, dimension (LDYQ,2*NB)
C             On exit, the leading N-by-(2*NB) part of this array
C             contains the matrix YQ.
C
C     LDYQ    INTEGER
C             The leading dimension of the array YQ.  LDYQ >= MAX(1,N).
C
C     CSL     (output) DOUBLE PRECISION array, dimension (2*NB)
C             On exit, the first 2NB elements of this array contain the
C             cosines and sines of the symplectic Givens rotations
C             applied from the left-hand side used to compute the
C             reduction.
C
C     CSR     (output) DOUBLE PRECISION array, dimension (2*NB)
C             On exit, the first 2NB-2 elements of this array contain
C             the cosines and sines of the symplectic Givens rotations
C             applied from the right-hand side used to compute the
C             reduction.
C
C     TAUL    (output) DOUBLE PRECISION array, dimension (NB)
C             On exit, the first NB elements of this array contain the
C             scalar factors of some of the elementary reflectors
C             applied form the left-hand side.
C
C     TAUR    (output) DOUBLE PRECISION array, dimension (NB)
C             On exit, the first NB-1 elements of this array contain the
C             scalar factors of some of the elementary reflectors
C             applied form the right-hand side.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (5*NB)
C
C     METHOD
C
C     For details regarding the representation of the orthogonal
C     symplectic matrices UU and VV within the arrays A, B, CSL, CSR, Q,
C     TAUL and TAUR see the description of MB04TB.
C
C     The contents of A, B, G and Q on exit are illustrated by the
C     following example with op(A) = A, op(B) = B, n = 5, k = 2 and
C     nb = 2:
C
C          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  )
C          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  )
C          ( r  r  r  r  r  )       ( r  r  r  r  r  r  r  )
C      A = ( u2 r  r  r  r  ),  G = ( r  r  r  r  r  r  r  ),
C          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  )
C          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  )
C          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  )
C
C          ( t  t  v1 v1 v1 )       ( r  r  r  r  r  v2 v2 )
C          ( u1 t  t  v1 v1 )       ( r  r  r  r  r  r  v2 )
C      Q = ( u1 u1 r  q  q  ),  B = ( b  b  b  r  r  b  b  ).
C          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  )
C          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  )
C
C     where a, b, g and q denote elements of the original matrices, r
C     denotes a modified element, t denotes a scalar factor of an
C     applied elementary reflector, ui and vi denote elements of the
C     matrices U and V, respectively.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires ( 16*K + 32*N + 42 )*N*NB +
C     ( 16*K + 112*N - 208/3*NB - 69 )*NB*NB - 29/3*NB floating point
C     operations and is numerically backward stable.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A numerically stable, structure preserving method for
C         computing the eigenvalues of real Hamiltonian or symplectic
C         pencils.
C         Numer. Math., Vol. 78 (3), pp. 329-358, 1998.
C
C     [2] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT Numerical Mathematics, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLASUB).
C
C     KEYWORDS
C
C     Elementary matrix operations, Matrix decompositions, Hamiltonian
C     matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           LTRA, LTRB
      INTEGER           K, LDA, LDB, LDG, LDQ, LDXA, LDXB, LDXG, LDXQ,
     $                  LDYA, LDYB, LDYG, LDYQ, N, NB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), CSL(*), CSR(*), DWORK(*),
     $                  G(LDG,*), Q(LDQ,*), TAUL(*), TAUR(*),
     $                  XA(LDXA,*), XB(LDXB,*), XG(LDXG,*), XQ(LDXQ,*),
     $                  YA(LDYA,*), YB(LDYB,*), YG(LDYG,*), YQ(LDYQ,*)
C     .. Local Scalars ..
      INTEGER           I, J, NB1, NB2, NB3, PDW
      DOUBLE PRECISION  ALPHA, C, S, TAUQ, TEMP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMV, DLARFG, DLARTG, DROT, DSCAL
C
C     .. Executable Statements ..
C
C     Quick return if possible.
C
      IF ( N+K.LE.0 ) THEN
         RETURN
      END IF
C
      NB1 = NB + 1
      NB2 = NB + NB
      NB3 = NB2 + NB
      PDW = NB3 + NB + 1
C
      IF ( LTRA.AND.LTRB ) THEN
         DO 90  I = 1, NB
C
C           Transform first row/column of A and Q. See routine MB04TS.
C
            ALPHA = Q(I,I)
            CALL DLARFG( N-I+1, ALPHA, Q(I+1,I), 1, TAUQ )
            Q(I,I) = ONE
            TEMP = -TAUQ*DDOT( N-I+1, Q(I,I), 1, A(I,K+I), LDA )
            CALL DAXPY( N-I+1, TEMP, Q(I,I), 1, A(I,K+I), LDA )
            TEMP = A(I,K+I)
            CALL DLARTG( TEMP, ALPHA, C, S, A(I,K+I) )
            CALL DLARFG( N-I+1, A(I,K+I), A(I,K+I+1), LDA, TAUL(I) )
            TEMP = A(I,K+I)
            A(I,K+I) = ONE
C
C           Update XQ with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  Q(I,I), 1, ZERO, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Q(I,1), LDQ,
     $                  Q(I,I), 1, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A(1,K+I), LDA,
     $                  Q(I,I), 1, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB+1), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,NB1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,1), LDQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, A(1,K+I), 1, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YQ(I,1), LDYQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YQ(I,NB1), LDYQ, ONE, Q(I,I+1), LDQ )
C
C           Update XA with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I+1, ONE, A(I+1,K+I),
     $                  LDA, Q(I,I), 1, ZERO, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB+1), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,NB1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update A(i+1:n,k+i).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,1), LDQ, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, A(1,K+I), 1, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YA(K+I,1), LDYA, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YA(K+I,NB1), LDYA, ONE, A(I+1,K+I), 1 )
C
C           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ].
C
            CALL DROT( N-I, A(I+1,K+I), 1, Q(I,I+1), LDQ, C, S )
C
C           Update XQ with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  A(I,K+I), LDA, ZERO, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK(NB2+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, A(1,K+I+1), LDA,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB3+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  A(I,K+I+1), LDA, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,NB1), LDYQ,
     $                  A(I,K+I+1), LDA, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XQ(I+1,I+NB), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XQ(I+1,I+NB), 1, Q(I,I+1), LDQ )
C
C           Update XA with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I+1, ONE, A(I+1,K+I),
     $                  LDA, A(I,K+I), LDA, ZERO, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  DWORK(NB2+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB3+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  A(I,K+I+1), LDA, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  A(I,K+I+1), LDA, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XA(I+1,I+NB), 1 )
C
C           Update A(i+1:n,k+i).
C
            CALL DAXPY( N-I, ONE, XA(I+1,I+NB), 1, A(I+1,K+I), 1 )
C
C           Update XG with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  Q(I,I), 1, ZERO, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB+1), 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,NB1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XG(1,I), 1 )
C
C           Update G(k+i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  Q(I,1), LDQ, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, A(1,K+I), 1, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YG(K+I,1), LDYG, ONE, G(K+I,K+I+1), LDG )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YG(K+I,NB1), LDYG, ONE, G(K+I,K+I+1), LDG )
C
C           Update XB with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I+1, ONE, B(1,I), LDB,
     $                  Q(I,I), 1, ZERO, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB, LDXB,
     $                  DWORK, 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB+1), 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,NB1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XB(1,I), 1 )
C
C           Update B(:,i).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  Q(I,1), LDQ, ONE, B(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, A(1,K+I), 1, ONE, B(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YB(I,1), LDYB, ONE, B(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YB(I,NB1), LDYB, ONE, B(K+I+1,I), 1 )
C
C           Apply rotation to [ G(k+i,:); B(:,i)' ].
C
            CALL DROT( K+N, G(K+I,1), LDG, B(1,I), 1, C, S )
C
            DO 10  J = 1, I-1
               YG(K+I,J) = ZERO
   10       CONTINUE
            DO 20  J = 1, I-1
               YG(K+I,NB+J) = ZERO
   20       CONTINUE
            DO 30  J = 1, I-1
               YA(K+I,J) = ZERO
   30       CONTINUE
            DO 40  J = 1, I-1
               YA(K+I,NB+J) = ZERO
   40       CONTINUE
C
C           Update XG with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  A(I,K+I), LDA, ZERO, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK(NB2+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB3+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,1), LDYG,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,NB1), LDYG,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XG(1,I+NB), 1 )
C
C           Update G(k+i,:).
C
            CALL DAXPY( K+N, ONE, XG(1,I+NB), 1, G(K+I,1), LDG )
C
C           Update XB with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I+1, ONE, B(1,I), LDB,
     $                  A(I,K+I), LDA, ZERO, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  DWORK(NB2+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB3+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,NB1), LDYB,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XB(1,I+NB), 1 )
C
C           Update B(:,i).
C
            CALL DAXPY( K+N, ONE, XB(1,I+NB), 1, B(1,I), 1 )
C
            A(I,K+I) = TEMP
            Q(I,I) = TAUQ
            CSL(2*I-1) = C
            CSL(2*I) = S
C
C           Transform first row/column of Q and B.
C
            ALPHA = Q(I,I+1)
            CALL DLARFG( N-I, ALPHA, Q(I,I+2), LDQ, TAUQ )
            Q(I,I+1) = ONE
            TEMP = -TAUQ*DDOT( N-I, Q(I,I+1), LDQ, B(K+I+1,I), 1 )
            CALL DAXPY( N-I, TEMP,  Q(I,I+1), LDQ, B(K+I+1,I), 1 )
            TEMP = B(K+I+1,I)
            CALL DLARTG( TEMP, ALPHA, C, S, B(K+I+1,I) )
            S = -S
            CALL DLARFG( N-I, B(K+I+1,I), B(K+I+2,I), 1, TAUR(I) )
            TEMP = B(K+I+1,I)
            B(K+I+1,I) = ONE
C
C           Update YB with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, B(K+I+1,I+1),
     $                  LDB, Q(I,I+1), LDQ, ZERO, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,NB1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  DWORK, 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB+1), 1, ONE, YB(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YB(I+1,I), 1 )
C
C           Update B(k+i+1,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XB(K+I+1,1), LDXB, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XB(K+I+1,NB1), LDXB, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  Q(1,I+1), 1, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, B(K+I+1,1), LDB, ONE, B(K+I+1,I+1), LDB )
C
C           Update YQ with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,NB1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK, 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB+1), 1, ONE, YQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YQ(I+1,I), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XQ(I+1,1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XQ(I+1,NB1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  Q(1,I+1), 1, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, B(K+I+1,1), LDB, ONE, Q(I+1,I+1), 1 )
C
C           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ].
C
            CALL DROT( N-I, Q(I+1,I+1), 1, B(K+I+1,I+1), LDB, C, S )
            DO 50  J = 1, I
               XB(K+I+1,J) = ZERO
   50       CONTINUE
            DO 60  J = 1, I
               XB(K+I+1,NB+J) = ZERO
   60       CONTINUE
C
C           Update YB with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, B(K+I+1,I+1),
     $                  LDB, B(K+I+1,I), 1, ZERO, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,1), LDXB,
     $                  B(K+I+2,I), 1, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,NB1), LDXB,
     $                  B(K+I+2,I), 1, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I, N-I-1, ONE, Q(1,I+2), LDQ,
     $                  B(K+I+2,I), 1, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  DWORK(NB2+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, B(K+I+2,1),
     $                  LDQ, B(K+I+2,I), 1, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB3+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YB(I+1,I+NB), 1 )
C
C           Update B(k+i+1,i+1:n).
C
            CALL DAXPY( N-I, ONE, YB(I+1,I+NB), 1, B(K+I+1,I+1), LDB )
C
C           Update YQ with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  B(K+I+1,I), 1, ZERO, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                  B(K+I+2,I), 1, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,NB1), LDXQ,
     $                  B(K+I+2,I), 1, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK(NB2+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB3+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YQ(I+1,I+NB), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DAXPY( N-I, ONE, YQ(I+1,I+NB), 1, Q(I+1,I+1), 1 )
C
C           Update YA with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, K+N, ONE, A(I+1,1), LDA,
     $                  Q(I,I+1), LDQ, ZERO, YA(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,NB1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB+1), 1, ONE, YA(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C
C           Update A(i+1,1:k+n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XA(I+1,1), LDXA, ONE, A(I+1,K+I+1), LDA )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(I+1,K+I+1), LDA )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  Q(1,I+1), 1, ONE, A(I+1,1), LDA )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  B(K+I+1,1), LDB, ONE, A(I+1,1), LDA )
C
C           Update YG with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  Q(I,I+1), LDQ, ZERO, YG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,NB1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG, LDYG,
     $                  DWORK, 1, ONE, YG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB+1), 1, ONE, YG(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YG(1,I), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XG(K+I+1,1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XG(K+I+1,NB1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  Q(1,I+1), 1, ONE, G(1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  B(K+I+1,1), LDB, ONE, G(1,K+I+1), 1 )
            DO 70  J = 1, I
               XG(K+I+1,J) = ZERO
   70       CONTINUE
            DO 80  J = 1, I
               XG(K+I+1,NB+J) = ZERO
   80       CONTINUE
C
C           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ].
C
            CALL DROT( K+N, A(I+1,1), LDA, G(1,K+I+1), 1, C, S )
C
C           Update YA with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, K+N, ONE, A(I+1,1), LDA,
     $                  B(K+I+1,I), 1, ZERO, YA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,NB1), LDXA,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  DWORK(NB2+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB3+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YA(1,I+NB), 1 )
C
C           Update A(i+1,1:k+n).
C
            CALL DAXPY( K+N, ONE, YA(1,I+NB), 1, A(I+1,1), LDA )
C
C           Update YG with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  B(K+I+1,I), 1, ZERO, YG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1), LDXG,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,NB1), LDXG,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  DWORK(NB2+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB3+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YG(1,I+NB), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DAXPY( K+N, ONE, YG(1,I+NB), 1, G(1,K+I+1), 1 )
C
            B(K+I+1,I) = TEMP
            Q(I,I+1) = TAUQ
            CSR(2*I-1) = C
            CSR(2*I) = S
   90    CONTINUE
      ELSE IF ( LTRA ) THEN
         DO 180  I = 1, NB
C
C           Transform first row/column of A and Q. See routine MB04TS.
C
            ALPHA = Q(I,I)
            CALL DLARFG( N-I+1, ALPHA, Q(I+1,I), 1, TAUQ )
            Q(I,I) = ONE
            TEMP = -TAUQ*DDOT( N-I+1, Q(I,I), 1, A(I,K+I), LDA )
            CALL DAXPY( N-I+1, TEMP, Q(I,I), 1, A(I,K+I), LDA )
            TEMP = A(I,K+I)
            CALL DLARTG( TEMP, ALPHA, C, S, A(I,K+I) )
            CALL DLARFG( N-I+1, A(I,K+I), A(I,K+I+1), LDA, TAUL(I) )
            TEMP = A(I,K+I)
            A(I,K+I) = ONE
C
C           Update XQ with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  Q(I,I), 1, ZERO, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Q(I,1), LDQ,
     $                  Q(I,I), 1, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A(1,K+I), LDA,
     $                  Q(I,I), 1, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB+1), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,NB1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,1), LDQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, A(1,K+I), 1, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YQ(I,1), LDYQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YQ(I,NB1), LDYQ, ONE, Q(I,I+1), LDQ )
C
C           Update XA with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I+1, ONE, A(I+1,K+I),
     $                  LDA, Q(I,I), 1, ZERO, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB+1), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,NB1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update A(i+1:n,k+i).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,1), LDQ, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, A(1,K+I), 1, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YA(K+I,1), LDYA, ONE, A(I+1,K+I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YA(K+I,NB1), LDYA, ONE, A(I+1,K+I), 1 )
C
C           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ].
C
            CALL DROT( N-I, A(I+1,K+I), 1, Q(I,I+1), LDQ, C, S )
C
C           Update XQ with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  A(I,K+I), LDA, ZERO, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK(NB2+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, A(1,K+I+1), LDA,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB3+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  A(I,K+I+1), LDA, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,NB1), LDYQ,
     $                  A(I,K+I+1), LDA, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XQ(I+1,I+NB), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XQ(I+1,I+NB), 1, Q(I,I+1), LDQ )
C
C           Update XA with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I+1, ONE, A(I+1,K+I),
     $                  LDA, A(I,K+I), LDA, ZERO, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  DWORK(NB2+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB3+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  A(I,K+I+1), LDA, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  A(I,K+I+1), LDA, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XA(I+1,I+NB), 1 )
C
C           Update A(i+1:n,k+i).
C
            CALL DAXPY( N-I, ONE, XA(I+1,I+NB), 1, A(I+1,K+I), 1 )
C
C           Update XG with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  Q(I,I), 1, ZERO, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB+1), 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,NB1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XG(1,I), 1 )
C
C           Update G(k+i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  Q(I,1), LDQ, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, A(1,K+I), 1, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YG(K+I,1), LDYG, ONE, G(K+I,K+I+1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YG(K+I,NB1), LDYG, ONE, G(K+I,K+I+1), LDG )
C
C           Update XB with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, B(I,1), LDB,
     $                  Q(I,I), 1, ZERO, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB, LDXB,
     $                  DWORK, 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB+1), 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,NB1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XB(1,I), 1 )
C
C           Update B(i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  Q(I,1), LDQ, ONE, B(I,1), LDB )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, A(1,K+I), 1, ONE, B(I,1), LDB )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YB(I,1), LDYB, ONE, B(I,K+I+1), LDB )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YB(I,NB1), LDYB, ONE, B(I,K+I+1), LDB )
C
C           Apply rotation to [ G(k+i,:); B(i,:) ].
C
            CALL DROT( K+N, G(K+I,1), LDG, B(I,1), LDB, C, S )
C
            DO 100  J = 1, I-1
               YG(K+I,J) = ZERO
  100       CONTINUE
            DO 110  J = 1, I-1
               YG(K+I,NB+J) = ZERO
  110       CONTINUE
            DO 120  J = 1, I-1
               YA(K+I,J) = ZERO
  120       CONTINUE
            DO 130  J = 1, I-1
               YA(K+I,NB+J) = ZERO
  130       CONTINUE
C
C           Update XG with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  A(I,K+I), LDA, ZERO, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK(NB2+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB3+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,1), LDYG,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,NB1), LDYG,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XG(1,I+NB), 1 )
C
C           Update G(k+i,:).
C
            CALL DAXPY( K+N, ONE, XG(1,I+NB), 1, G(K+I,1), LDG )
C
C           Update XB with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, B(I,1), LDB,
     $                  A(I,K+I), LDA, ZERO, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  DWORK(NB2+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB3+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,NB1), LDYB,
     $                  A(I,K+I+1), LDA, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XB(1,I+NB), 1 )
C
C           Update B(i,:).
C
            CALL DAXPY( K+N, ONE, XB(1,I+NB), 1, B(I,1), LDB )
C
            A(I,K+I) = TEMP
            Q(I,I) = TAUQ
            CSL(2*I-1) = C
            CSL(2*I) = S
C
C           Transform first rows of Q and B.
C
            ALPHA = Q(I,I+1)
            CALL DLARFG( N-I, ALPHA, Q(I,I+2), LDQ, TAUQ )
            Q(I,I+1) = ONE
            TEMP = -TAUQ*DDOT( N-I, Q(I,I+1), LDQ, B(I,K+I+1), LDB )
            CALL DAXPY( N-I, TEMP,  Q(I,I+1), LDQ, B(I,K+I+1), LDB )
            TEMP = B(I,K+I+1)
            CALL DLARTG( TEMP, ALPHA, C, S, B(I,K+I+1) )
            S = -S
            CALL DLARFG( N-I, B(I,K+I+1), B(I,K+I+2), LDB, TAUR(I) )
            TEMP = B(I,K+I+1)
            B(I,K+I+1) = ONE
C
C           Update YB with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, B(I+1,K+I+1),
     $                  LDB, Q(I,I+1), LDQ, ZERO, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,NB1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  DWORK, 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB+1), 1, ONE, YB(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YB(I+1,I), 1 )
C
C           Update B(i+1:n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XB(K+I+1,1), LDXB, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XB(K+I+1,NB1), LDXB, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  Q(1,I+1), 1, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, B(1,K+I+1), 1, ONE, B(I+1,K+I+1), 1 )
C
C           Update YQ with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,NB1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK, 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB+1), 1, ONE, YQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YQ(I+1,I), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XQ(I+1,1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XQ(I+1,NB1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  Q(1,I+1), 1, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, B(1,K+I+1), 1, ONE, Q(I+1,I+1), 1 )
C
C           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ].
C
            CALL DROT( N-I, Q(I+1,I+1), 1, B(I+1,K+I+1), 1, C, S )
            DO 140  J = 1, I
               XB(K+I+1,J) = ZERO
  140       CONTINUE
            DO 150  J = 1, I
               XB(K+I+1,NB+J) = ZERO
  150       CONTINUE
C
C           Update YB with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, B(I+1,K+I+1),
     $                  LDB, B(I,K+I+1), LDB, ZERO, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,1), LDXB,
     $                  B(I,K+I+2), LDB, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,NB1), LDXB,
     $                  B(I,K+I+2), LDB, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I, N-I-1, ONE, Q(1,I+2), LDQ,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  DWORK(NB2+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I-1, ONE, B(1,K+I+2),
     $                  LDQ, B(I,K+I+2), LDB, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB3+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YB(I+1,I+NB), 1 )
C
C           Update B(i+1:n,k+i+1).
C
            CALL DAXPY( N-I, ONE, YB(I+1,I+NB), 1, B(I+1,K+I+1), 1 )
C
C           Update YQ with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  B(I,K+I+1), LDB, ZERO, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                  B(I,K+I+2), LDB, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,NB1), LDXQ,
     $                  B(I,K+I+2), LDB, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK(NB2+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB3+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YQ(I+1,I+NB), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DAXPY( N-I, ONE, YQ(I+1,I+NB), 1, Q(I+1,I+1), 1 )
C
C           Update YA with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, K+N, ONE, A(I+1,1), LDA,
     $                  Q(I,I+1), LDQ, ZERO, YA(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,NB1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB+1), 1, ONE, YA(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C
C           Update A(i+1,1:k+n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XA(I+1,1), LDXA, ONE, A(I+1,K+I+1), LDA )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(I+1,K+I+1), LDA )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  Q(1,I+1), 1, ONE, A(I+1,1), LDA )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  B(1,K+I+1), 1, ONE, A(I+1,1), LDA )
C
C           Update YG with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  Q(I,I+1), LDQ, ZERO, YG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,NB1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG, LDYG,
     $                  DWORK, 1, ONE, YG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB+1), 1, ONE, YG(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YG(1,I), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XG(K+I+1,1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  XG(K+I+1,NB1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  Q(1,I+1), 1, ONE, G(1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  B(1,K+I+1), 1, ONE, G(1,K+I+1), 1 )
            DO 160  J = 1, I
               XG(K+I+1,J) = ZERO
  160       CONTINUE
            DO 170  J = 1, I
               XG(K+I+1,NB+J) = ZERO
  170       CONTINUE
C
C           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ].
C
            CALL DROT( K+N, A(I+1,1), LDA, G(1,K+I+1), 1, C, S )
C
C           Update YA with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, K+N, ONE, A(I+1,1), LDA,
     $                  B(I,K+I+1), LDB, ZERO, YA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,NB1), LDXA,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  DWORK(NB2+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB3+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YA(1,I+NB), 1 )
C
C           Update A(i+1,1:k+n).
C
            CALL DAXPY( K+N, ONE, YA(1,I+NB), 1, A(I+1,1), LDA )
C
C           Update YG with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  B(I,K+I+1), LDB, ZERO, YG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1), LDXG,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,NB1), LDXG,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I, N-I, ONE, A(1,K+I+1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  DWORK(NB2+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB3+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YG(1,I+NB), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DAXPY( K+N, ONE, YG(1,I+NB), 1, G(1,K+I+1), 1 )
C
            B(I,K+I+1) = TEMP
            Q(I,I+1) = TAUQ
            CSR(2*I-1) = C
            CSR(2*I) = S
  180    CONTINUE
C
      ELSE IF ( LTRB ) THEN
         DO 270  I = 1, NB
C
C           Transform first columns of A and Q. See routine MB04TS.
C
            ALPHA = Q(I,I)
            CALL DLARFG( N-I+1, ALPHA, Q(I+1,I), 1, TAUQ )
            Q(I,I) = ONE
            TEMP = -TAUQ*DDOT( N-I+1, Q(I,I), 1, A(K+I,I), 1 )
            CALL DAXPY( N-I+1, TEMP, Q(I,I), 1, A(K+I,I), 1 )
            TEMP = A(K+I,I)
            CALL DLARTG( TEMP, ALPHA, C, S, A(K+I,I) )
            CALL DLARFG( N-I+1, A(K+I,I), A(K+I+1,I), 1, TAUL(I) )
            TEMP = A(K+I,I)
            A(K+I,I) = ONE
C
C           Update XQ with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  Q(I,I), 1, ZERO, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Q(I,1), LDQ,
     $                  Q(I,I), 1, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, A(K+I,1), LDA,
     $                  Q(I,I), 1, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB+1), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,NB1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,1), LDQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, A(K+I,1), LDA, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YQ(I,1), LDYQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YQ(I,NB1), LDYQ, ONE, Q(I,I+1), LDQ )
C
C           Update XA with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, A(K+I,I+1), LDA,
     $                  Q(I,I), 1, ZERO, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB+1), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,NB1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update A(k+i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,1), LDQ, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, A(K+I,1), LDA, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YA(K+I,1), LDYA, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YA(K+I,NB1), LDYA, ONE, A(K+I,I+1), LDA )
C
C           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ].
C
            CALL DROT( N-I, A(K+I,I+1), LDA, Q(I,I+1), LDQ, C, S )
C
C           Update XQ with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  A(K+I,I), 1, ZERO, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  A(K+I+1,I), 1, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK(NB2+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  A(K+I+1,I), 1, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB3+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  A(K+I+1,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,NB1), LDYQ,
     $                  A(K+I+1,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XQ(I+1,I+NB), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XQ(I+1,I+NB), 1, Q(I,I+1), LDQ )
C
C           Update XA with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, A(K+I,I+1), LDA,
     $                  A(K+I,I), 1, ZERO, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  DWORK(NB2+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB3+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  A(K+I+1,I), 1, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  A(K+I+1,I), 1, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XA(I+1,I+NB), 1 )
C
C           Update A(k+i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XA(I+1,I+NB), 1, A(K+I,I+1), LDA )
C
C           Update XG with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  Q(I,I), 1, ZERO, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB+1), 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,NB1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XG(1,I), 1 )
C
C           Update G(k+i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  Q(I,1), LDQ, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, A(K+I,1), LDA, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YG(K+I,1), LDYG, ONE, G(K+I,K+I+1), LDG )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YG(K+I,NB1), LDYG, ONE, G(K+I,K+I+1), LDG )
C
C           Update XB with first Householder reflection.
C
            CALL DGEMV( 'No Transpose', K+N, N-I+1, ONE, B(1,I), LDB,
     $                  Q(I,I), 1, ZERO, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB, LDXB,
     $                  DWORK, 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB+1), 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,NB1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XB(1,I), 1 )
C
C           Update B(:,i).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  Q(I,1), LDQ, ONE, B(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, A(K+I,1), LDA, ONE, B(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YB(I,1), LDYB, ONE, B(K+I+1,I), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  YB(I,NB1), LDYB, ONE, B(K+I+1,I), 1 )
C
C           Apply rotation to [ G(k+i,:); B(:,i)' ].
C
            CALL DROT( K+N, G(K+I,1), LDG, B(1,I), 1, C, S )
C
            DO 190  J = 1, I-1
               YG(K+I,J) = ZERO
  190       CONTINUE
            DO 200  J = 1, I-1
               YG(K+I,NB+J) = ZERO
  200       CONTINUE
            DO 210  J = 1, I-1
               YA(K+I,J) = ZERO
  210       CONTINUE
            DO 220  J = 1, I-1
               YA(K+I,NB+J) = ZERO
  220       CONTINUE
C
C           Update XG with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  A(K+I,I), 1, ZERO, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK(NB2+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB3+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,1), LDYG,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,NB1), LDYG,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XG(1,I+NB), 1 )
C
C           Update G(k+i,:).
C
            CALL DAXPY( K+N, ONE, XG(1,I+NB), 1, G(K+I,1), LDG )
C
C           Update XB with second Householder reflection.
C
            CALL DGEMV( 'No Transpose', K+N, N-I+1, ONE, B(1,I), LDB,
     $                  A(K+I,I), 1, ZERO, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  DWORK(NB2+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB3+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,NB1), LDYB,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XB(1,I+NB), 1 )
C
C           Update B(:,i).
C
            CALL DAXPY( K+N, ONE, XB(1,I+NB), 1, B(1,I), 1 )
C
            A(K+I,I) = TEMP
            Q(I,I) = TAUQ
            CSL(2*I-1) = C
            CSL(2*I) = S
C
C           Transform first rows of Q and B.
C
            ALPHA = Q(I,I+1)
            CALL DLARFG( N-I, ALPHA, Q(I,I+2), LDQ, TAUQ )
            Q(I,I+1) = ONE
            TEMP = -TAUQ*DDOT( N-I, Q(I,I+1), LDQ, B(K+I+1,I), 1 )
            CALL DAXPY( N-I, TEMP,  Q(I,I+1), LDQ, B(K+I+1,I), 1 )
            TEMP = B(K+I+1,I)
            CALL DLARTG( TEMP, ALPHA, C, S, B(K+I+1,I) )
            S = -S
            CALL DLARFG( N-I, B(K+I+1,I), B(K+I+2,I), 1, TAUR(I) )
            TEMP = B(K+I+1,I)
            B(K+I+1,I) = ONE
C
C           Update YB with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, B(K+I+1,I+1),
     $                  LDB, Q(I,I+1), LDQ, ZERO, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,NB1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  DWORK, 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, B(K+I+1,1), LDB,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB+1), 1, ONE, YB(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YB(I+1,I), 1 )
C
C           Update B(k+i+1,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XB(K+I+1,1), LDXB, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XB(K+I+1,NB1), LDXB, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  Q(1,I+1), 1, ONE, B(K+I+1,I+1), LDB )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, B(K+I+1,1), LDB, ONE, B(K+I+1,I+1), LDB )
C
C           Update YQ with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,NB1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK, 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB+1), 1, ONE, YQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YQ(I+1,I), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XQ(I+1,1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XQ(I+1,NB1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  Q(1,I+1), 1, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, B(K+I+1,1), LDB, ONE, Q(I+1,I+1), 1 )
C
C           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ].
C
            CALL DROT( N-I, Q(I+1,I+1), 1, B(K+I+1,I+1), LDB, C, S )
            DO 230  J = 1, I
               XB(K+I+1,J) = ZERO
  230       CONTINUE
            DO 240  J = 1, I
               XB(K+I+1,NB+J) = ZERO
  240       CONTINUE
C
C           Update YB with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I, N-I, ONE, B(K+I+1,I+1),
     $                  LDB, B(K+I+1,I), 1, ZERO, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,1), LDXB,
     $                  B(K+I+2,I), 1, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,NB1), LDXB,
     $                  B(K+I+2,I), 1, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I, N-I-1, ONE, Q(1,I+2), LDQ,
     $                  B(K+I+2,I), 1, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  DWORK(NB2+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I-1, ONE, B(K+I+2,1),
     $                  LDQ, B(K+I+2,I), 1, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB3+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YB(I+1,I+NB), 1 )
C
C           Update B(k+i+1,i+1:n).
C
            CALL DAXPY( N-I, ONE, YB(I+1,I+NB), 1, B(K+I+1,I+1), LDB )
C
C           Update YQ with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  B(K+I+1,I), 1, ZERO, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                  B(K+I+2,I), 1, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,NB1), LDXQ,
     $                  B(K+I+2,I), 1, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK(NB2+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB3+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YQ(I+1,I+NB), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DAXPY( N-I, ONE, YQ(I+1,I+NB), 1, Q(I+1,I+1), 1 )
C
C           Update YA with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  Q(I,I+1), LDQ, ZERO, YA(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,NB1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB+1), 1, ONE, YA(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C
C           Update A(1:k+n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XA(I+1,1), LDXA, ONE, A(K+I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(K+I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  Q(1,I+1), 1, ONE, A(1,I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  B(K+I+1,1), LDB, ONE, A(1,I+1), 1 )
C
C           Update YG with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  Q(I,I+1), LDQ, ZERO, YG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,NB1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG, LDYG,
     $                  DWORK, 1, ONE, YG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB+1), 1, ONE, YG(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YG(1,I), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XG(K+I+1,1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XG(K+I+1,NB1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  Q(1,I+1), 1, ONE, G(1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  B(K+I+1,1), LDB, ONE, G(1,K+I+1), 1 )
            DO 250  J = 1, I
               XG(K+I+1,J) = ZERO
  250       CONTINUE
            DO 260  J = 1, I
               XG(K+I+1,NB+J) = ZERO
  260       CONTINUE
C
C           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ].
C
            CALL DROT( K+N, A(1,I+1), 1, G(1,K+I+1), 1, C, S )
C
C           Update YA with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  B(K+I+1,I), 1, ZERO, YA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,NB1), LDXA,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  DWORK(NB2+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB3+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YA(1,I+NB), 1 )
C
C           Update A(1:k+n,i+1).
C
            CALL DAXPY( K+N, ONE, YA(1,I+NB), 1, A(1,I+1), 1 )
C
C           Update YG with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  B(K+I+1,I), 1, ZERO, YG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1), LDXG,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,NB1), LDXG,
     $                  B(K+I+2,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  DWORK(NB2+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB3+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YG(1,I+NB), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DAXPY( K+N, ONE, YG(1,I+NB), 1, G(1,K+I+1), 1 )
C
            B(K+I+1,I) = TEMP
            Q(I,I+1) = TAUQ
            CSR(2*I-1) = C
            CSR(2*I) = S
  270    CONTINUE
C
      ELSE
         DO 360  I = 1, NB
C
C           Transform first columns of A and Q. See routine MB04TS.
C
            ALPHA = Q(I,I)
            CALL DLARFG( N-I+1, ALPHA, Q(I+1,I), 1, TAUQ )
            Q(I,I) = ONE
            TEMP = -TAUQ*DDOT( N-I+1, Q(I,I), 1, A(K+I,I), 1 )
            CALL DAXPY( N-I+1, TEMP, Q(I,I), 1, A(K+I,I), 1 )
            TEMP = A(K+I,I)
            CALL DLARTG( TEMP, ALPHA, C, S, A(K+I,I) )
            CALL DLARFG( N-I+1, A(K+I,I), A(K+I+1,I), 1, TAUL(I) )
            TEMP = A(K+I,I)
            A(K+I,I) = ONE
C
C           Update XQ with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  Q(I,I), 1, ZERO, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Q(I,1), LDQ,
     $                  Q(I,I), 1, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK, 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, A(K+I,1), LDA,
     $                  Q(I,I), 1, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB+1), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I), 1, ONE, XQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YQ(I,NB1), LDYQ,
     $                  Q(I,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XQ(I+1,I), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,1), LDQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, A(K+I,1), LDA, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YQ(I,1), LDYQ, ONE, Q(I,I+1), LDQ )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YQ(I,NB1), LDYQ, ONE, Q(I,I+1), LDQ )
C
C           Update XA with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, A(K+I,I+1), LDA,
     $                  Q(I,I), 1, ZERO, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,1), LDXA,
     $                  DWORK, 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB+1), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YA(K+I,NB1), LDYA,
     $                  Q(I,I), 1, ZERO, XA(1,I), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XA(1,I), 1, ONE, XA(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, XA(I+1,I), 1 )
C
C           Update A(k+i,i+1:n).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,1), LDQ, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, A(K+I,1), LDA, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YA(K+I,1), LDYA, ONE, A(K+I,I+1), LDA )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YA(K+I,NB1), LDYA, ONE, A(K+I,I+1), LDA )
C
C           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ].
C
            CALL DROT( N-I, A(K+I,I+1), LDA, Q(I,I+1), LDQ, C, S )
C
C           Update XQ with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, Q(I,I+1), LDQ,
     $                  A(K+I,I), 1, ZERO, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  A(K+I+1,I), 1, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  DWORK(NB2+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, A(K+I+1,1), LDA,
     $                  A(K+I+1,I), 1, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XQ(I+1,NB1),
     $                  LDXQ, DWORK(NB3+1), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  A(K+I+1,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YQ(I+1,NB1), LDYQ,
     $                  A(K+I+1,I), 1, ZERO, XQ(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XQ(1,I+NB), 1, ONE, XQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XQ(I+1,I+NB), 1 )
C
C           Update Q(i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XQ(I+1,I+NB), 1, Q(I,I+1), LDQ )
C
C           Update XA with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, N-I, ONE, A(K+I,I+1), LDA,
     $                  A(K+I,I), 1, ZERO, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  DWORK(NB2+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, XA(I+1,NB1),
     $                  LDXA, DWORK(NB3+1), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,1), LDYA,
     $                  A(K+I+1,I), 1, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YA(K+I+1,NB1), LDYA,
     $                  A(K+I+1,I), 1, ZERO, XA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  XA(1,I+NB), 1, ONE, XA(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUL(I), XA(I+1,I+NB), 1 )
C
C           Update A(k+i,i+1:n).
C
            CALL DAXPY( N-I, ONE, XA(I+1,I+NB), 1, A(K+I,I+1), LDA )
C
C           Update XG with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  Q(I,I), 1, ZERO, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG, LDXG,
     $                  DWORK, 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB+1), 1, ONE, XG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YG(K+I,NB1), LDYG,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XG(1,I), 1 )
C
C           Update G(k+i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  Q(I,1), LDQ, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, A(K+I,1), LDA, ONE, G(K+I,1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YG(K+I,1), LDYG, ONE, G(K+I,K+I+1), LDG )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YG(K+I,NB1), LDYG, ONE, G(K+I,K+I+1), LDG )
C
C           Update XB with first Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, B(I,1), LDB,
     $                  Q(I,I), 1, ZERO, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB, LDXB,
     $                  DWORK, 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB+1), 1, ONE, XB(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, YB(I,NB1), LDYB,
     $                  Q(I,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, XB(1,I), 1 )
C
C           Update B(i,:).
C
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  Q(I,1), LDQ, ONE, B(I,1), LDB )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, A(K+I,1), LDA, ONE, B(I,1), LDB )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  YB(I,1), LDYB, ONE, B(I,K+I+1), LDB )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  YB(I,NB1), LDYB, ONE, B(I,K+I+1), LDB )
C
C           Apply rotation to [ G(k+i,:); B(i,:) ].
C
            CALL DROT( K+N, G(K+I,1), LDG, B(I,1), LDB, C, S )
C
            DO 280  J = 1, I-1
               YG(K+I,J) = ZERO
  280       CONTINUE
            DO 290  J = 1, I-1
               YG(K+I,NB+J) = ZERO
  290       CONTINUE
            DO 300  J = 1, I-1
               YA(K+I,J) = ZERO
  300       CONTINUE
            DO 310  J = 1, I-1
               YA(K+I,NB+J) = ZERO
  310       CONTINUE
C
C           Update XG with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, G(K+I,1), LDG,
     $                  A(K+I,I), 1, ZERO, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XG, LDXG,
     $                  DWORK(NB2+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XG(1,NB1),
     $                  LDXG, DWORK(NB3+1), 1, ONE, XG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,1), LDYG,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YG(K+I+1,NB1), LDYG,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XG(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XG(1,I+NB), 1 )
C
C           Update G(k+i,:).
C
            CALL DAXPY( K+N, ONE, XG(1,I+NB), 1, G(K+I,1), LDG )
C
C           Update XB with second Householder reflection.
C
            CALL DGEMV( 'Transpose', N-I+1, K+N, ONE, B(I,1), LDB,
     $                  A(K+I,I), 1, ZERO, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, XB, LDXB,
     $                  DWORK(NB2+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, XB(1,NB1),
     $                  LDXB, DWORK(NB3+1), 1, ONE, XB(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I, I-1, ONE, YB(I+1,NB1), LDYB,
     $                  A(K+I+1,I), 1, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'Transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  DWORK(PDW), 1, ONE, XB(K+I+1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUL(I), XB(1,I+NB), 1 )
C
C           Update B(i,:).
C
            CALL DAXPY( K+N, ONE, XB(1,I+NB), 1, B(I,1), LDB )
C
            A(K+I,I) = TEMP
            Q(I,I) = TAUQ
            CSL(2*I-1) = C
            CSL(2*I) = S
C
C           Transform first rows of Q and B.
C
            ALPHA = Q(I,I+1)
            CALL DLARFG( N-I, ALPHA, Q(I,I+2), LDQ, TAUQ )
            Q(I,I+1) = ONE
            TEMP = -TAUQ*DDOT( N-I, Q(I,I+1), LDQ, B(I,K+I+1), LDB )
            CALL DAXPY( N-I, TEMP,  Q(I,I+1), LDQ, B(I,K+I+1), LDB )
            TEMP = B(I,K+I+1)
            CALL DLARTG( TEMP, ALPHA, C, S, B(I,K+I+1) )
            S = -S
            CALL DLARFG( N-I, B(I,K+I+1), B(I,K+I+2), LDB, TAUR(I) )
            TEMP = B(I,K+I+1)
            B(I,K+I+1) = ONE
C
C           Update YB with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, B(I+1,K+I+1),
     $                  LDB, Q(I,I+1), LDQ, ZERO, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XB(K+I+1,NB1), LDXB,
     $                  Q(I,I+1), LDQ, ZERO, YB(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YB(1,I), 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, Q(1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, DWORK, 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,1), LDYB,
     $                  DWORK, 1, ONE, YB(I+1,I), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I, ONE, B(1,K+I+1), LDB,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(NB+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB+1), 1, ONE, YB(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YB(I+1,I), 1 )
C
C           Update B(i+1:n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XB(K+I+1,1), LDXB, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XB(K+I+1,NB1), LDXB, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  Q(1,I+1), 1, ONE, B(I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, B(1,K+I+1), 1, ONE, B(I+1,K+I+1), 1 )
C
C           Update YQ with first Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XQ(I+1,NB1), LDXQ,
     $                  Q(I,I+1), LDQ, ZERO, YQ(1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YQ(1,I), 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK, 1, ONE, YQ(I+1,I), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB+1), 1, ONE, YQ(I+1,I), 1 )
            CALL DSCAL( N-I, -TAUQ, YQ(I+1,I), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XQ(I+1,1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XQ(I+1,NB1), LDXQ, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  Q(1,I+1), 1, ONE, Q(I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, B(1,K+I+1), 1, ONE, Q(I+1,I+1), 1 )
C
C           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ].
C
            CALL DROT( N-I, Q(I+1,I+1), 1, B(I+1,K+I+1), 1, C, S )
            DO 320  J = 1, I
               XB(K+I+1,J) = ZERO
  320       CONTINUE
            DO 330  J = 1, I
               XB(K+I+1,NB+J) = ZERO
  330       CONTINUE
C
C           Update YB with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, B(I+1,K+I+1),
     $                  LDB, B(I,K+I+1), LDB, ZERO, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,1), LDXB,
     $                  B(I,K+I+2), LDB, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XB(K+I+2,NB1), LDXB,
     $                  B(I,K+I+2), LDB, ZERO, YB(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YB(1,I+NB), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I, N-I-1, ONE, Q(1,I+2), LDQ,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(NB2+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YB(I+1,1), LDYB,
     $                  DWORK(NB2+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', I-1, N-I-1, ONE, B(1,K+I+2),
     $                  LDQ, B(I,K+I+2), LDB, ZERO, DWORK(NB3+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YB(I+1,NB1),
     $                  LDYB, DWORK(NB3+1), 1, ONE, YB(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YB(I+1,I+NB), 1 )
C
C           Update B(i+1:n,k+i+1).
C
            CALL DAXPY( N-I, ONE, YB(I+1,I+NB), 1, B(I+1,K+I+1), 1 )
C
C           Update YQ with second Householder reflection.
C
            CALL DGEMV( 'No transpose', N-I, N-I, ONE, Q(I+1,I+1), LDQ,
     $                  B(I,K+I+1), LDB, ZERO, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,1), LDXQ,
     $                  B(I,K+I+2), LDB, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XQ(I+2,NB1), LDXQ,
     $                  B(I,K+I+2), LDB, ZERO, YQ(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  YQ(1,I+NB), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, YQ(I+1,1), LDYQ,
     $                  DWORK(NB2+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', N-I, I-1, ONE, YQ(I+1,NB1),
     $                  LDYQ, DWORK(NB3+1), 1, ONE, YQ(I+1,I+NB), 1 )
            CALL DSCAL( N-I, -TAUR(I), YQ(I+1,I+NB), 1 )
C
C           Update Q(i+1:n,i+1).
C
            CALL DAXPY( N-I, ONE, YQ(I+1,I+NB), 1, Q(I+1,I+1), 1 )
C
C           Update YA with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  Q(I,I+1), LDQ, ZERO, YA(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XA(I+1,NB1), LDXA,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA, LDYA,
     $                  DWORK, 1, ONE, YA(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB+1), 1, ONE, YA(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YA(1,I), 1 )
C
C           Update A(1:k+n,i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XA(I+1,1), LDXA, ONE, A(K+I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XA(I+1,NB1), LDXA, ONE, A(K+I+1,I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  Q(1,I+1), 1, ONE, A(1,I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  B(1,K+I+1), 1, ONE, A(1,I+1), 1 )
C
C           Update YG with first Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  Q(I,I+1), LDQ, ZERO, YG(1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'Transpose', N-I, I, ONE, XG(K+I+1,NB1), LDXG,
     $                  Q(I,I+1), LDQ, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG, LDYG,
     $                  DWORK, 1, ONE, YG(1,I), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB+1), 1, ONE, YG(1,I), 1 )
            CALL DSCAL( K+N, -TAUQ, YG(1,I), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  XG(K+I+1,1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  XG(K+I+1,NB1), LDXG, ONE, G(K+I+1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  Q(1,I+1), 1, ONE, G(1,K+I+1), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  B(1,K+I+1), 1, ONE, G(1,K+I+1), 1 )
            DO 340  J = 1, I
               XG(K+I+1,J) = ZERO
  340       CONTINUE
            DO 350  J = 1, I
               XG(K+I+1,NB+J) = ZERO
  350       CONTINUE
C
C           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ].
C
            CALL DROT( K+N, A(1,I+1), 1, G(1,K+I+1), 1, C, S )
C
C           Update YA with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, A(1,I+1), LDA,
     $                  B(I,K+I+1), LDB, ZERO, YA(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,1), LDXA,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XA(I+2,NB1), LDXA,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YA(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YA, LDYA,
     $                  DWORK(NB2+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YA(1,NB1), LDYA,
     $                  DWORK(NB3+1), 1, ONE, YA(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YA(1,I+NB), 1 )
C
C           Update A(1:k+n,i+1).
C
            CALL DAXPY( K+N, ONE, YA(1,I+NB), 1, A(1,I+1), 1 )
C
C           Update YG with second Householder reflection.
C
            CALL DGEMV( 'No transpose', K+N, N-I, ONE, G(1,K+I+1), LDG,
     $                  B(I,K+I+1), LDB, ZERO, YG(1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,1), LDXG,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, Q(I+1,1), LDQ,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'Transpose', N-I-1, I, ONE, XG(K+I+2,NB1), LDXG,
     $                  B(I,K+I+2), LDB, ZERO, DWORK(PDW), 1 )
            CALL DGEMV( 'No transpose', N-I, I, ONE, A(K+I+1,1), LDA,
     $                  DWORK(PDW), 1, ONE, YG(K+I+1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I, ONE, YG, LDYG,
     $                  DWORK(NB2+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DGEMV( 'No transpose', K+N, I-1, ONE, YG(1,NB1), LDYG,
     $                  DWORK(NB3+1), 1, ONE, YG(1,I+NB), 1 )
            CALL DSCAL( K+N, -TAUR(I), YG(1,I+NB), 1 )
C
C           Update G(1:k+n,k+i+1).
C
            CALL DAXPY( K+N, ONE, YG(1,I+NB), 1, G(1,K+I+1), 1 )
C
            B(I,K+I+1) = TEMP
            Q(I,I+1) = TAUQ
            CSR(2*I-1) = C
            CSR(2*I) = S
  360    CONTINUE
      END IF
C
      RETURN
C *** Last line of MB03XU ***
      END
