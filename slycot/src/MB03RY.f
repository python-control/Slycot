      SUBROUTINE MB03RY( M, N, PMAX, A, LDA, B, LDB, C, LDC, INFO )
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
C     To solve the Sylvester equation -AX + XB = C, where A and B are
C     M-by-M and N-by-N matrices, respectively, in real Schur form.
C
C     This routine is intended to be called only by SLICOT Library
C     routine MB03RD. For efficiency purposes, the computations are
C     aborted when the infinity norm of an elementary submatrix of X is
C     greater than a given value PMAX.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of the matrix A and the number of rows of the
C             matrices C and X.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix B and the number of columns of the
C             matrices C and X.  N >= 0.
C
C     PMAX    (input) DOUBLE PRECISION
C             An upper bound for the infinity norm of an elementary
C             submatrix of X (see METHOD).
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,M)
C             The leading M-by-M part of this array must contain the
C             matrix A of the Sylvester equation, in real Schur form.
C             The elements below the real Schur form are not referenced.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,N)
C             The leading N-by-N part of this array must contain the
C             matrix B of the Sylvester equation, in real Schur form.
C             The elements below the real Schur form are not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix C of the Sylvester equation.
C             On exit, if INFO = 0, the leading M-by-N part of this
C             array contains the solution matrix X of the Sylvester
C             equation, and each elementary submatrix of X (see METHOD)
C             has the infinity norm less than or equal to PMAX.
C             On exit, if INFO = 1, the solution matrix X has not been
C             computed completely, because an elementary submatrix of X
C             had the infinity norm greater than PMAX. Part of the
C             matrix C has possibly been overwritten with the
C             corresponding part of X.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  an elementary submatrix of X had the infinity norm
C                   greater than the given value PMAX.
C
C     METHOD
C
C     The routine uses an adaptation of the standard method for solving
C     Sylvester equations [1], which controls the magnitude of the
C     individual elements of the computed solution [2]. The equation
C     -AX + XB = C can be rewritten as
C                                 p            l-1
C       -A  X   + X  B   = C   + sum  A  X   - sum  X  B
C         kk kl    kl ll    kl  i=k+1  ki il   j=1   kj jl
C
C     for l = 1:q, and k = p:-1:1, where A  , B  , C  , and X  , are
C                                         kk   ll   kl       kl
C     block submatrices defined by the partitioning induced by the Schur
C     form of A and B, and p and q are the numbers of the diagonal
C     blocks of A and B, respectively. So, the elementary submatrices of
C     X are found block column by block column, starting from the
C     bottom. If any such elementary submatrix has the infinity norm
C     greater than the given value PMAX, the calculations are ended.
C
C     REFERENCES
C
C     [1] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Bavely, C. and Stewart, G.W.
C         An Algorithm for Computing Reducing Subspaces by Block
C         Diagonalization.
C         SIAM J. Numer. Anal., 16, pp. 359-367, 1979.
C
C     NUMERICAL ASPECTS
C                               2      2
C     The algorithm requires 0(M N + MN ) operations.
C
C     FURTHER COMMENTS
C
C     Let
C
C            ( A   C )       ( I   X )
C        M = (       ),  Y = (       ).
C            ( 0   B )       ( 0   I )
C
C     Then
C
C         -1      ( A   0 )
C        Y  M Y = (       ),
C                 ( 0   B )
C
C     hence Y is an non-orthogonal transformation matrix which performs
C     the reduction of M to a block-diagonal form. Bounding a norm of
C     X is equivalent to setting an upper bound to the condition number
C     of the transformation matrix Y.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998.
C     Based on the RASP routine SYLSM by A. Varga, German Aerospace
C     Center, DLR Oberpfaffenhofen.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Diagonalization, real Schur form, Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, M, N
      DOUBLE PRECISION  PMAX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      INTEGER           DK, DL, I, IERR, J, K, KK, KK1, L, LL, LM1
      DOUBLE PRECISION  PNORM, SCALE
C     .. Local Arrays ..
      DOUBLE PRECISION  P(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGEMV, DLASY2
C     .. Executable Statements ..
C
C     For efficiency reasons, this routine does not check the input
C     parameters for errors.
C
      INFO = 0
C
C     Column loop indexed by L.
C
      L = 1
C     WHILE ( L.LE.N ) DO
   10 IF ( L.LE.N ) THEN
         LM1 = L - 1
         DL = 1
         IF ( L.LT.N ) THEN
            IF ( B(L+1,L).NE.ZERO )
     $         DL = 2
         ENDIF
         LL = LM1 + DL
         IF ( LM1.GT.0 ) THEN
C
C           Update one (or two) column(s) of C.
C
            IF ( DL.EQ.2 ) THEN
               CALL DGEMM( 'No transpose', 'No transpose', M, DL, LM1,
     $                     -ONE, C, LDC, B(1,L), LDB, ONE, C(1,L), LDC )
            ELSE
               CALL DGEMV( 'No transpose', M, LM1, -ONE, C, LDC, B(1,L),
     $                     1, ONE, C(1,L), 1 )
            END IF
         ENDIF
C
C        Row loop indexed by KK.
C
         KK = M
C        WHILE ( KK.GE.1 ) DO
   20    IF ( KK.GE.1 ) THEN
            KK1 = KK + 1
            DK = 1
            IF ( KK.GT.1 ) THEN
               IF ( A(KK,KK-1).NE.ZERO )
     $            DK = 2
            ENDIF
            K = KK1 - DK
            IF ( K.LT.M ) THEN
C
C              Update an elementary submatrix of C.
C
               DO 40 J = L, LL
C
                  DO 30 I = K, KK
                     C(I,J) = C(I,J) +
     $                        DDOT( M-KK, A(I,KK1), LDA, C(KK1,J), 1 )
   30             CONTINUE
C
   40          CONTINUE
C
            ENDIF
            CALL DLASY2( .FALSE., .FALSE., -1, DK, DL, A(K,K), LDA,
     $                   B(L,L), LDB, C(K,L), LDC, SCALE, P, DK, PNORM,
     $                   IERR )
            IF( SCALE.NE.ONE .OR. PNORM.GT.PMAX ) THEN
                INFO = 1
                RETURN
            END IF
            C(K,L) = -P(1)
            IF ( DL.EQ.1 ) THEN
               IF ( DK.EQ.2 )
     $            C(KK,L)  = -P(2)
            ELSE
               IF ( DK.EQ.1 ) THEN
                  C(K,LL)  = -P(2)
               ELSE
                  C(KK,L)  = -P(2)
                  C(K,LL)  = -P(3)
                  C(KK,LL) = -P(4)
               ENDIF
            ENDIF
            KK = KK - DK
            GO TO 20
         END IF
C        END WHILE 20
         L = L + DL
         GO TO 10
      END IF
C     END WHILE 10
      RETURN
C *** Last line of MB03RY ***
      END
