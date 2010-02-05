      SUBROUTINE MB02NY( UPDATU, UPDATV, M, N, I, K, Q, E, U, LDU, V,
     $                   LDV, DWORK )
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
C     To separate a zero singular value of a bidiagonal submatrix of
C     order k, k <= p, of the bidiagonal matrix
C
C               |Q(1) E(1)  0    ...   0   |
C               | 0   Q(2) E(2)        .   |
C           J = | .                    .   |
C               | .                  E(p-1)|
C               | 0   ...  ...   ...  Q(p) |
C
C     with p = MIN(M,N), by annihilating one or two superdiagonal
C     elements E(i-1) (if i > 1) and/or E(i) (if i < k).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPDATU  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix U the left-hand Givens rotations S, as follows:
C             = .FALSE.:  Do not form U;
C             = .TRUE. :  The given matrix U is updated (postmultiplied)
C                         by the left-hand Givens rotations S.
C
C     UPDATV  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix V the right-hand Givens rotations T, as follows:
C             = .FALSE.:  Do not form V;
C             = .TRUE. :  The given matrix V is updated (postmultiplied)
C                         by the right-hand Givens rotations T.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix U.  M >= 0.
C
C     N       (input) INTEGER
C             The number of rows of the matrix V.  N >= 0.
C
C     I       (input) INTEGER
C             The index of the negligible diagonal entry Q(I) of the
C             bidiagonal matrix J, I <= p.
C
C     K       (input) INTEGER
C             The index of the last diagonal entry of the considered
C             bidiagonal submatrix of J, i.e., E(K-1) is considered
C             negligible, K <= p.
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (p)
C             where p = MIN(M,N).
C             On entry, Q must contain the diagonal entries of the
C             bidiagonal matrix J.
C             On exit, Q contains the diagonal entries of the
C             transformed bidiagonal matrix S' J T.
C
C     E       (input/output) DOUBLE PRECISION array, dimension (p-1)
C             On entry, E must contain the superdiagonal entries of J.
C             On exit, E contains the superdiagonal entries of the
C             transformed bidiagonal matrix S' J T.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,p)
C             On entry, if UPDATU = .TRUE., U must contain the M-by-p
C             left transformation matrix.
C             On exit, if UPDATU = .TRUE., the Givens rotations S on the
C             left, annihilating E(i) if i < k, have been postmultiplied
C             into U.
C             U is not referenced if UPDATU = .FALSE..
C
C     LDU     INTEGER
C             The leading dimension of the array U.
C             LDU >= max(1,M) if UPDATU = .TRUE.;
C             LDU >= 1        if UPDATU = .FALSE..
C
C     V       (input/output) DOUBLE PRECISION array, dimension (LDV,p)
C             On entry, if UPDATV = .TRUE., V must contain the N-by-p
C             right transformation matrix.
C             On exit, if UPDATV = .TRUE., the Givens rotations T on the
C             right, annihilating E(i-1) if i > 1,  have been
C             postmultiplied into V.
C             V is not referenced if UPDATV = .FALSE..
C
C     LDV     INTEGER
C             The leading dimension of the array V.
C             LDV >= max(1,N) if UPDATV = .TRUE.;
C             LDV >= 1        if UPDATV = .FALSE..
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (MAX(1,LDWORK))
C             LDWORK >= 2*MAX(K-I,I-1),  if UPDATV = UPDATU = .TRUE.;
C             LDWORK >= 2*(K-I), if UPDATU = .TRUE., UPDATV = .FALSE.;
C             LDWORK >= 2*(I-1), if UPDATV = .TRUE., UPDATU = .FALSE.;
C             LDWORK >= 1,       if UPDATU = UPDATV = .FALSE..
C
C     METHOD
C
C     Let the considered bidiagonal submatrix be
C
C               |Q(1) E(1)  0                    ...   0   |
C               | 0   Q(2) E(2)                        .   |
C               | .                                    .   |
C               | .           Q(i-1) E(i-1)            .   |
C          Jk = | .                   Q(i) E(i)        .   |.
C               | .                       Q(i+1) .     .   |
C               | .                              ..    .   |
C               | .                                  E(k-1)|
C               | 0    ...                       ...  Q(k) |
C
C     A zero singular value of Jk manifests itself by a zero diagonal
C     entry Q(i) or in practice, a negligible value of Q(i).
C     When a negligible diagonal element Q(i) in Jk is present, the
C     bidiagonal submatrix Jk is split by the routine into 2 or 3
C     unreduced bidiagonal submatrices by annihilating E(i) (if i < k)
C     using Givens rotations S on the left and by annihilating E(i-1)
C     (if i > 1) using Givens rotations T on the right until Jk is
C     reduced to the form:
C
C               |Q(1) E(1)  0                ...   0   |
C               | 0         .                ...   .   |
C               | .                          ...   .   |
C               | .       Q(i-1) 0                 .   |
C     S' Jk T = | .              0   0             .   |.
C               | .                 Q(i+1)   .     .   |
C               | .                          ..    .   |
C               | .                              E(k-1)|
C               | 0    ...                   ...  Q(k) |
C
C     For more details, see [1, pp.11.12-11.14].
C
C     REFERENCES
C
C     [1] Dongarra, J.J., Bunch, J.R., Moler C.B. and Stewart, G.W.
C         LINPACK User's Guide.
C         SIAM, Philadelphia, 1979.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB02BZ by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Bidiagonal matrix, orthogonal transformation, singular values.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           UPDATU, UPDATV
      INTEGER           I, K, LDU, LDV, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), E(*), Q(*), U(LDU,*), V(LDV,*)
C     .. Local Scalars ..
      INTEGER           I1, IROT, L, L1, NROT
      DOUBLE PRECISION  C, F, G, R, S
C     .. External Subroutines ..
      EXTERNAL          DLARTG, DLASR
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     For speed, no tests of the input scalar arguments are done.
C
C     Quick return if possible.
C
      IF ( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
      IF ( I.LE.MIN( M, N ) ) Q(I) = ZERO
C
C     Annihilate E(I) (if I < K).
C
      IF ( I.LT.K ) THEN
         C = ZERO
         S = ONE
         IROT = 0
         NROT = K - I
C
         DO 20 L = I, K-1
            G = E(L)
            E(L) = C*G
            CALL DLARTG( Q(L+1), S*G, C, S, R )
            Q(L+1) = R
            IF ( UPDATU ) THEN
               IROT = IROT + 1
               DWORK(IROT) = C
               DWORK(IROT+NROT) = S
            END IF
   20    CONTINUE
C
         IF ( UPDATU )
     $      CALL DLASR( 'Right', 'Top', 'Forward', M, NROT+1, DWORK(1),
     $                  DWORK(NROT+1), U(1,I), LDU )
      END IF
C
C     Annihilate E(I-1) (if I > 1).
C
      IF ( I.GT.1 ) THEN
         I1 = I - 1
         F  = E(I1)
         E(I1) = ZERO
C
         DO 40 L1 = 1, I1 - 1
            L = I - L1
            CALL DLARTG( Q(L), F, C, S, R )
            Q(L) = R
            IF ( UPDATV ) THEN
               DWORK(L) = C
               DWORK(L+I1) = S
            END IF
            G =  E(L-1)
            F = -S*G
            E(L-1) = C*G
   40    CONTINUE
C
         CALL DLARTG( Q(1), F, C, S, R )
         Q(1) = R
         IF ( UPDATV ) THEN
            DWORK(1) = C
            DWORK(I) = S
            CALL DLASR( 'Right', 'Bottom', 'Backward', N, I, DWORK(1),
     $                  DWORK(I), V(1,1), LDV )
         END IF
      END IF
C
      RETURN
C *** Last line of MB02NY ***
      END
