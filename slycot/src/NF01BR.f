      SUBROUTINE NF01BR( COND, UPLO, TRANS, N, IPAR, LIPAR, R, LDR,
     $                   SDIAG, S, LDS, B, RANKS, TOL, DWORK, LDWORK,
     $                   INFO )
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
C     To solve one of the systems of linear equations
C
C           R*x = b ,  or  R'*x = b ,
C
C     in the least squares sense, where R is an n-by-n block upper
C     triangular matrix, with the structure
C
C         /   R_1    0    ..   0   |   L_1   \
C         |    0    R_2   ..   0   |   L_2   |
C         |    :     :    ..   :   |    :    | ,
C         |    0     0    ..  R_l  |   L_l   |
C         \    0     0    ..   0   |  R_l+1  /
C
C     with the upper triangular submatrices R_k, k = 1:l+1, square, and
C     the first l of the same order, BSN. The diagonal elements of each
C     block R_k have nonincreasing magnitude. The matrix R is stored in
C     the compressed form, as returned by SLICOT Library routine NF01BS,
C
C              /   R_1  |   L_1   \
C              |   R_2  |   L_2   |
C       Rc =   |    :   |    :    | ,
C              |   R_l  |   L_l   |
C              \    X   |  R_l+1  /
C
C     where the submatrix X is irrelevant. If the matrix R does not have
C     full rank, then a least squares solution is obtained. If l <= 1,
C     then R is an upper triangular matrix and its full upper triangle
C     is stored.
C
C     Optionally, the transpose of the matrix R can be stored in the
C     strict lower triangles of the submatrices R_k, k = 1:l+1, and in
C     the arrays SDIAG and S, as described at the parameter UPLO below.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COND    CHARACTER*1
C             Specifies whether the condition of submatrices R_k should
C             be estimated, as follows:
C             = 'E' :  use incremental condition estimation and store
C                      the numerical rank of R_k in the array entry
C                      RANKS(k), for k = 1:l+1;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of R_k for zero values;
C             = 'U' :  use the ranks already stored in RANKS(1:l+1).
C
C     UPLO    CHARACTER*1
C             Specifies the storage scheme for the matrix R, as follows:
C             = 'U' :  the upper triangular part is stored as in Rc;
C             = 'L' :  the lower triangular part is stored, namely,
C                      - the transpose of the strict upper triangle of
C                        R_k is stored in the strict lower triangle of
C                        R_k, for k = 1:l+1;
C                      - the diagonal elements of R_k, k = 1:l+1, are
C                        stored in the array SDIAG;
C                      - the transpose of the last block column in R
C                        (without R_l+1) is stored in the array S.
C
C     TRANS   CHARACTER*1
C             Specifies the form of the system of equations, as follows:
C             = 'N':  R*x  = b  (No transpose);
C             = 'T':  R'*x = b  (Transpose);
C             = 'C':  R'*x = b  (Transpose).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix R.  N = BN*BSN + ST >= 0.
C             (See parameter description below.)
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters describing the structure of the
C             matrix R, as follows:
C             IPAR(1) must contain ST, the number of columns of the
C                     submatrices L_k and the order of R_l+1.  ST >= 0.
C             IPAR(2) must contain BN, the number of blocks, l, in the
C                     block diagonal part of R.  BN >= 0.
C             IPAR(3) must contain BSM, the number of rows of the blocks
C                     R_k, k = 1:l.  BSM >= 0.
C             IPAR(4) must contain BSN, the number of columns of the
C                     blocks R_k, k = 1:l.  BSN >= 0.
C             BSM is not used by this routine, but assumed equal to BSN.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 4.
C
C     R       (input) DOUBLE PRECISION array, dimension (LDR, NC)
C             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1.
C             If UPLO = 'U', the leading N-by-NC part of this array must
C             contain the (compressed) representation (Rc) of the upper
C             triangular matrix R. The submatrix X in Rc and the strict
C             lower triangular parts of the diagonal blocks R_k,
C             k = 1:l+1, are not referenced. If BN <= 1 or BSN = 0, then
C             the full upper triangle of R must be stored.
C             If UPLO = 'L', BN > 1 and BSN > 0, the leading
C             (N-ST)-by-BSN part of this array must contain the
C             transposes of the strict upper triangles of R_k, k = 1:l,
C             stored in the strict lower triangles of R_k, and the
C             strict lower triangle of R_l+1 must contain the transpose
C             of the strict upper triangle of R_l+1. The submatrix X
C             in Rc is not referenced. The diagonal elements of R_k,
C             and, if COND = 'E', the upper triangular parts of R_k,
C             k = 1:l+1, are modified internally, but are restored
C             on exit.
C             If UPLO = 'L' and BN <= 1 or BSN = 0, the leading N-by-N
C             strict lower triangular part of this array must contain
C             the transpose of the strict upper triangular part of R.
C             The diagonal elements and, if COND = 'E', the upper
C             triangular elements are modified internally, but are
C             restored on exit.
C
C     LDR     INTEGER
C             The leading dimension of the array R.  LDR >= MAX(1,N).
C
C     SDIAG   (input) DOUBLE PRECISION array, dimension (N)
C             If UPLO = 'L', this array must contain the diagonal
C             entries of R_k, k = 1:l+1. This array is modified
C             internally, but is restored on exit.
C             This parameter is not referenced if UPLO = 'U'.
C
C     S       (input) DOUBLE PRECISION array, dimension (LDS,N-ST)
C             If UPLO = 'L', BN > 1, and BSN > 0, the leading
C             ST-by-(N-ST) part of this array must contain the transpose
C             of the rectangular part of the last block column in R,
C             that is [ L_1' L_2' ... L_l' ] . If COND = 'E', S is
C             modified internally, but is restored on exit.
C             This parameter is not referenced if UPLO = 'U', or
C             BN <= 1, or BSN = 0.
C
C     LDS     INTEGER
C             The leading dimension of the array S.
C             LDS >= 1,         if UPLO = 'U', or BN <= 1, or BSN = 0;
C             LDS >= MAX(1,ST), if UPLO = 'L', BN > 1, and BSN > 0.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the right hand side
C             vector b.
C             On exit, this array contains the (least squares) solution
C             of the system R*x = b or R'*x = b.
C
C     RANKS   (input or output) INTEGER array, dimension (r), where
C             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1;
C             r = BN,      if ST = 0 and BSN > 0;
C             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 );
C             r = 0,       if ST = 0 and BSN = 0.
C             On entry, if COND = 'U' and N > 0, this array must contain
C             the numerical ranks of the submatrices R_k, k = 1:l(+1).
C             On exit, if COND = 'E' or 'N' and N > 0, this array
C             contains the numerical ranks of the submatrices R_k,
C             k = 1:l(+1), estimated according to the value of COND.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             ranks of the submatrices R_k. If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             the reciprocal condition number;  a (sub)matrix whose
C             estimated condition number is less than 1/TOL is
C             considered to be of full rank. If the user sets TOL <= 0,
C             then an implicitly computed, default tolerance, defined by
C             TOLDEF = N*EPS,  is used instead, where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH).
C             This parameter is not relevant if COND = 'U' or 'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             Denote Full = ( BN <= 1 or  BSN = 0 );
C                    Comp = ( BN >  1 and BSN > 0 ).
C             LDWORK >= 2*N,           if Full and COND = 'E';
C             LDWORK >= 2*MAX(BSN,ST), if Comp and COND = 'E';
C             LDWORK >= 0,   in the remaining cases.
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
C     Block back or forward substitution is used (depending on TRANS
C     and UPLO), exploiting the special structure and storage scheme of
C     the matrix R. If a submatrix R_k, k = 1:l+1, is singular, a local
C     basic least squares solution is computed. Therefore, the returned
C     result is not the basic least squares solution for the whole
C     problem, but a concatenation of (least squares) solutions of the
C     individual subproblems involving R_k, k = 1:l+1 (with adapted
C     right hand sides).
C
C     NUMERICAL ASPECTS
C                                    2    2
C     The algorithm requires 0(BN*BSN + ST + N*ST) operations and is
C     backward stable, if R is nonsingular.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005.
C
C     KEYWORDS
C
C     Linear system of equations, matrix operations, plane rotations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, SVLMAX
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, SVLMAX = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COND, TRANS, UPLO
      INTEGER           INFO, LDR, LDS, LDWORK, LIPAR, N
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IPAR(*), RANKS(*)
      DOUBLE PRECISION  B(*), DWORK(*), R(LDR,*), S(LDS,*), SDIAG(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TOLDEF
      INTEGER           BN, BSM, BSN, I, I1, J, K, L, NC, NTHS, RANK, ST
      CHARACTER         TRANSL, UPLOL
      LOGICAL           ECOND, FULL, LOWER, NCOND, TRANR
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(3)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSWAP, DTRSV, MB03OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      ECOND = LSAME( COND,  'E' )
      NCOND = LSAME( COND,  'N' )
      LOWER = LSAME( UPLO,  'L' )
      TRANR = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
C
      INFO = 0
      IF( .NOT.( ECOND .OR. NCOND .OR. LSAME( COND,  'U' ) ) ) THEN
         INFO = -1
      ELSEIF( .NOT.( LOWER .OR. LSAME( UPLO,  'U' ) ) ) THEN
         INFO = -2
      ELSEIF( .NOT.( TRANR .OR. LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -3
      ELSEIF( N.LT.0 ) THEN
         INFO = -4
      ELSEIF( LIPAR.LT.4 ) THEN
         INFO = -6
      ELSE
         ST   = IPAR(1)
         BN   = IPAR(2)
         BSM  = IPAR(3)
         BSN  = IPAR(4)
         NTHS = BN*BSN
         FULL = BN.LE.1 .OR. BSN.EQ.0
         IF ( MIN( ST, BN, BSM, BSN ).LT.0 ) THEN
            INFO = -5
         ELSEIF ( N.NE.NTHS + ST ) THEN
            INFO = -4
         ELSEIF ( LDR.LT.MAX( 1, N ) ) THEN
            INFO = -8
         ELSEIF ( LDS.LT.1 .OR. ( LOWER .AND. .NOT.FULL .AND.
     $            LDS.LT.ST ) ) THEN
            INFO = -11
         ELSE
            IF ( ECOND ) THEN
               IF ( FULL ) THEN
                  L = 2*N
               ELSE
                  L = 2*MAX( BSN, ST )
               END IF
            ELSE
               L = 0
            END IF
            IF ( LDWORK.LT.L )
     $         INFO = -16
         END IF
      END IF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BR', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      IF ( ECOND ) THEN
         TOLDEF = TOL
         IF ( TOLDEF.LE.ZERO ) THEN
C
C           Use the default tolerance in rank determination.
C
            TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
         END IF
      END IF
C
      NC = BSN + ST
      IF ( FULL ) THEN
C
C        Special case: l <= 1 or BSN = 0; R is just an upper triangular
C        matrix.
C
         IF ( LOWER ) THEN
C
C           Swap the diagonal elements of R and the elements of SDIAG
C           and, if COND = 'E', swap the upper and lower triangular
C           parts of R, in order to find the numerical rank.
C
            CALL DSWAP( N, R, LDR+1, SDIAG, 1 )
            IF ( ECOND ) THEN
               UPLOL  = 'U'
               TRANSL = TRANS
C
               DO 10 J = 1, N
                  CALL DSWAP( N-J+1, R(J,J), LDR, R(J,J), 1 )
   10          CONTINUE
C
            ELSE
               UPLOL = UPLO
               IF ( TRANR ) THEN
                  TRANSL = 'N'
               ELSE
                  TRANSL = 'T'
               END IF
            END IF
         ELSE
            UPLOL  = UPLO
            TRANSL = TRANS
         END IF
C
         IF ( ECOND ) THEN
C
C           Estimate the reciprocal condition number and set the rank.
C           Workspace: 2*N.
C
            CALL MB03OD( 'No QR', N, N, R, LDR, IPAR, TOLDEF, SVLMAX,
     $                   DWORK, RANK, DUM, DWORK, LDWORK, INFO )
            RANKS(1) = RANK
C
         ELSEIF ( NCOND ) THEN
C
C           Determine rank(R) by checking zero diagonal entries.
C
            RANK = N
C
            DO 20 J = 1, N
               IF ( R(J,J).EQ.ZERO .AND. RANK.EQ.N )
     $            RANK = J - 1
   20       CONTINUE
C
            RANKS(1) = RANK
C
         ELSE
C
C           Use the stored rank.
C
            RANK = RANKS(1)
         END IF
C
C        Solve R*x = b, or R'*x = b using back or forward substitution.
C
         DUM(1) = ZERO
         IF ( RANK.LT.N )
     $      CALL DCOPY( N-RANK, DUM, 0, B(RANK+1), 1 )
         CALL DTRSV( UPLOL, TRANSL, 'NonUnit', RANK, R, LDR, B, 1 )
C
         IF ( LOWER ) THEN
C
C           Swap the diagonal elements of R and the elements of SDIAG
C           and, if COND = 'E', swap back the upper and lower triangular
C           parts of R.
C
            CALL DSWAP( N, R, LDR+1, SDIAG, 1 )
            IF ( ECOND ) THEN
C
               DO 30 J = 1, N
                  CALL DSWAP( N-J+1, R(J,J), LDR, R(J,J), 1 )
   30          CONTINUE
C
            END IF
C
         END IF
         RETURN
      END IF
C
C     General case: l > 1 and BSN > 0.
C
      I = 1
      L = BN
      IF ( ECOND ) THEN
C
C        Estimate the reciprocal condition numbers and set the ranks.
C
         IF ( LOWER ) THEN
C
C           Swap the diagonal elements of R and the elements of SDIAG
C           and swap the upper and lower triangular parts of R, in order
C           to find the numerical rank. Swap S and the transpose of the
C           rectangular part of the last block column of R.
C
            DO 50 K = 1, BN
               CALL DSWAP( BSN, R(I,1), LDR+1, SDIAG(I), 1 )
C
               DO 40 J = 1, BSN
                  CALL DSWAP( BSN-J+1, R(I,J), LDR, R(I,J), 1 )
                  I = I + 1
   40          CONTINUE
C
   50       CONTINUE
C
            IF ( ST.GT.0 ) THEN
               CALL DSWAP( ST, R(I,BSN+1), LDR+1, SDIAG(I), 1 )
C
               DO 60 J = BSN + 1, NC
                  CALL DSWAP( NTHS, R(1,J), 1, S(J-BSN,1), LDS )
                  CALL DSWAP( NC-J+1, R(I,J), LDR, R(I,J), 1 )
                  I = I + 1
   60          CONTINUE
C
            END IF
C
         END IF
C
         I1 = 1
C
C        Determine rank(R_k) using incremental condition estimation.
C        Workspace 2*MAX(BSN,ST).
C
         DO 70 K = 1, BN
            CALL MB03OD( 'No QR', BSN, BSN, R(I1,1), LDR, IPAR, TOLDEF,
     $                   SVLMAX, DWORK, RANKS(K), DUM, DWORK, LDWORK,
     $                   INFO )
            I1 = I1 + BSN
   70    CONTINUE
C
         IF ( ST.GT.0 ) THEN
            L = L + 1
            CALL MB03OD( 'No QR', ST, ST, R(I1,BSN+1), LDR, IPAR,
     $                   TOLDEF, SVLMAX, DWORK, RANKS(L), DUM, DWORK,
     $                   LDWORK, INFO )
         END IF
C
      ELSEIF ( NCOND ) THEN
C
C        Determine rank(R_k) by checking zero diagonal entries.
C
         IF ( LOWER ) THEN
C
            DO 90 K = 1, BN
               RANK = BSN
C
               DO 80 J = 1, BSN
                  IF ( SDIAG(I).EQ.ZERO .AND. RANK.EQ.BSN )
     $               RANK = J - 1
                  I = I + 1
   80          CONTINUE
C
               RANKS(K) = RANK
   90       CONTINUE
C
            IF ( ST.GT.0 ) THEN
               L = L + 1
               RANK = ST
C
               DO 100 J = 1, ST
                  IF ( SDIAG(I).EQ.ZERO .AND. RANK.EQ.ST )
     $               RANK = J - 1
                  I = I + 1
  100          CONTINUE
C
               RANKS(L) = RANK
            END IF
C
         ELSE
C
            DO 120 K = 1, BN
               RANK = BSN
C
               DO 110 J = 1, BSN
                  IF ( R(I,J).EQ.ZERO .AND. RANK.EQ.BSN )
     $               RANK = J - 1
                  I = I + 1
  110          CONTINUE
C
               RANKS(K) = RANK
  120       CONTINUE
C
            IF ( ST.GT.0 ) THEN
               L = L + 1
               RANK = ST
C
               DO 130 J = BSN + 1, NC
                  IF ( R(I,J).EQ.ZERO .AND. RANK.EQ.ST )
     $               RANK = J - BSN - 1
                  I = I + 1
  130          CONTINUE
C
               RANKS(L) = RANK
            END IF
         END IF
C
      ELSE
C
C        Set the number of elements of RANKS. Then use the stored ranks.
C
         IF ( ST.GT.0 )
     $      L = L + 1
      END IF
C
C     Solve the triangular system for x. If the system is singular,
C     then obtain a basic least squares solution.
C
      DUM(1) = ZERO
      IF ( LOWER .AND. .NOT.ECOND ) THEN
C
         IF ( .NOT.TRANR ) THEN
C
C           Solve R*x = b using back substitution, with R' stored in
C           the arrays R, SDIAG and S. Swap diag(R) and SDIAG.
C
            I1 = NTHS + 1
            IF ( ST.GT.0 ) THEN
               RANK = RANKS(L)
               IF ( RANK.LT.ST )
     $            CALL DCOPY( ST-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DSWAP( ST, R(I1,BSN+1), LDR+1, SDIAG(I1), 1 )
               CALL DTRSV( 'Lower', 'Transpose', 'NonUnit', RANK,
     $                     R(I1,BSN+1), LDR, B(I1), 1 )
               CALL DSWAP( ST, R(I1,BSN+1), LDR+1, SDIAG(I1), 1 )
               CALL DGEMV( 'Transpose', ST, NTHS, -ONE, S, LDS,
     $                     B(NTHS+1), 1, ONE, B, 1 )
            END IF
C
            DO 140 K = BN, 1, -1
               I1   = I1 - BSN
               RANK = RANKS(K)
               IF ( RANK.LT.BSN )
     $            CALL DCOPY( BSN-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DSWAP( BSN, R(I1,1), LDR+1, SDIAG(I1), 1 )
               CALL DTRSV( 'Lower', 'Transpose', 'NonUnit', RANK,
     $                     R(I1,1), LDR, B(I1), 1 )
               CALL DSWAP( BSN, R(I1,1), LDR+1, SDIAG(I1), 1 )
  140       CONTINUE
C
         ELSE
C
C           Solve R'*x = b using forward substitution, with R' stored in
C           the arrays R, SDIAG and S. Swap diag(R) and SDIAG.
C
            I1 = 1
            IF ( TRANR ) THEN
               TRANSL = 'N'
            ELSE
               TRANSL = 'T'
            END IF
C
            DO 150 K = 1, BN
               RANK = RANKS(K)
               IF ( RANK.LT.BSN )
     $            CALL DCOPY( BSN-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DSWAP( BSN, R(I1,1), LDR+1, SDIAG(I1), 1 )
               CALL DTRSV( 'Lower', TRANSL, 'NonUnit', RANK, R(I1,1),
     $                     LDR, B(I1), 1 )
               CALL DSWAP( BSN, R(I1,1), LDR+1, SDIAG(I1), 1 )
               I1 = I1 + BSN
  150       CONTINUE
C
            IF ( ST.GT.0 ) THEN
               RANK = RANKS(L)
               IF ( RANK.LT.ST )
     $            CALL DCOPY( ST-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DGEMV( 'NoTranspose', ST, NTHS, -ONE, S, LDS, B, 1,
     $                     ONE, B(I1), 1 )
               CALL DSWAP( ST, R(I1,BSN+1), LDR+1, SDIAG(I1), 1 )
               CALL DTRSV( 'Lower', TRANSL, 'NonUnit', RANK,
     $                     R(I1,BSN+1), LDR, B(I1), 1 )
               CALL DSWAP( ST, R(I1,BSN+1), LDR+1, SDIAG(I1), 1 )
            END IF
C
         END IF
C
      ELSE
C
         IF ( .NOT.TRANR ) THEN
C
C           Solve R*x = b using back substitution.
C
            I1 = NTHS + 1
            IF ( ST.GT.0 ) THEN
               RANK = RANKS(L)
               IF ( RANK.LT.ST )
     $            CALL DCOPY( ST-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DTRSV( 'Upper', TRANS, 'NonUnit', RANK, R(I1,BSN+1),
     $                     LDR, B(I1), 1 )
               CALL DGEMV( TRANS, NTHS, ST, -ONE, R(1,BSN+1), LDR,
     $                     B(NTHS+1), 1, ONE, B, 1 )
            END IF
C
            DO 160 K = BN, 1, -1
               I1   = I1 - BSN
               RANK = RANKS(K)
               IF ( RANK.LT.BSN )
     $            CALL DCOPY( BSN-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DTRSV( 'Upper', TRANS, 'NonUnit', RANK, R(I1,1),
     $                     LDR, B(I1), 1 )
  160       CONTINUE
C
         ELSE
C
C           Solve R'*x = b using forward substitution.
C
            I1 = 1
C
            DO 170 K = 1, BN
               RANK = RANKS(K)
               IF ( RANK.LT.BSN )
     $            CALL DCOPY( BSN-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DTRSV( 'Upper', TRANS, 'NonUnit', RANK, R(I1,1),
     $                     LDR, B(I1), 1 )
               I1 = I1 + BSN
  170       CONTINUE
C
            IF ( ST.GT.0 ) THEN
               RANK = RANKS(L)
               IF ( RANK.LT.ST )
     $            CALL DCOPY( ST-RANK, DUM, 0, B(I1+RANK), 1 )
               CALL DGEMV( TRANS, NTHS, ST, -ONE, R(1,BSN+1), LDR, B, 1,
     $                     ONE, B(I1), 1 )
               CALL DTRSV( 'Upper', TRANS, 'NonUnit', RANK, R(I1,BSN+1),
     $                     LDR, B(I1), 1 )
            END IF
C
         END IF
      END IF
C
      IF ( ECOND .AND. LOWER ) THEN
         I = 1
C
C        If COND = 'E' and UPLO = 'L', swap the diagonal elements of R
C        and the elements of SDIAG and swap back the upper and lower
C        triangular parts of R, including the part corresponding to S.
C
         DO 190 K = 1, BN
            CALL DSWAP( BSN, R(I,1), LDR+1, SDIAG(I), 1 )
C
            DO 180 J = 1, BSN
               CALL DSWAP( BSN-J+1, R(I,J), LDR, R(I,J), 1 )
               I = I + 1
  180       CONTINUE
C
  190    CONTINUE
C
         IF ( ST.GT.0 ) THEN
            CALL DSWAP( ST, R(I,BSN+1), LDR+1, SDIAG(I), 1 )
C
            DO 200 J = BSN + 1, NC
               CALL DSWAP( NTHS, R(1,J), 1, S(J-BSN,1), LDS )
               CALL DSWAP( NC-J+1, R(I,J), LDR, R(I,J), 1 )
               I = I + 1
  200       CONTINUE
C
         END IF
C
      END IF
C
      RETURN
C
C *** Last line of NF01BR ***
      END
