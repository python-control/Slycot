      SUBROUTINE NF01BQ( COND, N, IPAR, LIPAR, R, LDR, IPVT, DIAG, QTB,
     $                   RANKS, X, TOL, DWORK, LDWORK, INFO )
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
C     To determine a vector x which solves the system of linear
C     equations
C
C           J*x = b ,     D*x = 0 ,
C
C     in the least squares sense, where J is an m-by-n matrix,
C     D is an n-by-n diagonal matrix, and b is an m-vector. The matrix J
C     is the current Jacobian of a nonlinear least squares problem,
C     provided in a compressed form by SLICOT Library routine NF01BD.
C     It is assumed that a block QR factorization, with column pivoting,
C     of J is available, that is, J*P = Q*R, where P is a permutation
C     matrix, Q has orthogonal columns, and R is an upper triangular
C     matrix with diagonal elements of nonincreasing magnitude for each
C     block, as returned by SLICOT Library routine NF01BS. The routine
C     NF01BQ needs the upper triangle of R in compressed form, the
C     permutation matrix P, and the first n components of Q'*b
C     (' denotes the transpose). The system J*x = b, D*x = 0, is then
C     equivalent to
C
C           R*z = Q'*b ,  P'*D*P*z = 0 ,                             (1)
C
C     where x = P*z. If this system does not have full rank, then an
C     approximate least squares solution is obtained (see METHOD).
C     On output, NF01BQ also provides an upper triangular matrix S
C     such that
C
C           P'*(J'*J + D*D)*P = S'*S .
C
C     The system (1) is equivalent to S*z = c , where c contains the
C     first n components of the vector obtained by applying to
C     [ (Q'*b)'  0 ]' the transformations which triangularized
C     [ R'  P'*D*P ]', getting S.
C
C     The matrix R has the following structure
C
C         /   R_1    0    ..   0   |   L_1   \
C         |    0    R_2   ..   0   |   L_2   |
C         |    :     :    ..   :   |    :    | ,
C         |    0     0    ..  R_l  |   L_l   |
C         \    0     0    ..   0   |  R_l+1  /
C
C     where the submatrices R_k, k = 1:l, have the same order BSN,
C     and R_k, k = 1:l+1, are square and upper triangular. This matrix
C     is stored in the compressed form
C
C              /   R_1  |   L_1   \
C              |   R_2  |   L_2   |
C       Rc =   |    :   |    :    | ,
C              |   R_l  |   L_l   |
C              \    X   |  R_l+1  /
C
C     where the submatrix X is irrelevant. The matrix S has the same
C     structure as R, and its diagonal blocks are denoted by S_k,
C     k = 1:l+1.
C
C     If l <= 1, then the full upper triangle of the matrix R is stored.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COND    CHARACTER*1
C             Specifies whether the condition of the matrices S_k should
C             be estimated, as follows:
C             = 'E' :  use incremental condition estimation and store
C                      the numerical rank of S_k in the array entry
C                      RANKS(k), for k = 1:l+1;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of S_k for zero values;
C             = 'U' :  use the ranks already stored in RANKS(1:l+1).
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
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR, NC)
C             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1.
C             On entry, the leading N-by-NC part of this array must
C             contain the (compressed) representation (Rc) of the upper
C             triangular matrix R. If BN > 1, the submatrix X in Rc is
C             not referenced. The zero strict lower triangles of R_k,
C             k = 1:l+1, need not be set. If BN <= 1 or BSN = 0, then
C             the full upper triangle of R must be stored.
C             On exit, the full upper triangles of R_k, k = 1:l+1, and
C             L_k, k = 1:l, are unaltered, and the strict lower
C             triangles of R_k, k = 1:l+1, contain the corresponding
C             strict upper triangles (transposed) of the upper
C             triangular matrix S.
C             If BN <= 1 or BSN = 0, then the transpose of the strict
C             upper triangle of S is stored in the strict lower triangle
C             of R.
C
C     LDR     INTEGER
C             The leading dimension of the array R.  LDR >= MAX(1,N).
C
C     IPVT    (input) INTEGER array, dimension (N)
C             This array must define the permutation matrix P such that
C             J*P = Q*R. Column j of P is column IPVT(j) of the identity
C             matrix.
C
C     DIAG    (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the diagonal elements of the
C             matrix D.
C
C     QTB     (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the first n elements of the
C             vector Q'*b.
C
C     RANKS   (input or output) INTEGER array, dimension (r), where
C             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1;
C             r = BN,      if ST = 0 and BSN > 0;
C             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 );
C             r = 0,       if ST = 0 and BSN = 0.
C             On entry, if COND = 'U' and N > 0, this array must contain
C             the numerical ranks of the submatrices S_k, k = 1:l(+1).
C             On exit, if COND = 'E' or 'N' and N > 0, this array
C             contains the numerical ranks of the submatrices S_k,
C             k = 1:l(+1), estimated according to the value of COND.
C
C     X       (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the least squares solution of the
C             system J*x = b, D*x = 0.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             ranks of the submatrices S_k. If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             the reciprocal condition number;  a (sub)matrix whose
C             estimated condition number is less than 1/TOL is
C             considered to be of full rank.  If the user sets TOL <= 0,
C             then an implicitly computed, default tolerance, defined by
C             TOLDEF = N*EPS,  is used instead, where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH).
C             This parameter is not relevant if COND = 'U' or 'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, the first N elements of this array contain the
C             diagonal elements of the upper triangular matrix S, and
C             the next N elements contain the solution z.
C             If BN > 1 and BSN > 0, the elements 2*N+1 : 2*N+ST*(N-ST)
C             contain the submatrix (S(1:N-ST,N-ST+1:N))' of the
C             matrix S.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 2*N,              if BN <= 1 or  BSN = 0 and
C                                                        COND <> 'E';
C             LDWORK >= 4*N,              if BN <= 1 or  BSN = 0 and
C                                                        COND =  'E';
C             LDWORK >= ST*(N-ST) + 2*N,  if BN >  1 and BSN > 0 and
C                                                        COND <> 'E';
C             LDWORK >= ST*(N-ST) + 2*N + 2*MAX(BSN,ST),
C                                         if BN >  1 and BSN > 0 and
C                                                        COND = 'E'.
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
C     Standard plane rotations are used to annihilate the elements of
C     the diagonal matrix D, updating the upper triangular matrix R
C     and the first n elements of the vector Q'*b. A basic least squares
C     solution is computed. The computations exploit the special
C     structure and storage scheme of the matrix R. If one or more of
C     the submatrices S_k, k = 1:l+1, is singular, then the computed
C     result is not the basic least squares solution for the whole
C     problem, but a concatenation of (least squares) solutions of the
C     individual subproblems involving R_k, k = 1:l+1 (with adapted
C     right hand sides).
C
C     REFERENCES
C
C     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E.
C         User's Guide for MINPACK-1.
C         Applied Math. Division, Argonne National Laboratory, Argonne,
C         Illinois, Report ANL-80-74, 1980.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0(N*(BSN+ST)) operations and is backward
C     stable, if R is nonsingular.
C
C     FURTHER COMMENTS
C
C     This routine is a structure-exploiting, LAPACK-based modification
C     of QRSOLV from the MINPACK package [1], and with optional
C     condition estimation.
C     The option COND = 'U' is useful when dealing with several
C     right-hand side vectors.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Linear system of equations, matrix operations, plane rotations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COND
      INTEGER           INFO, LDR, LDWORK, LIPAR, N
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IPAR(*), IPVT(*), RANKS(*)
      DOUBLE PRECISION  DIAG(*), DWORK(*), QTB(*), R(LDR,*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  QTBPJ
      INTEGER           BN, BSM, BSN, I, IB, IBSN, IS, ITC, ITR, J,
     $                  JW, K, KF, L, NC, NTHS, ST
      LOGICAL           ECOND
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, MB02YD, MB04OW, NF01BR, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      ECOND = LSAME( COND, 'E' )
      INFO  = 0
      IF( .NOT.( ECOND .OR. LSAME( COND, 'N' ) .OR.
     $                      LSAME( COND, 'U' ) ) ) THEN
         INFO = -1
      ELSEIF( N.LT.0 ) THEN
         INFO = -2
      ELSEIF( LIPAR.LT.4 ) THEN
         INFO = -4
      ELSE
         ST   = IPAR(1)
         BN   = IPAR(2)
         BSM  = IPAR(3)
         BSN  = IPAR(4)
         NTHS = BN*BSN
         IF ( MIN( ST, BN, BSM, BSN ).LT.0 ) THEN
            INFO = -3
         ELSEIF ( N.NE.NTHS + ST ) THEN
            INFO = -2
         ELSEIF ( LDR.LT.MAX( 1, N ) ) THEN
            INFO = -6
         ELSE
            JW = 2*N
            IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
               IF ( ECOND )
     $            JW = 4*N
            ELSE
               JW = ST*NTHS + JW
               IF ( ECOND )
     $            JW = 2*MAX( BSN, ST ) + JW
            END IF
            IF ( LDWORK.LT.JW )
     $         INFO = -14
         ENDIF
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BQ', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
C
C        Special case: R is an upper triangular matrix.
C        Workspace: 4*N, if COND =  'E';
C                   2*N, if COND <> 'E'.
C
         CALL MB02YD( COND, N, R, LDR, IPVT, DIAG, QTB, RANKS(1), X,
     $                TOL, DWORK, LDWORK, INFO )
         RETURN
      END IF
C
C     General case: BN > 1 and BSN > 0.
C     Copy R and Q'*b to preserve input and initialize S.
C     In particular, save the diagonal elements of R in X.
C
      IB = N  + 1
      IS = IB + N
      JW = IS + ST*NTHS
      I  = 1
      L  = IS
      NC = BSN + ST
      KF = NC
C
      DO 20 K = 1, BN
C
         DO 10 J = 1, BSN
            X(I) = R(I,J)
            CALL DCOPY( BSN-J+1, R(I,J), LDR, R(I,J), 1 )
            I = I + 1
   10    CONTINUE
C
   20 CONTINUE
C
C     DWORK(IS) contains a copy of [ L_1' ... L_l' ].
C     Workspace:  ST*(N-ST)+2*N;
C
      DO 30 J = BSN + 1, NC
         CALL DCOPY( NTHS, R(1,J), 1, DWORK(L), ST )
         X(I) = R(I,J)
         CALL DCOPY( NC-J+1, R(I,J), LDR, R(I,J), 1 )
         I = I + 1
         L = L + 1
   30 CONTINUE
C
      CALL DCOPY( N, QTB, 1, DWORK(IB), 1 )
      IF ( ST.GT.0 ) THEN
         ITR = NTHS + 1
         ITC = BSN  + 1
      ELSE
         ITR = 1
         ITC = 1
      END IF
      IBSN = 0
C
C     Eliminate the diagonal matrix D using Givens rotations.
C
      DO 50 J = 1, N
         IBSN = IBSN + 1
         I    = IBSN
C
C        Prepare the row of D to be eliminated, locating the
C        diagonal element using P from the QR factorization.
C
         L = IPVT(J)
         IF ( DIAG(L).NE.ZERO ) THEN
            QTBPJ = ZERO
            DWORK(J) = DIAG(L)
C
            DO 40 K = J + 1, MIN( J + KF - 1, N )
               DWORK(K) = ZERO
   40       CONTINUE
C
C           The transformations to eliminate the row of D modify only
C           a single element of Q'*b beyond the first n, which is
C           initially zero.
C
            IF ( J.LT.NTHS ) THEN
               CALL MB04OW( BSN-IBSN+1, ST, 1, R(J,IBSN), LDR,
     $                      R(ITR,ITC), LDR, DWORK(J), 1, DWORK(IB+J-1),
     $                      BSN, DWORK(IB+NTHS), ST, QTBPJ, 1 )
               IF ( IBSN.EQ.BSN )
     $            IBSN = 0
            ELSE IF ( J.EQ.NTHS ) THEN
               CALL MB04OW( 1, ST, 1, R(J,IBSN), LDR, R(ITR,ITC), LDR,
     $                      DWORK(J), 1, DWORK(IB+J-1), BSN,
     $                      DWORK(IB+NTHS), ST, QTBPJ, 1 )
               KF = ST
            ELSE
               CALL MB04OW( 0, N-J+1, 1, R(J,IBSN), LDR, R(J,IBSN), LDR,
     $                      DWORK(J), 1, DWORK(IB+J-1), 1,
     $                      DWORK(IB+J-1), ST, QTBPJ, 1 )
            END IF
         ELSE
            IF ( J.LT.NTHS ) THEN
               IF ( IBSN.EQ.BSN )
     $            IBSN = 0
               ELSE IF ( J.EQ.NTHS ) THEN
                  KF = ST
            END IF
         END IF
C
C        Store the diagonal element of S.
C
         DWORK(J) = R(J,I)
   50 CONTINUE
C
C     Solve the triangular system for z. If the system is singular,
C     then obtain an approximate least squares solution.
C     Additional workspace:   2*MAX(BSN,ST), if COND =  'E';
C                             0,             if COND <> 'E'.
C
      CALL NF01BR( COND, 'Upper', 'NoTranspose', N, IPAR, LIPAR, R, LDR,
     $             DWORK, DWORK(IS), 1, DWORK(IB), RANKS, TOL,
     $             DWORK(JW), LDWORK-JW+1, INFO )
      I = 1
C
C     Restore the diagonal elements of R from X and interchange
C     the upper and lower triangular parts of R.
C
      DO 70 K = 1, BN
C
         DO 60 J = 1, BSN
            R(I,J) = X(I)
            CALL DSWAP( BSN-J+1, R(I,J), LDR, R(I,J), 1 )
            I = I + 1
   60    CONTINUE
C
   70 CONTINUE
C
      DO 80 J = BSN + 1, NC
         CALL DSWAP( NTHS, R(1,J), 1, DWORK(IS), ST )
         R(I,J) = X(I)
         CALL DSWAP( NC-J+1, R(I,J), LDR, R(I,J), 1 )
         I  = I  + 1
         IS = IS + 1
   80 CONTINUE
C
C     Permute the components of z back to components of x.
C
      DO 90 J = 1, N
         L    = IPVT(J)
         X(L) = DWORK(N+J)
   90 CONTINUE
C
      RETURN
C
C *** Last line of NF01BQ ***
      END
