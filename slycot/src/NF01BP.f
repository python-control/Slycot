      SUBROUTINE NF01BP( COND, N, IPAR, LIPAR, R, LDR, IPVT, DIAG, QTB,
     $                   DELTA, PAR, RANKS, X, RX, TOL, DWORK, LDWORK,
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
C     To determine a value for the Levenberg-Marquardt parameter PAR
C     such that if x solves the system
C
C           J*x = b ,     sqrt(PAR)*D*x = 0 ,
C
C     in the least squares sense, where J is an m-by-n matrix, D is an
C     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if
C     DELTA is a positive number, DXNORM is the Euclidean norm of D*x,
C     then either PAR is zero and
C
C           ( DXNORM - DELTA ) .LE. 0.1*DELTA ,
C
C     or PAR is positive and
C
C           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA .
C
C     The matrix J is the current Jacobian matrix of a nonlinear least
C     squares problem, provided in a compressed form by SLICOT Library
C     routine NF01BD. It is assumed that a block QR factorization, with
C     column pivoting, of J is available, that is, J*P = Q*R, where P is
C     a permutation matrix, Q has orthogonal columns, and R is an upper
C     triangular matrix with diagonal elements of nonincreasing
C     magnitude for each block, as returned by SLICOT Library
C     routine NF01BS. The routine NF01BP needs the upper triangle of R
C     in compressed form, the permutation matrix P, and the first
C     n components of Q'*b (' denotes the transpose). On output,
C     NF01BP also provides a compressed representation of an upper
C     triangular matrix S, such that
C
C           P'*(J'*J + PAR*D*D)*P = S'*S .
C
C     Matrix S is used in the solution process. The matrix R has the
C     following structure
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
C             Specifies whether the condition of the diagonal blocks R_k
C             and S_k of the matrices R and S should be estimated,
C             as follows:
C             = 'E' :  use incremental condition estimation for each
C                      diagonal block of R_k and S_k to find its
C                      numerical rank;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of R_k and S_k for zero values;
C             = 'U' :  use the ranks already stored in RANKS (for R).
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
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     IPVT    (input) INTEGER array, dimension (N)
C             This array must define the permutation matrix P such that
C             J*P = Q*R. Column j of P is column IPVT(j) of the identity
C             matrix.
C
C     DIAG    (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the diagonal elements of the
C             matrix D.  DIAG(I) <> 0, I = 1,...,N.
C
C     QTB     (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the first n elements of the
C             vector Q'*b.
C
C     DELTA   (input) DOUBLE PRECISION
C             An upper bound on the Euclidean norm of D*x.  DELTA > 0.
C
C     PAR     (input/output) DOUBLE PRECISION
C             On entry, PAR must contain an initial estimate of the
C             Levenberg-Marquardt parameter.  PAR >= 0.
C             On exit, it contains the final estimate of this parameter.
C
C     RANKS   (input or output) INTEGER array, dimension (r), where
C             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1;
C             r = BN,      if ST = 0 and BSN > 0;
C             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 );
C             r = 0,       if ST = 0 and BSN = 0.
C             On entry, if COND = 'U' and N > 0, this array must contain
C             the numerical ranks of the submatrices R_k, k = 1:l(+1).
C             On exit, if N > 0, this array contains the numerical ranks
C             of the submatrices S_k, k = 1:l(+1).
C
C     X       (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the least squares solution of the
C             system J*x = b, sqrt(PAR)*D*x = 0.
C
C     RX      (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the matrix-vector product -R*P'*x.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             ranks of the submatrices R_k and S_k. If the user sets
C             TOL > 0, then the given value of TOL is used as a lower
C             bound for the reciprocal condition number;  a (sub)matrix
C             whose estimated condition number is less than 1/TOL is
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
C             diagonal elements of the upper triangular matrix S.
C             If BN > 1 and BSN > 0, the elements N+1 : N+ST*(N-ST)
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
C                                                        COND =  'E'.
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
C     The algorithm computes the Gauss-Newton direction. An approximate
C     basic least squares solution is found if the Jacobian is rank
C     deficient. The computations exploit the special structure and
C     storage scheme of the matrix R. If one or more of the submatrices
C     R_k or S_k, k = 1:l+1, is singular, then the computed result is
C     not the basic least squares solution for the whole problem, but a
C     concatenation of (least squares) solutions of the individual
C     subproblems involving R_k or S_k, k = 1:l+1 (with adapted right
C     hand sides).
C
C     If the Gauss-Newton direction is not acceptable, then an iterative
C     algorithm obtains improved lower and upper bounds for the
C     Levenberg-Marquardt parameter PAR. Only a few iterations are
C     generally needed for convergence of the algorithm. If, however,
C     the limit of ITMAX = 10 iterations is reached, then the output PAR
C     will contain the best value obtained so far. If the Gauss-Newton
C     step is acceptable, it is stored in x, and PAR is set to zero,
C     hence S = R.
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
C     of LMPAR from the MINPACK package [1], and with optional condition
C     estimation. The option COND = 'U' is useful when dealing with
C     several right-hand side vectors, but RANKS array should be reset.
C     If COND = 'E', but the matrix S is guaranteed to be nonsingular
C     and well conditioned relative to TOL, i.e., rank(R) = N, and
C     min(DIAG) > 0, then its condition is not estimated.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     V. Sima, Feb. 2004.
C
C     KEYWORDS
C
C     Linear system of equations, matrix operations, plane rotations.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           ITMAX
      PARAMETER         ( ITMAX = 10 )
      DOUBLE PRECISION  P1, P001, ZERO, ONE
      PARAMETER         ( P1  = 1.0D-1, P001 = 1.0D-3, ZERO = 0.0D0,
     $                    ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COND
      INTEGER           INFO, LDR, LDWORK, LIPAR, N
      DOUBLE PRECISION  DELTA, PAR, TOL
C     .. Array Arguments ..
      INTEGER           IPAR(*), IPVT(*), RANKS(*)
      DOUBLE PRECISION  DIAG(*), DWORK(*), QTB(*), R(LDR,*), RX(*), X(*)
C     .. Local Scalars ..
      INTEGER           BN, BSM, BSN, I, IBSN, ITER, J, JW, K, L, LDS,
     $                  N2, NTHS, RANK, ST
      DOUBLE PRECISION  DMINO, DWARF, DXNORM, FP, GNORM, PARC, PARL,
     $                  PARU, SUM, TEMP, TOLDEF
      LOGICAL           BADRK, ECOND, NCOND, SING, UCOND
      CHARACTER         CONDL
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DDOT, DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DTRMV, MD03BY, NF01BQ, NF01BR,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      ECOND = LSAME( COND, 'E' )
      NCOND = LSAME( COND, 'N' )
      UCOND = LSAME( COND, 'U' )
      INFO  = 0
      N2    = 2*N
      IF( .NOT.( ECOND .OR. NCOND .OR. UCOND ) ) THEN
         INFO = -1
      ELSEIF( N.LT.0 ) THEN
         INFO = -2
      ELSEIF( LIPAR.LT.4 ) THEN
         INFO = -4
      ELSEIF ( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSEIF( DELTA.LE.ZERO ) THEN
         INFO = -10
      ELSEIF( PAR.LT.ZERO ) THEN
         INFO = -11
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
         ELSE
            IF ( N.GT.0 )
     $         DMINO = DIAG(1)
            SING = .FALSE.
C
            DO 10 J = 1, N
               IF ( DIAG(J).LT.DMINO )
     $            DMINO = DIAG(J)
               SING = SING .OR. DIAG(J).EQ.ZERO
   10       CONTINUE
C
            IF ( SING ) THEN
               INFO = -8
            ELSEIF ( UCOND ) THEN
               BADRK = .FALSE.
               IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
                  IF ( N.GT.0 )
     $               BADRK = RANKS(1).LT.0 .OR. RANKS(1).GT.N
               ELSE
                  RANK  = 0
C
                  DO 20 K = 1, BN
                     BADRK = BADRK .OR. RANKS(K).LT.0
     $                             .OR. RANKS(K).GT.BSN
                     RANK  = RANK + RANKS(K)
   20             CONTINUE
C
                  IF ( ST.GT.0 ) THEN
                     BADRK = BADRK .OR. RANKS(BN+1).LT.0 .OR.
     $                                  RANKS(BN+1).GT.ST
                     RANK  = RANK + RANKS(BN+1)
                  END IF
               END IF
               IF ( BADRK )
     $            INFO = -12
            ELSE
               JW = N2
               IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
                  IF ( ECOND )
     $               JW = 4*N
               ELSE
                  JW = ST*NTHS + JW
                  IF ( ECOND )
     $               JW = 2*MAX( BSN, ST ) + JW
               END IF
               IF ( LDWORK.LT.JW )
     $            INFO = -17
            ENDIF
         ENDIF
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BP', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         PAR = ZERO
         RETURN
      END IF
C
      IF ( BN.LE.1 .OR. BSN.EQ.0 ) THEN
C
C        Special case: R is just an upper triangular matrix.
C        Workspace: 4*N, if COND =  'E';
C                   2*N, if COND <> 'E'.
C
         CALL MD03BY( COND, N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR,
     $                RANKS(1), X, RX, TOL, DWORK, LDWORK, INFO )
         RETURN
      END IF
C
C     General case: l > 1 and BSN > 0.
C     DWARF is the smallest positive magnitude.
C
      DWARF = DLAMCH( 'Underflow' )
C
C     Compute and store in x the Gauss-Newton direction. If the
C     Jacobian is rank-deficient, obtain a least squares solution.
C     The array RX is used as workspace.
C     Workspace: 2*MAX(BSN,ST), if COND =  'E';
C                0,             if COND <> 'E'.
C
      CALL DCOPY( N, QTB, 1, RX, 1 )
      CALL NF01BR( COND, 'Upper', 'No transpose', N, IPAR, LIPAR, R,
     $             LDR, DWORK, DWORK, 1, RX, RANKS, TOL, DWORK, LDWORK,
     $             INFO )
C
      DO 30 J = 1, N
         L    = IPVT(J)
         X(L) = RX(J)
   30 CONTINUE
C
C     Initialize the iteration counter.
C     Evaluate the function at the origin, and test
C     for acceptance of the Gauss-Newton direction.
C
      ITER = 0
C
      DO 40 J = 1, N
         DWORK(J) = DIAG(J)*X(J)
   40 CONTINUE
C
      DXNORM = DNRM2( N, DWORK, 1 )
      FP     = DXNORM - DELTA
      IF ( FP.GT.P1*DELTA ) THEN
C
C        Set an appropriate option for estimating the condition of
C        the matrix S.
C
         LDS = MAX( 1, ST )
         JW  = N2 + ST*NTHS
         IF ( UCOND ) THEN
            IF ( LDWORK.GE.JW + 2*MAX( BSN, ST ) ) THEN
               CONDL  = 'E'
               TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
            ELSE
               CONDL  = 'N'
               TOLDEF = TOL
            END IF
         ELSE
            RANK = 0
C
            DO 50 K = 1, BN
               RANK = RANK + RANKS(K)
   50       CONTINUE
C
            IF ( ST.GT.0 )
     $         RANK = RANK + RANKS(BN+1)
            CONDL  = COND
            TOLDEF = TOL
         END IF
C
C        If the Jacobian is not rank deficient, the Newton
C        step provides a lower bound, PARL, for the zero of
C        the function. Otherwise set this bound to zero.
C
         IF ( RANK.EQ.N ) THEN
C
            DO 60 J = 1, N
               L = IPVT(J)
               RX(J) = DIAG(L)*( DWORK(L)/DXNORM )
   60       CONTINUE
C
            CALL NF01BR( 'Use ranks', 'Upper', 'Transpose', N, IPAR,
     $                   LIPAR, R, LDR, DWORK, DWORK, 1, RX, RANKS, TOL,
     $                   DWORK, LDWORK, INFO )
            TEMP = DNRM2( N, RX, 1 )
            PARL = ( ( FP/DELTA )/TEMP )/TEMP
C
C           For efficiency, use CONDL = 'U', if possible.
C
            IF ( .NOT.LSAME( CONDL, 'U' ) .AND. DMINO.GT.ZERO )
     $         CONDL = 'U'
         ELSE
            PARL = ZERO
         END IF
C
         IBSN = 0
         K    = 1
C
C        Calculate an upper bound, PARU, for the zero of the function.
C
         DO 70 J = 1, N
            IBSN = IBSN + 1
            IF ( J.LT.NTHS ) THEN
               SUM = DDOT( IBSN, R(K,IBSN), 1, QTB(K), 1 )
               IF ( IBSN.EQ.BSN ) THEN
                  IBSN = 0
                  K    = K + BSN
               END IF
            ELSE IF ( J.EQ.NTHS ) THEN
               SUM = DDOT( IBSN, R(K,IBSN), 1, QTB(K), 1 )
            ELSE
               SUM = DDOT( J, R(1,IBSN), 1, QTB, 1 )
            END IF
            L = IPVT(J)
            RX(J) = SUM/DIAG(L)
   70    CONTINUE
C
         GNORM = DNRM2( N, RX, 1 )
         PARU  = GNORM/DELTA
         IF ( PARU.EQ.ZERO )
     $      PARU = DWARF/MIN( DELTA, P1 )/P001
C
C        If the input PAR lies outside of the interval (PARL,PARU),
C        set PAR to the closer endpoint.
C
         PAR = MAX( PAR, PARL )
         PAR = MIN( PAR, PARU )
         IF ( PAR.EQ.ZERO )
     $      PAR = GNORM/DXNORM
C
C        Beginning of an iteration.
C
   80    CONTINUE
            ITER = ITER + 1
C
C           Evaluate the function at the current value of PAR.
C
            IF ( PAR.EQ.ZERO )
     $         PAR = MAX( DWARF, P001*PARU )
            TEMP = SQRT( PAR )
C
            DO 90 J = 1, N
               RX(J) = TEMP*DIAG(J)
   90       CONTINUE
C
C           Solve the system J*x = b , sqrt(PAR)*D*x = 0 , in a least
C           square sense.
C           The first N elements of DWORK contain the diagonal elements
C           of the upper triangular matrix S, and the next N elements
C           contain the the vector z, so that x = P*z (see NF01BQ).
C           The vector z is not preserved, to reduce the workspace.
C           The elements 2*N+1 : 2*N+ST*(N-ST) contain the
C           submatrix (S(1:N-ST,N-ST+1:N))' of the matrix S.
C           Workspace: ST*(N-ST) + 2*N,                 if CONDL <> 'E';
C                      ST*(N-ST) + 2*N + 2*MAX(BSN,ST), if CONDL =  'E'.
C
            CALL NF01BQ( CONDL, N, IPAR, LIPAR, R, LDR, IPVT, RX, QTB,
     $                   RANKS, X, TOLDEF, DWORK, LDWORK, INFO )
C
            DO 100 J = 1, N
               DWORK(N+J) = DIAG(J)*X(J)
  100       CONTINUE
C
            DXNORM = DNRM2( N, DWORK(N+1), 1 )
            TEMP   = FP
            FP     = DXNORM - DELTA
C
C           If the function is small enough, accept the current value
C           of PAR. Also test for the exceptional cases where PARL
C           is zero or the number of iterations has reached ITMAX.
C
            IF ( ABS( FP ).GT.P1*DELTA .AND.
     $         ( PARL.NE.ZERO .OR. FP.GT.TEMP .OR. TEMP.GE.ZERO ) .AND.
     $           ITER.LT.ITMAX ) THEN
C
C              Compute the Newton correction.
C
               DO 110 J = 1, N
                  L = IPVT(J)
                  RX(J) = DIAG(L)*( DWORK(N+L)/DXNORM )
  110          CONTINUE
C
               CALL NF01BR( 'Use ranks', 'Lower', 'Transpose', N, IPAR,
     $                      LIPAR, R, LDR, DWORK, DWORK(N2+1), LDS, RX,
     $                      RANKS, TOL, DWORK(JW), LDWORK-JW, INFO )
               TEMP = DNRM2( N, RX, 1 )
               PARC = ( ( FP/DELTA )/TEMP )/TEMP
C
C              Depending on the sign of the function, update PARL
C              or PARU.
C
               IF ( FP.GT.ZERO ) THEN
                  PARL = MAX( PARL, PAR )
               ELSE IF ( FP.LT.ZERO ) THEN
                  PARU = MIN( PARU, PAR )
               END IF
C
C              Compute an improved estimate for PAR.
C
               PAR = MAX( PARL, PAR + PARC )
C
C              End of an iteration.
C
               GO TO 80
            END IF
      END IF
C
C     Compute -R*P'*x = -R*z.
C
      DO 120 J = 1, N
         L     = IPVT(J)
         RX(J) = -X(L)
  120 CONTINUE
C
      DO 130 I = 1, NTHS, BSN
         CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', BSN, R(I,1),
     $               LDR, RX(I), 1 )
  130 CONTINUE
C
      IF ( ST.GT.0 ) THEN
         CALL DGEMV( 'NoTranspose', NTHS, ST, ONE, R(1,BSN+1), LDR,
     $               RX(NTHS+1), 1, ONE, RX, 1 )
         CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', ST,
     $               R(NTHS+1,BSN+1), LDR, RX(NTHS+1), 1 )
      END IF
C
C     Termination. If PAR = 0, set S.
C
      IF ( ITER.EQ.0 ) THEN
         PAR = ZERO
         I   = 1
C
         DO 150 K = 1, BN
C
            DO 140 J = 1, BSN
               DWORK(I) = R(I,J)
               CALL DCOPY( BSN-J+1, R(I,J), LDR, R(I,J), 1 )
               I = I + 1
  140       CONTINUE
C
  150    CONTINUE
C
         IF ( ST.GT.0 ) THEN
C
            DO 160 J = BSN + 1, BSN + ST
               CALL DCOPY( NTHS, R(1,J), 1, DWORK(N+J-BSN), ST )
               DWORK(I) = R(I,J)
               CALL DCOPY( BSN+ST-J+1, R(I,J), LDR, R(I,J), 1 )
               I = I + 1
  160       CONTINUE
C
         END IF
      ELSE
C
         DO 170 K = N + 1, N + ST*NTHS
            DWORK(K) = DWORK(K+N)
  170    CONTINUE
C
      END IF
C
      RETURN
C
C *** Last line of NF01BP ***
      END
