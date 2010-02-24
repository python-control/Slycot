      SUBROUTINE MD03BY( COND, N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR,
     $                   RANK, X, RX, TOL, DWORK, LDWORK, INFO )
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
C     To determine a value for the parameter PAR such that if x solves
C     the system
C
C           A*x = b ,     sqrt(PAR)*D*x = 0 ,
C
C     in the least squares sense, where A is an m-by-n matrix, D is an
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
C     It is assumed that a QR factorization, with column pivoting, of A
C     is available, that is, A*P = Q*R, where P is a permutation matrix,
C     Q has orthogonal columns, and R is an upper triangular matrix
C     with diagonal elements of nonincreasing magnitude.
C     The routine needs the full upper triangle of R, the permutation
C     matrix P, and the first n components of Q'*b (' denotes the
C     transpose). On output, MD03BY also provides an upper triangular
C     matrix S such that
C
C           P'*(A'*A + PAR*D*D)*P = S'*S .
C
C     Matrix S is used in the solution process.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COND    CHARACTER*1
C             Specifies whether the condition of the matrices R and S
C             should be estimated, as follows:
C             = 'E' :  use incremental condition estimation for R and S;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of R and S for zero values;
C             = 'U' :  use the rank already stored in RANK (for R).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix R.  N >= 0.
C
C     R       (input/output) DOUBLE PRECISION array, dimension (LDR, N)
C             On entry, the leading N-by-N upper triangular part of this
C             array must contain the upper triangular matrix R.
C             On exit, the full upper triangle is unaltered, and the
C             strict lower triangle contains the strict upper triangle
C             (transposed) of the upper triangular matrix S.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     IPVT    (input) INTEGER array, dimension (N)
C             This array must define the permutation matrix P such that
C             A*P = Q*R. Column j of P is column IPVT(j) of the identity
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
C     RANK    (input or output) INTEGER
C             On entry, if COND = 'U', this parameter must contain the
C             (numerical) rank of the matrix R.
C             On exit, this parameter contains the numerical rank of
C             the matrix S.
C
C     X       (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the least squares solution of the
C             system A*x = b, sqrt(PAR)*D*x = 0.
C
C     RX      (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the matrix-vector product -R*P'*x.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             rank of the matrices R and S. If the user sets TOL > 0,
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
C             diagonal elements of the upper triangular matrix S.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 4*N, if COND =  'E';
C             LDWORK >= 2*N, if COND <> 'E'.
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
C     The algorithm computes the Gauss-Newton direction. A least squares
C     solution is found if the Jacobian is rank deficient. If the Gauss-
C     Newton direction is not acceptable, then an iterative algorithm
C     obtains improved lower and upper bounds for the parameter PAR.
C     Only a few iterations are generally needed for convergence of the
C     algorithm. If, however, the limit of ITMAX = 10 iterations is
C     reached, then the output PAR will contain the best value obtained
C     so far. If the Gauss-Newton step is acceptable, it is stored in x,
C     and PAR is set to zero, hence S = R.
C
C     REFERENCES
C
C     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E.
C         User's Guide for MINPACK-1.
C         Applied Math. Division, Argonne National Laboratory, Argonne,
C         Illinois, Report ANL-80-74, 1980.
C
C     NUMERICAL ASPECTS
C                               2
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     FURTHER COMMENTS
C
C     This routine is a LAPACK-based modification of LMPAR from the
C     MINPACK package [1], and with optional condition estimation.
C     The option COND = 'U' is useful when dealing with several
C     right-hand side vectors, but RANK should be reset.
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
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005.
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
      DOUBLE PRECISION  P1, P001, ZERO, SVLMAX
      PARAMETER         ( P1  = 1.0D-1, P001 = 1.0D-3, ZERO = 0.0D0,
     $                    SVLMAX = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COND
      INTEGER           INFO, LDR, LDWORK, N, RANK
      DOUBLE PRECISION  DELTA, PAR, TOL
C     .. Array Arguments ..
      INTEGER           IPVT(*)
      DOUBLE PRECISION  DIAG(*), DWORK(*), QTB(*), R(LDR,*), RX(*), X(*)
C     .. Local Scalars ..
      INTEGER           ITER, J, L, N2
      DOUBLE PRECISION  DMINO, DWARF, DXNORM, FP, GNORM, PARC, PARL,
     $                  PARU, TEMP, TOLDEF
      LOGICAL           ECOND, NCOND, SING, UCOND
      CHARACTER         CONDL
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DLAMCH, DNRM2
      LOGICAL           LSAME
      EXTERNAL          DDOT, DLAMCH, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSWAP, DTRMV, DTRSV, MB02YD,
     $                  MB03OD, XERBLA
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
      IF( .NOT.( ECOND .OR. NCOND .OR. UCOND ) ) THEN
         INFO = -1
      ELSEIF( N.LT.0 ) THEN
         INFO = -2
      ELSEIF ( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSEIF ( DELTA.LE.ZERO ) THEN
         INFO = -8
      ELSEIF( PAR.LT.ZERO ) THEN
         INFO = -9
      ELSEIF ( UCOND .AND. ( RANK.LT.0 .OR. RANK.GT.N ) ) THEN
         INFO = -10
      ELSEIF ( LDWORK.LT.2*N .OR. ( ECOND .AND. LDWORK.LT.4*N ) ) THEN
         INFO = -15
      ELSEIF ( N.GT.0 ) THEN
         DMINO = DIAG(1)
         SING  = .FALSE.
C
         DO 10 J = 1, N
            IF ( DIAG(J).LT.DMINO )
     $         DMINO = DIAG(J)
            SING = SING .OR. DIAG(J).EQ.ZERO
   10    CONTINUE
C
         IF ( SING )
     $      INFO = -6
      END IF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MD03BY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         PAR  = ZERO
         RANK = 0
         RETURN
      END IF
C
C     DWARF is the smallest positive magnitude.
C
      DWARF = DLAMCH( 'Underflow' )
      N2    = N
C
C     Estimate the rank of R, if required.
C
      IF ( ECOND ) THEN
         N2   = 2*N
         TEMP = TOL
         IF ( TEMP.LE.ZERO ) THEN
C
C           Use the default tolerance in rank determination.
C
            TEMP = DBLE( N )*DLAMCH( 'Epsilon' )
         END IF
C
C        Estimate the reciprocal condition number of R and set the rank.
C        Workspace: 2*N.
C
         CALL MB03OD( 'No QR', N, N, R, LDR, IPVT, TEMP, SVLMAX, DWORK,
     $                RANK, DUM, DWORK, LDWORK, INFO )
C
      ELSEIF ( NCOND ) THEN
         J = 1
C
   20    CONTINUE
            IF ( R(J,J).NE.ZERO ) THEN
               J = J + 1
               IF ( J.LE.N )
     $            GO TO 20
            END IF
C
         RANK = J - 1
      END IF
C
C     Compute and store in x the Gauss-Newton direction. If the
C     Jacobian is rank-deficient, obtain a least squares solution.
C     The array RX is used as workspace.
C
      CALL DCOPY( RANK, QTB, 1, RX, 1 )
      DUM(1) = ZERO
      IF ( RANK.LT.N )
     $   CALL DCOPY( N-RANK, DUM, 0, RX(RANK+1), 1 )
      CALL DTRSV( 'Upper', 'No transpose', 'Non unit', RANK, R, LDR,
     $            RX, 1 )
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
         IF ( UCOND ) THEN
            IF ( LDWORK.GE.4*N ) THEN
               CONDL  = 'E'
               TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
            ELSE
               CONDL  = 'N'
               TOLDEF = TOL
            END IF
         ELSE
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
            DO 50 J = 1, N
               L = IPVT(J)
               RX(J) = DIAG(L)*( DWORK(L)/DXNORM )
   50       CONTINUE
C
            CALL DTRSV( 'Upper', 'Transpose', 'Non unit', N, R, LDR,
     $                  RX, 1 )
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
C        Calculate an upper bound, PARU, for the zero of the function.
C
         DO 60 J = 1, N
            L = IPVT(J)
            RX(J) = DDOT( J, R(1,J), 1, QTB, 1 )/DIAG(L)
   60    CONTINUE
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
   70    CONTINUE
            ITER = ITER + 1
C
C           Evaluate the function at the current value of PAR.
C
            IF ( PAR.EQ.ZERO )
     $         PAR = MAX( DWARF, P001*PARU )
            TEMP = SQRT( PAR )
C
            DO 80 J = 1, N
               RX(J) = TEMP*DIAG(J)
   80       CONTINUE
C
C           Solve the system A*x = b , sqrt(PAR)*D*x = 0 , in a least
C           square sense. The first N elements of DWORK contain the
C           diagonal elements of the upper triangular matrix S, and
C           the next N elements contain the vector z, so that x = P*z.
C           The vector z is preserved if COND = 'E'.
C           Workspace:   4*N, if CONDL =  'E';
C                        2*N, if CONDL <> 'E'.
C
            CALL MB02YD( CONDL, N, R, LDR, IPVT, RX, QTB, RANK, X,
     $                   TOLDEF, DWORK, LDWORK, INFO )
C
            DO 90 J = 1, N
               DWORK(N2+J) = DIAG(J)*X(J)
   90       CONTINUE
C
            DXNORM = DNRM2( N, DWORK(N2+1), 1 )
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
               DO 100 J = 1, RANK
                  L = IPVT(J)
                  RX(J) = DIAG(L)*( DWORK(N2+L)/DXNORM )
  100          CONTINUE
C
               IF ( RANK.LT.N )
     $            CALL DCOPY( N-RANK, DUM, 0, RX(RANK+1), 1 )
               CALL DSWAP( N, R, LDR+1, DWORK, 1 )
               CALL DTRSV( 'Lower', 'No transpose', 'Non Unit', RANK,
     $                     R, LDR, RX, 1 )
               CALL DSWAP( N, R, LDR+1, DWORK, 1 )
               TEMP = DNRM2( RANK, RX, 1 )
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
               GO TO 70
            END IF
      END IF
C
C     Compute -R*P'*x = -R*z.
C
      IF ( ECOND .AND. ITER.GT.0 ) THEN
C
         DO 110 J = 1, N
            RX(J) = -DWORK(N+J)
  110    CONTINUE
C
         CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', N, R, LDR,
     $               RX, 1 )
      ELSE
C
         DO 120 J = 1, N
            RX(J) = ZERO
            L     = IPVT(J)
            CALL DAXPY( J, -X(L), R(1,J), 1, RX, 1 )
  120    CONTINUE
C
      END IF
C
C     Termination. If PAR = 0, set S.
C
      IF ( ITER.EQ.0 ) THEN
         PAR = ZERO
C
         DO 130 J = 1, N - 1
            DWORK(J) = R(J,J)
            CALL DCOPY( N-J, R(J,J+1), LDR, R(J+1,J), 1 )
  130    CONTINUE
C
         DWORK(N) = R(N,N)
      END IF
C
      RETURN
C
C *** Last line of MD03BY ***
      END
