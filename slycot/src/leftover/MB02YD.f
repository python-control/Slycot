      SUBROUTINE MB02YD( COND, N, R, LDR, IPVT, DIAG, QTB, RANK, X, TOL,
     $                   DWORK, LDWORK, INFO )
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
C           A*x = b ,     D*x = 0 ,
C
C     in the least squares sense, where A is an m-by-n matrix,
C     D is an n-by-n diagonal matrix, and b is an m-vector.
C     It is assumed that a QR factorization, with column pivoting, of A
C     is available, that is, A*P = Q*R, where P is a permutation matrix,
C     Q has orthogonal columns, and R is an upper triangular matrix
C     with diagonal elements of nonincreasing magnitude.
C     The routine needs the full upper triangle of R, the permutation
C     matrix P, and the first n components of Q'*b (' denotes the
C     transpose). The system A*x = b, D*x = 0, is then equivalent to
C
C           R*z = Q'*b ,  P'*D*P*z = 0 ,                             (1)
C
C     where x = P*z. If this system does not have full rank, then a
C     least squares solution is obtained. On output, MB02YD also
C     provides an upper triangular matrix S such that
C
C           P'*(A'*A + D*D)*P = S'*S .
C
C     The system (1) is equivalent to S*z = c , where c contains the
C     first n components of the vector obtained by applying to
C     [ (Q'*b)'  0 ]' the transformations which triangularized
C     [ R'  P'*D*P ]', getting S.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COND    CHARACTER*1
C             Specifies whether the condition of the matrix S should be
C             estimated, as follows:
C             = 'E' :  use incremental condition estimation and store
C                      the numerical rank of S in RANK;
C             = 'N' :  do not use condition estimation, but check the
C                      diagonal entries of S for zero values;
C             = 'U' :  use the rank already stored in RANK.
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
C             matrix D.
C
C     QTB     (input) DOUBLE PRECISION array, dimension (N)
C             This array must contain the first n elements of the
C             vector Q'*b.
C
C     RANK    (input or output) INTEGER
C             On entry, if COND = 'U', this parameter must contain the
C             (numerical) rank of the matrix S.
C             On exit, if COND = 'E' or 'N', this parameter contains
C             the numerical rank of the matrix S, estimated according
C             to the value of COND.
C
C     X       (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the least squares solution of the
C             system A*x = b, D*x = 0.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If COND = 'E', the tolerance to be used for finding the
C             rank of the matrix S. If the user sets TOL > 0, then the
C             given value of TOL is used as a lower bound for the
C             reciprocal condition number;  a (sub)matrix whose
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
C     Standard plane rotations are used to annihilate the elements of
C     the diagonal matrix D, updating the upper triangular matrix R
C     and the first n elements of the vector Q'*b. A basic least squares
C     solution is computed.
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
C     This routine is a LAPACK-based modification of QRSOLV from the
C     MINPACK package [1], and with optional condition estimation.
C     The option COND = 'U' is useful when dealing with several
C     right-hand side vectors.
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
      DOUBLE PRECISION  ZERO, SVLMAX
      PARAMETER         ( ZERO = 0.0D0, SVLMAX = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COND
      INTEGER           INFO, LDR, LDWORK, N, RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IPVT(*)
      DOUBLE PRECISION  DIAG(*), DWORK(*), QTB(*), R(LDR,*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, QTBPJ, SN, TEMP, TOLDEF
      INTEGER           I, J, K, L
      LOGICAL           ECOND, NCOND, UCOND
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(3)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLARTG, DROT, DSWAP, MB03OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
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
      ELSEIF ( UCOND .AND. ( RANK.LT.0 .OR. RANK.GT.N ) ) THEN
         INFO = -8
      ELSEIF ( LDWORK.LT.2*N .OR. ( ECOND .AND. LDWORK.LT.4*N ) ) THEN
         INFO = -12
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02YD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         IF ( .NOT.UCOND )
     $      RANK = 0
         RETURN
      END IF
C
C     Copy R and Q'*b to preserve input and initialize S.
C     In particular, save the diagonal elements of R in X.
C
      DO 20 J = 1, N
         X(J) = R(J,J)
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10    CONTINUE
   20 CONTINUE
C
      CALL DCOPY( N, QTB, 1, DWORK(N+1), 1 )
C
C     Eliminate the diagonal matrix D using Givens rotations.
C
      DO 50 J = 1, N
C
C        Prepare the row of D to be eliminated, locating the
C        diagonal element using P from the QR factorization.
C
         L = IPVT(J)
         IF ( DIAG(L).NE.ZERO ) THEN
            QTBPJ = ZERO
            DWORK(J) = DIAG(L)
C
            DO 30 K = J + 1, N
               DWORK(K) = ZERO
   30       CONTINUE
C
C           The transformations to eliminate the row of D modify only
C           a single element of Q'*b beyond the first n, which is
C           initially zero.
C
            DO 40 K = J, N
C
C              Determine a Givens rotation which eliminates the
C              appropriate element in the current row of D.
C
               IF ( DWORK(K).NE.ZERO ) THEN
C
                  CALL DLARTG( R(K,K), DWORK(K), CS, SN, TEMP )
C
C                 Compute the modified diagonal element of R and
C                 the modified elements of (Q'*b,0).
C                 Accumulate the tranformation in the row of S.
C
                  TEMP  =  CS*DWORK(N+K) + SN*QTBPJ
                  QTBPJ = -SN*DWORK(N+K) + CS*QTBPJ
                  DWORK(N+K) = TEMP
                  CALL DROT( N-K+1, R(K,K), 1, DWORK(K), 1, CS, SN )
C
               END IF
   40       CONTINUE
C
         END IF
C
C        Store the diagonal element of S and, if COND <> 'E', restore
C        the corresponding diagonal element of R.
C
         DWORK(J) = R(J,J)
         IF ( .NOT.ECOND )
     $      R(J,J) = X(J)
   50 CONTINUE
C
C     Solve the triangular system for z. If the system is singular,
C     then obtain a least squares solution.
C
      IF ( ECOND ) THEN
         TOLDEF = TOL
         IF ( TOLDEF.LE.ZERO ) THEN
C
C           Use the default tolerance in rank determination.
C
            TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
         END IF
C
C        Interchange the strict upper and lower triangular parts of R.
C
         DO 60 J = 2, N
            CALL DSWAP( J-1, R(1,J), 1, R(J,1), LDR )
   60    CONTINUE
C
C        Estimate the reciprocal condition number of S and set the rank.
C        Additional workspace: 2*N.
C
         CALL MB03OD( 'No QR', N, N, R, LDR, IPVT, TOLDEF, SVLMAX,
     $                DWORK, RANK, DUM, DWORK(2*N+1), LDWORK-2*N,
     $                INFO )
         R(1,1) = X(1)
C
C        Restore the strict upper and lower triangular parts of R.
C
         DO 70 J = 2, N
            CALL DSWAP( J-1, R(1,J), 1, R(J,1), LDR )
            R(J,J) = X(J)
   70    CONTINUE
C
      ELSEIF ( NCOND ) THEN
C
C        Determine rank(S) by checking zero diagonal entries.
C
         RANK = N
C
         DO 80 J = 1, N
            IF ( DWORK(J).EQ.ZERO .AND. RANK.EQ.N )
     $         RANK = J - 1
   80    CONTINUE
C
      END IF
C
      DUM(1) = ZERO
      IF ( RANK.LT.N )
     $   CALL DCOPY( N-RANK, DUM, 0, DWORK(N+RANK+1), 1 )
C
C     Solve S*z = c using back substitution.
C
      DO 100 J = RANK, 1, -1
         TEMP = ZERO
C
         DO 90 I = J + 1, RANK
            TEMP = TEMP + R(I,J)*DWORK(N+I)
   90    CONTINUE
C
         DWORK(N+J) = ( DWORK(N+J) - TEMP )/DWORK(J)
  100 CONTINUE
C
C     Permute the components of z back to components of x.
C
      DO 110 J = 1, N
         L    = IPVT(J)
         X(L) = DWORK(N+J)
  110 CONTINUE
C
      RETURN
C
C *** Last line of MB02YD ***
      END
