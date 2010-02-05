      SUBROUTINE SB02MU( DICO, HINV, UPLO, N, A, LDA, G, LDG, Q, LDQ, S,
     $                   LDS, IWORK, DWORK, LDWORK, INFO )
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
C     To construct the 2n-by-2n Hamiltonian or symplectic matrix S
C     associated to the linear-quadratic optimization problem, used to
C     solve the continuous- or discrete-time algebraic Riccati equation,
C     respectively.
C
C     For a continuous-time problem, S is defined by
C
C             (  A  -G )
C         S = (        ),                                       (1)
C             ( -Q  -A')
C
C     and for a discrete-time problem by
C
C                 -1       -1
C             (  A        A  *G     )
C         S = (   -1           -1   ),                          (2)
C             ( QA     A' + Q*A  *G )
C
C     or
C
C                       -T         -T
C             (  A + G*A  *Q   -G*A   )
C         S = (      -T            -T ),                        (3)
C             (    -A  *Q         A   )
C
C     where A, G, and Q are N-by-N matrices, with G and Q symmetric.
C     Matrix A must be nonsingular in the discrete-time case.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system as follows:
C             = 'C':  Continuous-time system;
C             = 'D':  Discrete-time system.
C
C     HINV    CHARACTER*1
C             If DICO = 'D', specifies which of the matrices (2) or (3)
C             is constructed, as follows:
C             = 'D':  The matrix S in (2) is constructed;
C             = 'I':  The (inverse) matrix S in (3) is constructed.
C             HINV is not referenced if DICO = 'C'.
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the matrices G and Q is
C             stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C             On exit, if DICO = 'D', and INFO = 0, the leading N-by-N
C                                                     -1
C             part of this array contains the matrix A  .
C             Otherwise, the array A is unchanged on exit.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     G       (input) DOUBLE PRECISION array, dimension (LDG,N)
C             The leading N-by-N upper triangular part (if UPLO = 'U')
C             or lower triangular part (if UPLO = 'L') of this array
C             must contain the upper triangular part or lower triangular
C             part, respectively, of the symmetric matrix G. The stricly
C             lower triangular part (if UPLO = 'U') or stricly upper
C             triangular part (if UPLO = 'L') is not referenced.
C
C     LDG     INTEGER
C             The leading dimension of array G.  LDG >= MAX(1,N).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N)
C             The leading N-by-N upper triangular part (if UPLO = 'U')
C             or lower triangular part (if UPLO = 'L') of this array
C             must contain the upper triangular part or lower triangular
C             part, respectively, of the symmetric matrix Q. The stricly
C             lower triangular part (if UPLO = 'U') or stricly upper
C             triangular part (if UPLO = 'L') is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,N).
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N)
C             If INFO = 0, the leading 2N-by-2N part of this array
C             contains the Hamiltonian or symplectic matrix of the
C             problem.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,2*N).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK; if DICO = 'D', DWORK(2) returns the reciprocal
C             condition number of the given matrix  A.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 1          if DICO = 'C';
C             LDWORK >= MAX(2,4*N) if DICO = 'D'.
C             For optimum performance LDWORK should be larger, if
C             DICO = 'D'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = i:  if the leading i-by-i (1 <= i <= N) upper triangular
C                   submatrix of A is singular in discrete-time case;
C             = N+1:  if matrix A is numerically singular in discrete-
C                   time case.
C
C     METHOD
C
C     For a continuous-time problem, the 2n-by-2n Hamiltonian matrix (1)
C     is constructed.
C     For a discrete-time problem, the 2n-by-2n symplectic matrix (2) or
C     (3) - the inverse of the matrix in (2) - is constructed.
C
C     NUMERICAL ASPECTS
C
C     The discrete-time case needs the inverse of the matrix A, hence
C     the routine should not be used when A is ill-conditioned.
C                               3
C     The algorithm requires 0(n ) floating point operations in the
C     discrete-time case.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, discrete-time system, optimal regulator, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDWORK, N
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*),
     $                  S(LDS,*)
C     .. Local Scalars ..
      LOGICAL           DISCR, LHINV, LUPLO
      INTEGER           I, J, MAXWRK, N2, NJ, NP1
      DOUBLE PRECISION  ANORM, RCOND
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGECON, DGEMM, DGETRF, DGETRI, DGETRS,
     $                  DLACPY, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO  = 0
      N2 = N + N
      DISCR = LSAME( DICO,  'D' )
      LUPLO = LSAME( UPLO,  'U' )
      IF( DISCR ) THEN
         LHINV = LSAME( HINV, 'D' )
      ELSE
         LHINV = .FALSE.
      END IF
C
C     Test the input scalar arguments.
C
      IF( .NOT.DISCR .AND. .NOT.LSAME( DICO, 'C' ) ) THEN
         INFO = -1
      ELSE IF( DISCR ) THEN
         IF( .NOT.LHINV .AND. .NOT.LSAME( HINV, 'I' ) )
     $      INFO = -2
      END IF
      IF( .NOT.LUPLO .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDS.LT.MAX( 1, N2 ) ) THEN
         INFO = -12
      ELSE IF( ( LDWORK.LT.1 ) .OR.
     $         ( DISCR .AND. LDWORK.LT.MAX( 2, 4*N ) ) ) THEN
         INFO = -15
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB02MU', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         IF ( DISCR ) DWORK(2) = ONE
         RETURN
      END IF
C
C     The code tries to exploit data locality as much as possible.
C
      IF ( .NOT.LHINV ) THEN
         CALL DLACPY( 'Full', N, N, A, LDA, S, LDS )
C
C        Construct Hamiltonian matrix in the continuous-time case, or
C        prepare symplectic matrix in (3) in the discrete-time case:
C
C        Construct full Q in S(N+1:2*N,1:N) and change the sign, and
C        construct full G in S(1:N,N+1:2*N) and change the sign.
C
         DO 200 J = 1, N
            NJ = N + J
            IF ( LUPLO ) THEN
C
               DO 20 I = 1, J
                  S(N+I,J) = -Q(I,J)
   20          CONTINUE
C
               DO 40 I = J + 1, N
                  S(N+I,J) = -Q(J,I)
   40          CONTINUE
C
               DO 60 I = 1, J
                  S(I,NJ) = -G(I,J)
   60          CONTINUE
C
               DO 80 I = J + 1, N
                  S(I,NJ) = -G(J,I)
   80          CONTINUE
C
            ELSE
C
               DO 100 I = 1, J - 1
                  S(N+I,J) = -Q(J,I)
  100          CONTINUE
C
               DO 120 I = J, N
                  S(N+I,J) = -Q(I,J)
  120          CONTINUE
C
               DO 140 I = 1, J - 1
                  S(I,NJ) = -G(J,I)
  140          CONTINUE
C
               DO 180 I = J, N
                  S(I,NJ) = -G(I,J)
  180          CONTINUE
C
            END IF
  200    CONTINUE
C
         IF ( .NOT.DISCR ) THEN
C
            DO 240 J = 1, N
               NJ = N + J
C
               DO 220 I = 1, N
                  S(N+I,NJ) = -A(J,I)
  220          CONTINUE
C
  240       CONTINUE
C
            DWORK(1) = ONE
         END IF
      END IF
C
      IF ( DISCR ) THEN
C
C        Construct the symplectic matrix (2) or (3) in the discrete-time
C        case.
C
C        Compute workspace.
C        (Note: Comments in the code beginning "Workspace:" describe the
C        minimal amount of workspace needed at that point in the code,
C        as well as the preferred amount for good performance.
C        NB refers to the optimal block size for the immediately
C        following subroutine, as returned by ILAENV.)
C
         MAXWRK = MAX( 4*N,
     $                 N*ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 ) )
         NP1 = N + 1
C
         IF ( LHINV ) THEN
C
C           Put  A'  in  S(N+1:2*N,N+1:2*N).
C
            DO 260 I = 1, N
               CALL DCOPY( N, A(I, 1), LDA, S(NP1,N+I), 1 )
  260       CONTINUE
C
         END IF
C
C        Compute the norm of the matrix A.
C
         ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK )
C
C        Compute the LU factorization of A.
C
         CALL DGETRF( N, N, A, LDA, IWORK, INFO )
C
C        Return if INFO is non-zero.
C
         IF( INFO.GT.0 ) THEN
            DWORK(2) = ZERO
            RETURN
         END IF
C
C        Compute the reciprocal of the condition number of A.
C        Workspace: need 4*N.
C
         CALL DGECON( '1-norm', N, A, LDA, ANORM, RCOND, DWORK,
     $                IWORK(NP1), INFO )
C
C        Return if the matrix is singular to working precision.
C
         IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
            INFO = N + 1
            DWORK(2) = RCOND
            RETURN
         END IF
C
         IF ( LHINV ) THEN
C
C           Compute S in (2).
C
C           Construct full Q in S(N+1:2*N,1:N).
C
            IF ( LUPLO ) THEN
               DO 270 J = 1, N - 1
                  CALL DCOPY( J, Q(1,J), 1, S(NP1,J), 1 )
                  CALL DCOPY( N-J, Q(J,J+1), LDQ, S(NP1+J,J), 1 )
  270          CONTINUE
               CALL DCOPY( N, Q(1,N), 1, S(NP1,N), 1 )
            ELSE
               CALL DCOPY( N, Q(1,1), 1, S(NP1,1), 1 )
               DO 280 J = 2, N
                  CALL DCOPY( J-1, Q(J,1), LDQ, S(NP1,J), 1 )
                  CALL DCOPY( N-J+1, Q(J,J), 1, S(N+J,J), 1 )
  280          CONTINUE
            END IF
C
C           Compute the solution matrix  X  of the system  X*A = Q  by
C                                                                    -1
C           solving  A'*X' = Q and transposing the result to get  Q*A  .
C
            CALL DGETRS( 'Transpose', N, N, A, LDA, IWORK, S(NP1,1),
     $                   LDS, INFO )
C
            DO 300 J = 1, N - 1
               CALL DSWAP( N-J, S(NP1+J,J), 1, S(N+J,J+1), LDS )
  300       CONTINUE
C
C           Construct full G in S(1:N,N+1:2*N).
C
            IF ( LUPLO ) THEN
               DO 310 J = 1, N - 1
                  CALL DCOPY( J, G(1,J), 1, S(1,N+J), 1 )
                  CALL DCOPY( N-J, G(J,J+1), LDG, S(J+1,N+J), 1 )
  310          CONTINUE
               CALL DCOPY( N, G(1,N), 1, S(1,N2), 1 )
            ELSE
               CALL DCOPY( N, G(1,1), 1, S(1,NP1), 1 )
               DO 320 J = 2, N
                  CALL DCOPY( J-1, G(J,1), LDG, S(1,N+J), 1 )
                  CALL DCOPY( N-J+1, G(J,J), 1, S(J,N+J), 1 )
  320          CONTINUE
            END IF
C                            -1
C           Compute  A' + Q*A  *G  in  S(N+1:2N,N+1:2N).
C
            CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                  S(NP1,1), LDS, S(1,NP1), LDS, ONE, S(NP1,NP1),
     $                  LDS )
C
C           Compute the solution matrix  Y  of the system  A*Y = G.
C
            CALL DGETRS( 'No transpose', N, N, A, LDA, IWORK, S(1,NP1),
     $                   LDS, INFO )
C
C           Compute the inverse of  A  in situ.
C           Workspace: need N;  prefer N*NB.
C
            CALL DGETRI( N, A, LDA, IWORK, DWORK, LDWORK, INFO )
C                  -1
C           Copy  A    in  S(1:N,1:N).
C
            CALL DLACPY( 'Full', N, N, A, LDA, S, LDS )
C
         ELSE
C
C           Compute S in (3) using the already prepared part.
C
C           Compute the solution matrix  X'  of the system  A*X' = -G
C                                                       -T
C           and transpose the result to obtain  X = -G*A  .
C
            CALL DGETRS( 'No transpose', N, N, A, LDA, IWORK, S(1,NP1),
     $                   LDS, INFO )
C
            DO 340 J = 1, N - 1
               CALL DSWAP( N-J, S(J+1,N+J), 1, S(J,NP1+J), LDS )
  340       CONTINUE
C                           -T
C           Compute  A + G*A  *Q  in  S(1:N,1:N).
C
            CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                  S(1,NP1), LDS, S(NP1, 1), LDS, ONE, S, LDS )
C
C           Compute the solution matrix  Y  of the system  A'*Y = -Q.
C
            CALL DGETRS( 'Transpose', N, N, A, LDA, IWORK, S(NP1,1),
     $                   LDS, INFO )
C
C           Compute the inverse of  A  in situ.
C           Workspace: need N;  prefer N*NB.
C
            CALL DGETRI( N, A, LDA, IWORK, DWORK, LDWORK, INFO )
C                  -T
C           Copy  A    in  S(N+1:2N,N+1:2N).
C
            DO 360 J = 1, N
               CALL DCOPY( N, A(J,1), LDA, S(NP1,N+J), 1 )
  360       CONTINUE
C
         END IF
         DWORK(1) = MAXWRK
         DWORK(2) = RCOND
      END IF
C
C *** Last line of SB02MU ***
      RETURN
      END
