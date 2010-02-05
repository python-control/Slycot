      SUBROUTINE SB02MD( DICO, HINV, UPLO, SCAL, SORT, N, A, LDA, G,
     $                   LDG, Q, LDQ, RCOND, WR, WI, S, LDS, U, LDU,
     $                   IWORK, DWORK, LDWORK, BWORK, INFO )
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
C     To solve for X either the continuous-time algebraic Riccati
C     equation
C                              -1
C        Q + A'*X + X*A - X*B*R  B'*X = 0                            (1)
C
C     or the discrete-time algebraic Riccati equation
C                                        -1
C        X = A'*X*A - A'*X*B*(R + B'*X*B)  B'*X*A + Q                (2)
C
C     where A, B, Q and R are N-by-N, N-by-M, N-by-N and M-by-M matrices
C     respectively, with Q symmetric and R symmetric nonsingular; X is
C     an N-by-N symmetric matrix.
C                       -1
C     The matrix G = B*R  B' must be provided on input, instead of B and
C     R, that is, for instance, the continuous-time equation
C
C        Q + A'*X + X*A - X*G*X = 0                                  (3)
C
C     is solved, where G is an N-by-N symmetric matrix. SLICOT Library
C     routine SB02MT should be used to compute G, given B and R. SB02MT
C     also enables to solve Riccati equations corresponding to optimal
C     problems with coupling terms.
C
C     The routine also returns the computed values of the closed-loop
C     spectrum of the optimal system, i.e., the stable eigenvalues
C     lambda(1),...,lambda(N) of the corresponding Hamiltonian or
C     symplectic matrix associated to the optimal problem.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of Riccati equation to be solved as
C             follows:
C             = 'C':  Equation (3), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     HINV    CHARACTER*1
C             If DICO = 'D', specifies which symplectic matrix is to be
C             constructed, as follows:
C             = 'D':  The matrix H in (5) (see METHOD) is constructed;
C             = 'I':  The inverse of the matrix H in (5) is constructed.
C             HINV is not used if DICO = 'C'.
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the matrices G and Q is
C             stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
C
C     SCAL    CHARACTER*1
C             Specifies whether or not a scaling strategy should be
C             used, as follows:
C             = 'G':  General scaling should be used;
C             = 'N':  No scaling should be used.
C
C     SORT    CHARACTER*1
C             Specifies which eigenvalues should be obtained in the top
C             of the Schur form, as follows:
C             = 'S':  Stable   eigenvalues come first;
C             = 'U':  Unstable eigenvalues come first.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, Q, G and X.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the coefficient matrix A of the equation.
C             On exit, if DICO = 'D', and INFO = 0 or INFO > 1, the
C                                                                    -1
C             leading N-by-N part of this array contains the matrix A  .
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
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, the leading N-by-N upper triangular part (if
C             UPLO = 'U') or lower triangular part (if UPLO = 'L') of
C             this array must contain the upper triangular part or lower
C             triangular part, respectively, of the symmetric matrix Q.
C             The stricly lower triangular part (if UPLO = 'U') or
C             stricly upper triangular part (if UPLO = 'L') is not used.
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains the solution matrix X of the problem.
C
C     LDQ     INTEGER
C             The leading dimension of array N.  LDQ >= MAX(1,N).
C
C     RCOND   (output) DOUBLE PRECISION
C             An estimate of the reciprocal of the condition number (in
C             the 1-norm) of the N-th order system of algebraic
C             equations from which the solution matrix X is obtained.
C
C     WR      (output) DOUBLE PRECISION array, dimension (2*N)
C     WI      (output) DOUBLE PRECISION array, dimension (2*N)
C             If INFO = 0 or INFO = 5, these arrays contain the real and
C             imaginary parts, respectively, of the eigenvalues of the
C             2N-by-2N matrix S, ordered as specified by SORT (except
C             for the case HINV = 'D', when the order is opposite to
C             that specified by SORT). The leading N elements of these
C             arrays contain the closed-loop spectrum of the system
C                           -1
C             matrix A - B*R  *B'*X, if DICO = 'C', or of the matrix
C                               -1
C             A - B*(R + B'*X*B)  B'*X*A, if DICO = 'D'. Specifically,
C                lambda(k) = WR(k) + j*WI(k), for k = 1,2,...,N.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N)
C             If INFO = 0 or INFO = 5, the leading 2N-by-2N part of this
C             array contains the ordered real Schur form S of the
C             Hamiltonian or symplectic matrix H. That is,
C
C                    (S   S  )
C                    ( 11  12)
C                S = (       ),
C                    (0   S  )
C                    (     22)
C
C             where S  , S   and S   are N-by-N matrices.
C                    11   12      22
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,2*N).
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N)
C             If INFO = 0 or INFO = 5, the leading 2N-by-2N part of this
C             array contains the transformation matrix U which reduces
C             the Hamiltonian or symplectic matrix H to the ordered real
C             Schur form S. That is,
C
C                    (U   U  )
C                    ( 11  12)
C                U = (       ),
C                    (U   U  )
C                    ( 21  22)
C
C             where U  , U  , U   and U   are N-by-N matrices.
C                    11   12   21      22
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,2*N).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) returns the scaling factor used
C             (set to 1 if SCAL = 'N'), also set if INFO = 5;
C             if DICO = 'D', DWORK(3) returns the reciprocal condition
C             number of the given matrix  A.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(2,6*N) if DICO = 'C';
C             LDWORK >= MAX(3,6*N) if DICO = 'D'.
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if matrix A is (numerically) singular in discrete-
C                   time case;
C             = 2:  if the Hamiltonian or symplectic matrix H cannot be
C                   reduced to real Schur form;
C             = 3:  if the real Schur form of the Hamiltonian or
C                   symplectic matrix H cannot be appropriately ordered;
C             = 4:  if the Hamiltonian or symplectic matrix H has less
C                   than N stable eigenvalues;
C             = 5:  if the N-th order system of linear algebraic
C                   equations, from which the solution matrix X would
C                   be obtained, is singular to working precision.
C
C     METHOD
C
C     The method used is the Schur vector approach proposed by Laub.
C     It is assumed that [A,B] is a stabilizable pair (where for (3) B
C     is any matrix such that B*B' = G with rank(B) = rank(G)), and
C     [E,A] is a detectable pair, where E is any matrix such that
C     E*E' = Q with rank(E) = rank(Q). Under these assumptions, any of
C     the algebraic Riccati equations (1)-(3) is known to have a unique
C     non-negative definite solution. See [2].
C     Now consider the 2N-by-2N Hamiltonian or symplectic matrix
C
C                 ( A   -G )
C            H =  (        ),                                    (4)
C                 (-Q   -A'),
C
C     for continuous-time equation, and
C                    -1        -1
C                 ( A         A  *G   )
C            H =  (   -1          -1  ),                         (5)
C                 (Q*A    A' + Q*A  *G)
C                                                            -1
C     for discrete-time equation, respectively, where G = B*R  *B'.
C     The assumptions guarantee that H in (4) has no pure imaginary
C     eigenvalues, and H in (5) has no eigenvalues on the unit circle.
C     If Y is an N-by-N matrix then there exists an orthogonal matrix U
C     such that U'*Y*U is an upper quasi-triangular matrix. Moreover, U
C     can be chosen so that the 2-by-2 and 1-by-1 diagonal blocks
C     (corresponding to the complex conjugate eigenvalues and real
C     eigenvalues respectively) appear in any desired order. This is the
C     ordered real Schur form. Thus, we can find an orthogonal
C     similarity transformation U which puts (4) or (5) in ordered real
C     Schur form
C
C            U'*H*U = S = (S(1,1)  S(1,2))
C                         (  0     S(2,2))
C
C     where S(i,j) is an N-by-N matrix and the eigenvalues of S(1,1)
C     have negative real parts in case of (4), or moduli greater than
C     one in case of (5). If U is conformably partitioned into four
C     N-by-N blocks
C
C               U = (U(1,1)  U(1,2))
C                   (U(2,1)  U(2,2))
C
C     with respect to the assumptions we then have
C     (a) U(1,1) is invertible and X = U(2,1)*inv(U(1,1)) solves (1),
C         (2), or (3) with X = X' and non-negative definite;
C     (b) the eigenvalues of S(1,1) (if DICO = 'C') or S(2,2) (if
C         DICO = 'D') are equal to the eigenvalues of optimal system
C         (the 'closed-loop' spectrum).
C
C     [A,B] is stabilizable if there exists a matrix F such that (A-BF)
C     is stable. [E,A] is detectable if [A',E'] is stabilizable.
C
C     REFERENCES
C
C     [1] Laub, A.J.
C         A Schur Method for Solving Algebraic Riccati equations.
C         IEEE Trans. Auto. Contr., AC-24, pp. 913-921, 1979.
C
C     [2] Wonham, W.M.
C         On a matrix Riccati equation of stochastic control.
C         SIAM J. Contr., 6, pp. 681-697, 1968.
C
C     [3] Sima, V.
C         Algorithms for Linear-Quadratic Optimization.
C         Pure and Applied Mathematics: A Series of Monographs and
C         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     FURTHER COMMENTS
C
C     To obtain a stabilizing solution of the algebraic Riccati
C     equation for DICO = 'D', set SORT = 'U', if HINV = 'D', or set
C     SORT = 'S', if HINV = 'I'.
C
C     The routine can also compute the anti-stabilizing solutions of
C     the algebraic Riccati equations, by specifying
C         SORT = 'U' if DICO = 'D' and HINV = 'I', or DICO = 'C', or
C         SORT = 'S' if DICO = 'D' and HINV = 'D'.
C
C     Usually, the combinations HINV = 'D' and SORT = 'U', or HINV = 'I'
C     and SORT = 'U', will be faster then the other combinations [3].
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB02AD by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, March 1982.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2002.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, discrete-time system, optimal regulator, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, HINV, SCAL, SORT, UPLO
      INTEGER           INFO, LDA, LDG, LDQ, LDS, LDU, LDWORK, N
      DOUBLE PRECISION  RCOND
C     .. Array Arguments ..
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), Q(LDQ,*),
     $                  S(LDS,*), U(LDU,*), WR(*), WI(*)
C     .. Local Scalars ..
      LOGICAL           DISCR, LHINV, LSCAL, LSORT, LUPLO
      INTEGER           I, IERR, ISCL, N2, NP1, NROT
      DOUBLE PRECISION  GNORM, QNORM, RCONDA, UNORM, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME, SB02MR, SB02MS, SB02MV, SB02MW
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANSY
      EXTERNAL          DLAMCH, DLANGE, DLANSY, LSAME, SB02MR, SB02MS,
     $                  SB02MV, SB02MW
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGECON, DGEES, DGETRF, DGETRS,
     $                  DLACPY, DLASCL, DLASET, DSCAL, DSWAP, SB02MU,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C     .. Executable Statements ..
C
      INFO = 0
      N2  = N + N
      NP1 = N + 1
      DISCR = LSAME( DICO, 'D' )
      LSCAL = LSAME( SCAL, 'G' )
      LSORT = LSAME( SORT, 'S' )
      LUPLO = LSAME( UPLO, 'U' )
      IF ( DISCR ) LHINV = LSAME( HINV, 'D' )
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
      ELSE IF( .NOT.LSCAL .AND. .NOT.LSAME( SCAL, 'N' ) ) THEN
         INFO = -4
      ELSE IF( .NOT.LSORT .AND. .NOT.LSAME( SORT, 'U' ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDS.LT.MAX( 1, N2 ) ) THEN
         INFO = -17
      ELSE IF( LDU.LT.MAX( 1, N2 ) ) THEN
         INFO = -19
      ELSE IF( ( .NOT.DISCR .AND. LDWORK.LT.MAX( 2, 6*N ) ) .OR.
     $         (      DISCR .AND. LDWORK.LT.MAX( 3, 6*N ) ) ) THEN
         INFO = -22
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB02MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         RCOND = ONE
         DWORK(1) = ONE
         DWORK(2) = ONE
         IF ( DISCR ) DWORK(3) = ONE
         RETURN
      END IF
C
      IF ( LSCAL ) THEN
C
C        Compute the norms of the matrices Q and G.
C
         QNORM = DLANSY( '1-norm', UPLO, N, Q, LDQ, DWORK )
         GNORM = DLANSY( '1-norm', UPLO, N, G, LDG, DWORK )
      END IF
C
C     Initialise the Hamiltonian or symplectic matrix associated with
C     the problem.
C     Workspace:  need   1          if DICO = 'C';
C                        max(2,4*N) if DICO = 'D';
C                 prefer larger if DICO = 'D'.
C
      CALL SB02MU( DICO, HINV, UPLO, N, A, LDA, G, LDG, Q, LDQ, S, LDS,
     $             IWORK, DWORK, LDWORK, INFO )
      IF ( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
      WRKOPT = DWORK(1)
      IF ( DISCR ) RCONDA = DWORK(2)
C
      ISCL = 0
      IF ( LSCAL ) THEN
C
C        Scale the Hamiltonian or symplectic matrix.
C
         IF( QNORM.GT.GNORM .AND. GNORM.GT.ZERO ) THEN
            CALL DLASCL( 'G', 0, 0, QNORM, GNORM, N, N, S(NP1,1), N2,
     $                   IERR )
            CALL DLASCL( 'G', 0, 0, GNORM, QNORM, N, N, S(1,NP1), N2,
     $                   IERR )
            ISCL = 1
         END IF
      END IF
C
C     Find the ordered Schur factorization of S,   S = U*H*U'.
C     Workspace:  need   6*N;
C                 prefer larger.
C
      IF ( .NOT.DISCR ) THEN
         IF ( LSORT ) THEN
            CALL DGEES( 'Vectors', 'Sorted', SB02MV, N2, S, LDS, NROT,
     $                  WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         ELSE
            CALL DGEES( 'Vectors', 'Sorted', SB02MR, N2, S, LDS, NROT,
     $                  WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         END IF
      ELSE
         IF ( LSORT ) THEN
            CALL DGEES( 'Vectors', 'Sorted', SB02MW, N2, S, LDS, NROT,
     $                  WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         ELSE
            CALL DGEES( 'Vectors', 'Sorted', SB02MS, N2, S, LDS, NROT,
     $                  WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         END IF
         IF ( LHINV ) THEN
            CALL DSWAP( N, WR, 1, WR(NP1), 1 )
            CALL DSWAP( N, WI, 1, WI(NP1), 1 )
         END IF
      END IF
      IF ( INFO.GT.N2 ) THEN
         INFO = 3
      ELSE IF ( INFO.GT.0 ) THEN
         INFO = 2
      ELSE IF ( NROT.NE.N ) THEN
         INFO = 4
      END IF
      IF ( INFO.NE.0 )
     $   RETURN
C
      WRKOPT = MAX( WRKOPT, DWORK(1) )
C
C     Check if U(1,1) is singular.  Use the (2,1) block of S as a
C     workspace for factoring U(1,1).
C
      UNORM = DLANGE( '1-norm', N, N, U, LDU, DWORK )
C
      CALL DLACPY( 'Full', N, N, U, LDU, S(NP1,1), LDS )
      CALL DGETRF( N, N, S(NP1,1), LDS, IWORK, INFO )
C
      IF ( INFO.GT.0 ) THEN
C
C        Singular matrix.  Set INFO and RCOND for error return.
C
         INFO  = 5
         RCOND = ZERO
         GO TO 100
      END IF
C
C     Estimate the reciprocal condition of U(1,1).
C     Workspace: 6*N.
C
      CALL DGECON( '1-norm', N, S(NP1,1), LDS, UNORM, RCOND,
     $             DWORK, IWORK(NP1), INFO )
C
      IF ( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
C
C        Nearly singular matrix.  Set INFO for error return.
C
         INFO = 5
         RETURN
      END IF
C
C     Transpose U(2,1) in Q and compute the solution.
C
      DO 60 I = 1, N
         CALL DCOPY( N, U(NP1,I), 1, Q(I,1), LDQ )
   60 CONTINUE
C
      CALL DGETRS( 'Transpose', N, N, S(NP1,1), LDS, IWORK, Q, LDQ,
     $             INFO )
C
C     Set S(2,1) to zero.
C
      CALL DLASET( 'Full', N, N, ZERO, ZERO, S(NP1,1), LDS )
C
C     Make sure the solution matrix X is symmetric.
C
      DO 80 I = 1, N - 1
         CALL DAXPY( N-I, ONE, Q(I,I+1), LDQ, Q(I+1,I), 1 )
         CALL DSCAL( N-I, HALF, Q(I+1,I), 1 )
         CALL DCOPY( N-I, Q(I+1,I), 1, Q(I,I+1), LDQ )
   80 CONTINUE
C
      IF( LSCAL ) THEN
C
C        Undo scaling for the solution matrix.
C
         IF( ISCL.EQ.1 )
     $      CALL DLASCL( 'G', 0, 0, GNORM, QNORM, N, N, Q, LDQ, IERR )
      END IF
C
C     Set the optimal workspace, the scaling factor, and reciprocal
C     condition number (if any).
C
      DWORK(1) = WRKOPT
  100 CONTINUE
      IF( ISCL.EQ.1 ) THEN
         DWORK(2) = QNORM / GNORM
      ELSE
         DWORK(2) = ONE
      END IF
      IF ( DISCR ) DWORK(3) = RCONDA
C
      RETURN
C *** Last line of SB02MD ***
      END
