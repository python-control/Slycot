      SUBROUTINE SB02PD( JOB, TRANA, UPLO, N, A, LDA, G, LDG, Q, LDQ, X,
     $                   LDX, RCOND, FERR, WR, WI, IWORK, DWORK, LDWORK,
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
C     To solve the real continuous-time matrix algebraic Riccati
C     equation
C
C        op(A)'*X + X*op(A) + Q - X*G*X = 0,
C
C     where op(A) = A or A' = A**T and G, Q are symmetric (G = G**T,
C     Q = Q**T). The matrices A, G and Q are N-by-N and the solution X
C     is an N-by-N symmetric matrix.
C
C     An error bound on the solution and a condition estimate are also
C     optionally provided.
C
C     It is assumed that the matrices A, G and Q are such that the
C     corresponding Hamiltonian matrix has N eigenvalues with negative
C     real parts.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the computation to be performed, as follows:
C             = 'X':  Compute the solution only;
C             = 'A':  Compute all: the solution, reciprocal condition
C                     number, and the error bound.
C
C     TRANA   CHARACTER*1
C             Specifies the option op(A):
C             = 'N':  op(A) = A    (No transpose);
C             = 'T':  op(A) = A**T (Transpose);
C             = 'C':  op(A) = A**T (Conjugate transpose = Transpose).
C
C     UPLO    CHARACTER*1
C             Specifies which triangle of the matrices G and Q is
C             stored, as follows:
C             = 'U':  Upper triangles of G and Q are stored;
C             = 'L':  Lower triangles of G and Q are stored.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, Q, and X.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             coefficient matrix A of the equation.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     G       (input) DOUBLE PRECISION array, dimension (LDG,N)
C             If UPLO = 'U', the leading N-by-N upper triangular part of
C             this array must contain the upper triangular part of the
C             matrix G.
C             If UPLO = 'L', the leading N-by-N lower triangular part of
C             this array must contain the lower triangular part of the
C             matrix G.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= max(1,N).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N)
C             If UPLO = 'U', the leading N-by-N upper triangular part of
C             this array must contain the upper triangular part of the
C             matrix Q.
C             If UPLO = 'L', the leading N-by-N lower triangular part of
C             this array must contain the lower triangular part of the
C             matrix Q.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= max(1,N).
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N)
C             If INFO = 0, INFO = 2, or INFO = 4, the leading N-by-N
C             part of this array contains the symmetric solution matrix
C             X of the algebraic Riccati equation.
C
C     LDX     INTEGER
C             The leading dimension of the array X.  LDX >= max(1,N).
C
C     RCOND   (output) DOUBLE PRECISION
C             If JOB = 'A', the estimate of the reciprocal condition
C             number of the Riccati equation.
C
C     FERR    (output) DOUBLE PRECISION
C             If JOB = 'A', the estimated forward error bound for the
C             solution X. If XTRUE is the true solution, FERR bounds the
C             magnitude of the largest entry in (X - XTRUE) divided by
C             the magnitude of the largest entry in X.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             If JOB = 'A' and TRANA = 'N', WR and WI contain the real
C             and imaginary parts, respectively, of the eigenvalues of
C             the matrix A - G*X, i.e., the closed-loop system poles.
C             If JOB = 'A' and TRANA = 'T' or 'C', WR and WI contain the
C             real and imaginary parts, respectively, of the eigenvalues
C             of the matrix A - X*G, i.e., the closed-loop system poles.
C             If JOB = 'X', these arrays are not referenced.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK), where
C             LIWORK >= 2*N,          if JOB = 'X';
C             LIWORK >= max(2*N,N*N), if JOB = 'A'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or INFO = 2, DWORK(1) contains the
C             optimal value of LDWORK. If JOB = 'A', then DWORK(2:N*N+1)
C             and DWORK(N*N+2:2*N*N+1) contain a real Schur form of the
C             closed-loop system matrix, Ac = A - G*X (if TRANA = 'N')
C             or Ac = A - X*G (if TRANA = 'T' or 'C'), and the
C             orthogonal matrix which reduced Ac to real Schur form,
C             respectively.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 4*N*N + 8*N + 1,               if JOB = 'X';
C             LDWORK >= max( 4*N*N + 8*N, 6*N*N ) + 1, if JOB = 'A'.
C             For good performance, LDWORK should be larger, e.g.,
C             LDWORK >= 4*N*N + 6*N +( 2*N+1 )*NB,     if JOB = 'X',
C             where NB is the optimal blocksize.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the Hamiltonian matrix has eigenvalues on the
C                   imaginary axis, so the solution and error bounds
C                   could not be computed;
C             = 2:  the iteration for the matrix sign function failed to
C                   converge after 50 iterations, but an approximate
C                   solution and error bounds (if JOB = 'A') have been
C                   computed;
C             = 3:  the system of linear equations for the solution is
C                   singular to working precision, so the solution and
C                   error bounds could not be computed;
C             = 4:  the matrix A-G*X (or A-X*G) cannot be reduced to
C                   Schur canonical form and condition number estimate
C                   and forward error estimate have not been computed.
C
C     METHOD
C
C     The Riccati equation is solved by the matrix sign function
C     approach [1], [2], implementing a scaling which enhances the
C     numerical stability [4].
C
C     REFERENCES
C
C     [1] Bai, Z., Demmel, J., Dongarra, J., Petitet, A., Robinson, H.,
C         and Stanley, K.
C         The spectral decomposition of nonsymmetric matrices on
C         distributed memory parallel computers.
C         SIAM J. Sci. Comput., vol. 18, pp. 1446-1461, 1997.
C
C     [2] Byers, R., He, C., and Mehrmann, V.
C         The matrix sign function method and the computation of
C         invariant subspaces.
C         SIAM J. Matrix Anal. Appl., vol. 18, pp. 615-632, 1997.
C
C     [3] Higham, N.J.
C         Perturbation theory and backward error for AX-XB=C.
C         BIT, vol. 33, pp. 124-136, 1993.
C
C     [4] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V.,
C         DGRSVX and DMSRIC: Fortran 77 subroutines for solving
C         continuous-time matrix algebraic Riccati equations with
C         condition and accuracy estimates.
C         Preprint SFB393/98-16, Fak. f. Mathematik, Technical
C         University Chemnitz, May 1998.
C
C     NUMERICAL ASPECTS
C
C     The solution accuracy can be controlled by the output parameter
C     FERR.
C
C     FURTHER COMMENTS
C
C     The condition number of the Riccati equation is estimated as
C
C     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) +
C                 norm(Pi)*norm(G) ) / norm(X),
C
C     where Omega, Theta and Pi are linear operators defined by
C
C     Omega(W) = op(Ac)'*W + W*op(Ac),
C     Theta(W) = inv(Omega(op(W)'*X + X*op(W))),
C        Pi(W) = inv(Omega(X*W*X)),
C
C     and the matrix Ac (the closed-loop system matrix) is given by
C        Ac = A - G*X, if TRANA = 'N', or
C        Ac = A - X*G, if TRANA = 'T' or 'C'.
C
C     The program estimates the quantities
C
C     sep(op(Ac),-op(Ac)') = 1 / norm(inv(Omega)),
C
C     norm(Theta) and norm(Pi) using 1-norm condition estimator.
C
C     The forward error bound is estimated using a practical error bound
C     similar to the one proposed in [3].
C
C     CONTRIBUTOR
C
C     P. Petkov, Tech. University of Sofia, March 2000.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, June 2000.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, continuous-time system,
C     optimal control, optimal regulator.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 50 )
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
     $                     TWO  = 2.0D+0, TEN  = 10.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          JOB, TRANA, UPLO
      INTEGER            INFO, LDA, LDG, LDQ, LDWORK, LDX, N
      DOUBLE PRECISION   FERR, RCOND
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), DWORK( * ), G( LDG, * ),
     $                   Q( LDQ, * ), WI( * ), WR( * ), X( LDX, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            ALL, LOWER, NOTRNA
      CHARACTER          EQUED, LOUP
      INTEGER            I, IAF, IB, IBR, IC, IFR, IJ, IJ1, IJ2, INFO2,
     $                   INI, IR, ISCL, ISV, IT, ITAU, ITER, IU, IWRK,
     $                   J, JI, LWAMAX, MINWRK, N2, SDIM
      DOUBLE PRECISION   CONV, GNORM2, EPS, HNORM, HINNRM, QNORM2,
     $                   SCALE, SEP, TEMP, TOL
C     ..
C     .. Local Arrays ..
      LOGICAL            BWORK( 1 )
C     ..
C     .. External Functions ..
      LOGICAL            LSAME, SELECT
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           DLAMCH, DLANSY, ILAENV, LSAME, SELECT
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEES, DGEQP3, DGESVX, DLACPY, DLASCL,
     $                   DLASET, DORMQR, DSCAL, DSWAP, DSYMM, DSYTRF,
     $                   DSYTRI, MA02AD, MA02ED, SB02QD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      ALL    = LSAME( JOB,   'A' )
      NOTRNA = LSAME( TRANA, 'N' )
      LOWER  = LSAME( UPLO,  'L' )
C
      INFO = 0
      IF( .NOT.ALL .AND. .NOT.LSAME( JOB, 'X' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANA, 'T' ) .AND.
     $         .NOT.LSAME( TRANA, 'C' ) .AND. .NOT.NOTRNA ) THEN
         INFO = -2
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE
C
C        Compute workspace.
C
         IF( ALL ) THEN
            MINWRK = MAX( 4*N*N + 8*N + 1, 6*N*N )
         ELSE
            MINWRK = 4*N*N + 8*N + 1
         END IF
         IF( LDWORK.LT.MINWRK ) THEN
            INFO = -19
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB02PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         IF( ALL ) THEN
            RCOND = ONE
            FERR  = ZERO
         END IF
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Set tol.
C
      EPS = DLAMCH( 'P' )
      TOL = TEN*DBLE( N )*EPS
C
C     Compute the square-roots of the norms of the matrices Q and G .
C
      QNORM2 = SQRT( DLANSY( '1', UPLO, N, Q, LDQ, DWORK ) )
      GNORM2 = SQRT( DLANSY( '1', UPLO, N, G, LDG, DWORK ) )
C
      N2 = 2*N
C
C     Construct the lower (if UPLO = 'L') or upper (if UPLO = 'U')
C     triangle of the symmetric block-permuted Hamiltonian matrix.
C     During iteration, both the current iterate corresponding to the
C     Hamiltonian matrix, and its inverse are needed. To reduce the
C     workspace length, the transpose of the triangle specified by UPLO
C     of the current iterate H is saved in the opposite triangle,
C     suitably shifted with one column, and then the inverse of H
C     overwrites H. The triangles of the saved iterate and its inverse
C     are stored together in an 2*N-by-(2*N+1) matrix. For instance, if
C     UPLO = 'U', then the upper triangle is built starting from the
C     location 2*N+1 of the array DWORK, so that its transpose can be
C     stored in the lower triangle of DWORK.
C     Workspace: need   4*N*N,        if UPLO = 'L';
C                       4*N*N + 2*N,  if UPLO = 'U'.
C
      IF ( LOWER ) THEN
         INI  = 0
         ISV  = N2
         LOUP = 'U'
C
         DO 40 J = 1, N
            IJ = ( J - 1 )*N2 + J
C
            DO 10 I = J, N
               DWORK(IJ) = -Q(I,J)
               IJ = IJ + 1
   10       CONTINUE
C
            IF( NOTRNA ) THEN
C
               DO 20 I = 1, N
                  DWORK( IJ ) = -A( I, J )
                  IJ = IJ + 1
   20          CONTINUE
C
            ELSE
C
               DO 30 I = 1, N
                  DWORK( IJ ) = -A( J, I )
                  IJ = IJ + 1
   30          CONTINUE
C
            END IF
   40    CONTINUE
C
         DO 60 J = 1, N
            IJ = ( N + J - 1 )*N2 + N + J
C
            DO 50 I = J, N
               DWORK( IJ ) = G( I, J )
               IJ = IJ + 1
   50       CONTINUE
C
   60    CONTINUE
C
      ELSE
         INI  = N2
         ISV  = 0
         LOUP = 'L'
C
         DO 80 J = 1, N
            IJ = J*N2 + 1
C
            DO 70 I = 1, J
               DWORK(IJ) = -Q(I,J)
               IJ = IJ + 1
   70       CONTINUE
C
   80    CONTINUE
C
         DO 120 J = 1, N
            IJ = ( N + J )*N2 + 1
C
            IF( NOTRNA ) THEN
C
               DO 90 I = 1, N
                  DWORK( IJ ) = -A( J, I )
                  IJ = IJ + 1
   90          CONTINUE
C
            ELSE
C
               DO 100 I = 1, N
                  DWORK( IJ ) = -A( I, J )
                  IJ = IJ + 1
  100          CONTINUE
C
            END IF
C
            DO 110 I = 1, J
               DWORK( IJ ) = G( I, J )
               IJ = IJ + 1
  110       CONTINUE
C
  120    CONTINUE
C
      END IF
C
C     Block-scaling.
C
      ISCL = 0
      IF( QNORM2.GT.GNORM2 .AND. GNORM2.GT.ZERO ) THEN
         CALL DLASCL( UPLO, 0, 0, QNORM2, GNORM2, N, N, DWORK( INI+1 ),
     $                N2, INFO2 )
         CALL DLASCL( UPLO, 0, 0, GNORM2, QNORM2, N, N,
     $                DWORK( N2*N+N+INI+1 ), N2, INFO2 )
         ISCL = 1
      END IF
C
C     Workspace usage.
C
      ITAU = N2*N2
      IWRK = ITAU + N2
C
      LWAMAX = N2*ILAENV( 1, 'DSYTRF', UPLO, N2, -1, -1, -1 )
C
C     Compute the matrix sign function.
C
      DO 230 ITER = 1, MAXIT
C
C        Save the transpose of the corresponding triangle of the
C        current iterate in the free locations of the shifted opposite
C        triangle.
C        Workspace: need   4*N*N + 2*N.
C
         IF( LOWER ) THEN
C
            DO 130 I = 1, N2
               CALL DCOPY( I, DWORK( I ), N2, DWORK( I*N2+1 ), 1 )
  130       CONTINUE
C
         ELSE
C
            DO 140 I = 1, N2
               CALL DCOPY( I, DWORK( I*N2+1 ), 1, DWORK( I ), N2 )
  140       CONTINUE
C
         END IF
C
C        Store the norm of the Hamiltonian matrix.
C
         HNORM = DLANSY( 'F', UPLO, N2, DWORK( INI+1 ), N2, DWORK )
C
C        Compute the inverse of the block-permuted Hamiltonian matrix.
C        Workspace: need   4*N*N + 2*N + 1;
C                   prefer 4*N*N + 2*N + 2*N*NB.
C
         CALL DSYTRF( UPLO, N2, DWORK( INI+1 ), N2, IWORK,
     $                DWORK( IWRK+1 ), LDWORK-IWRK, INFO2 )
         IF( INFO2.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
C
C        Workspace: need   4*N*N + 4*N.
C
         CALL DSYTRI( UPLO, N2, DWORK( INI+1 ), N2, IWORK,
     $                DWORK( IWRK+1 ), INFO2 )
C
C        Block-permutation of the inverse matrix.
C
         IF( LOWER ) THEN
C
            DO 160 J = 1, N
               IJ2 = ( N + J - 1 )*N2 + N + J
C
               DO 150 IJ1 = ( J - 1 )*N2 + J, ( J - 1 )*N2 + N
                  TEMP = DWORK( IJ1 )
                  DWORK( IJ1 ) = -DWORK( IJ2 )
                  DWORK( IJ2 ) = -TEMP
                  IJ2 = IJ2 + 1
  150          CONTINUE
C
               CALL DSWAP( J-1, DWORK( N+J ), N2, DWORK( (J-1)*N2+N+1 ),
     $                     1 )
  160       CONTINUE
C
         ELSE
C
            DO 180 J = 1, N
               IJ2 = ( N + J )*N2 + N + 1
C
               DO 170 IJ1 = J*N2 + 1, J*N2 + J
                  TEMP = DWORK( IJ1 )
                  DWORK( IJ1 ) = -DWORK( IJ2 )
                  DWORK( IJ2 ) = -TEMP
                  IJ2 = IJ2 + 1
  170          CONTINUE
C
               CALL DSWAP( J-1, DWORK( (N+1)*N2+J ), N2,
     $                     DWORK( (N+J)*N2+1 ), 1 )
  180       CONTINUE
C
         END IF
C
C        Scale the Hamiltonian matrix and its inverse and compute
C        the next iterate.
C
         HINNRM = DLANSY( 'F', UPLO, N2, DWORK( INI+1 ), N2, DWORK )
         SCALE  = SQRT( HINNRM / HNORM )
C
         IF( LOWER ) THEN
C
            DO 200 J = 1, N2
               JI = ( J - 1 )*N2 + J
C
               DO 190 IJ = JI, J*N2
                  JI = JI + N2
                  DWORK( IJ ) = ( DWORK( IJ ) / SCALE +
     $                            DWORK( JI )*SCALE ) / TWO
                  DWORK( JI ) =   DWORK( JI ) - DWORK( IJ )
  190          CONTINUE
C
  200       CONTINUE
C
         ELSE
C
            DO 220 J = 1, N2
               JI = J
C
               DO 210 IJ = J*N2 + 1, J*N2 + J
                  DWORK( IJ ) = ( DWORK( IJ ) / SCALE +
     $                            DWORK( JI )*SCALE ) / TWO
                  DWORK( JI ) =   DWORK( JI ) - DWORK( IJ )
                  JI = JI + N2
  210          CONTINUE
C
  220       CONTINUE
C
         END IF
C
C        Test for convergence.
C
         CONV = DLANSY( 'F', LOUP, N2, DWORK( ISV+1 ), N2, DWORK )
         IF( CONV.LE.TOL*HNORM ) GO TO 240
  230 CONTINUE
C
C     No convergence after MAXIT iterations, but an approximate solution
C     has been found.
C
      INFO = 2
C
  240 CONTINUE
C
C     If UPLO = 'U', shift the upper triangle one column to the left.
C
      IF( .NOT.LOWER )
     $   CALL DLACPY( 'U', N2, N2, DWORK( INI+1 ), N2, DWORK, N2 )
C
C     Divide the triangle elements by -2 and then fill-in the other
C     triangle by symmetry.
C
      IF( LOWER ) THEN
C
         DO 250 I = 1, N2
            CALL DSCAL( N2-I+1, -HALF, DWORK( (I-1)*N2+I ), 1 )
  250    CONTINUE
C
      ELSE
C
         DO 260 I = 1, N2
            CALL DSCAL( I, -HALF, DWORK( (I-1)*N2+1 ), 1 )
  260    CONTINUE
C
      END IF
      CALL MA02ED( UPLO, N2, DWORK, N2 )
C
C     Back block-permutation.
C
      DO 280 J = 1, N2
C
         DO 270 I = ( J - 1 )*N2 + 1, ( J - 1 )*N2 + N
            TEMP = DWORK( I )
            DWORK( I )   = -DWORK( I+N )
            DWORK( I+N ) = TEMP
  270    CONTINUE
C
  280 CONTINUE
C
C     Compute the QR decomposition of the projector onto the stable
C     invariant subspace.
C     Workspace: need   4*N*N + 8*N + 1.
C                prefer 4*N*N + 6*N + ( 2*N+1 )*NB.
C
      DO 290 I = 1, N2
         IWORK( I ) = 0
         DWORK( ( I-1 )*N2 + I ) = DWORK( ( I-1 )*N2 + I ) + HALF
  290 CONTINUE
C
      CALL DGEQP3( N2, N2, DWORK, N2, IWORK, DWORK( ITAU+1 ),
     $             DWORK( IWRK+1 ), LDWORK-IWRK, INFO2 )
      LWAMAX = MAX( INT( DWORK( IWRK+1 ) ), LWAMAX )
C
C     Accumulate the orthogonal transformations. Note that only the
C     first N columns of the array DWORK, returned by DGEQP3, are
C     needed, so that the last N columns of DWORK are used to get the
C     orthogonal basis for the stable invariant subspace.
C     Workspace: need   4*N*N + 3*N.
C                prefer 4*N*N + 2*N + N*NB.
C
      IB  = N*N
      IAF = N2*N
      CALL DLASET( 'F', N2, N, ZERO, ONE, DWORK( IAF+1 ), N2 )
      CALL DORMQR( 'L', 'N', N2, N, N, DWORK, N2, DWORK( ITAU+1 ),
     $             DWORK( IAF+1 ), N2, DWORK( IWRK+1 ), LDWORK-IWRK,
     $             INFO2 )
      LWAMAX = IWRK + MAX( INT( DWORK( IWRK+1 ) ), LWAMAX )
C
C     Store the matrices V11 and V21' .
C
      CALL DLACPY( 'F', N, N, DWORK( IAF+1 ), N2, DWORK, N )
      CALL MA02AD( 'F', N, N, DWORK( IAF+N+1 ), N2, DWORK( IB+1 ), N )
C
      IR   = IAF + IB
      IC   = IR  + N
      IFR  = IC  + N
      IBR  = IFR + N
      IWRK = IBR + N
C
C     Compute the solution matrix X .
C     Workspace: need   3*N*N + 8*N.
C
      CALL DGESVX( 'E', 'T', N, N, DWORK, N, DWORK( IAF+1 ), N,
     $             IWORK, EQUED, DWORK( IR+1 ), DWORK( IC+1 ),
     $             DWORK( IB+1 ), N, X, LDX, RCOND, DWORK( IFR+1 ),
     $             DWORK( IBR+1 ), DWORK( IWRK+1 ), IWORK( N+1 ),
     $             INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 3
         RETURN
      END IF
C
C     Symmetrize the solution.
C
      DO 310 I = 1, N - 1
C
         DO 300 J = I + 1, N
            TEMP = ( X( I, J ) + X( J, I ) ) / TWO
            X( I, J ) = TEMP
            X( J, I ) = TEMP
  300    CONTINUE
C
  310 CONTINUE
C
C     Undo scaling for the solution matrix.
C
      IF( ISCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, GNORM2, QNORM2, N, N, X, LDX, INFO2 )
      END IF
C
      IF( ALL ) THEN
C
C        Compute the estimates of the reciprocal condition number and
C        error bound.
C        Workspace usage.
C
         IT   = 1
         IU   = IT + N*N
         IWRK = IU + N*N
C
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IT+1 ), N )
         IF( NOTRNA ) THEN
C
C           Compute Ac = A-G*X .
C
            CALL DSYMM( 'L', UPLO, N, N, -ONE, G, LDG, X, LDX, ONE,
     $                  DWORK( IT+1 ), N )
         ELSE
C
C           Compute Ac = A-X*G .
C
            CALL DSYMM( 'R', UPLO, N, N, -ONE, G, LDG, X, LDX, ONE,
     $                  DWORK( IT+1 ), N )
         END IF
C
C        Compute the Schur factorization of Ac .
C        Workspace: need   2*N*N + 5*N + 1;
C                   prefer larger.
C
         CALL DGEES( 'V', 'N', SELECT, N, DWORK( IT+1 ), N, SDIM, WR,
     $               WI, DWORK( IU+1 ), N, DWORK( IWRK+1 ), LDWORK-IWRK,
     $               BWORK, INFO2 )
         IF( INFO2.GT.0 ) THEN
            INFO = 4
            RETURN
         END IF
         LWAMAX = IWRK + MAX( INT( DWORK( IWRK+1 ) ), LWAMAX )
C
C        Estimate the reciprocal condition number and the forward error.
C        Workspace: need   6*N*N + 1;
C                   prefer larger.
C
         CALL SB02QD( 'B', 'F', TRANA, UPLO, 'O', N, A, LDA,
     $                DWORK( IT+1 ), N, DWORK( IU+1 ), N, G, LDG, Q,
     $                LDQ, X, LDX, SEP, RCOND, FERR, IWORK,
     $                DWORK( IWRK+1 ), LDWORK-IWRK, INFO2 )
         LWAMAX = IWRK + MAX( INT( DWORK( IWRK+1 ) ), LWAMAX )
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB02PD
      END
