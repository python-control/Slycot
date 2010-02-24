      SUBROUTINE MB04ZD( COMPU, N, A, LDA, QG, LDQG, U, LDU, DWORK, INFO
     $                 )
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
C     To transform a Hamiltonian matrix
C
C               ( A   G  )
C           H = (      T )                                           (1)
C               ( Q  -A  )
C
C     into a square-reduced Hamiltonian matrix
C
C                ( A'  G'  )
C           H' = (       T )                                         (2)
C                ( Q' -A'  )
C                                                                 T
C     by an orthogonal symplectic similarity transformation H' = U H U,
C     where
C               (  U1   U2 )
C           U = (          ).                                        (3)
C               ( -U2   U1 )
C                                                              T
C     The square-reduced Hamiltonian matrix satisfies Q'A' - A' Q' = 0,
C     and
C
C           2       T     2     ( A''   G''  )
C         H'  :=  (U  H U)   =  (          T ).
C                               ( 0     A''  )
C
C     In addition, A'' is upper Hessenberg and G'' is skew symmetric.
C     The square roots of the eigenvalues of A'' = A'*A' + G'*Q' are the
C     eigenvalues of H.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPU   CHARACTER*1
C             Indicates whether the orthogonal symplectic similarity
C             transformation matrix U in (3) is returned or
C             accumulated into an orthogonal symplectic matrix, or if
C             the transformation matrix is not required, as follows:
C             = 'N':         U is not required;
C             = 'I' or 'F':  on entry, U need not be set;
C                            on exit, U contains the orthogonal
C                            symplectic matrix U from (3);
C             = 'V' or 'A':  the orthogonal symplectic similarity
C                            transformations are accumulated into U;
C                            on input, U must contain an orthogonal
C                            symplectic matrix S;
C                            on exit, U contains S*U with U from (3).
C             See the description of U below for details.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, the leading N-by-N part of this array must
C             contain the upper left block A of the Hamiltonian matrix H
C             in (1).
C             On output, the leading N-by-N part of this array contains
C             the upper left block A' of the square-reduced Hamiltonian
C             matrix H' in (2).
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C             (LDQG,N+1)
C             On input, the leading N-by-N lower triangular part of this
C             array must contain the lower triangle of the lower left
C             symmetric block Q of the Hamiltonian matrix H in (1), and
C             the N-by-N upper triangular part of the submatrix in the
C             columns 2 to N+1 of this array must contain the upper
C             triangle of the upper right symmetric block G of H in (1).
C             So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C             and G(i,j) = G(j,i) is stored in QG(j,i+1).
C             On output, the leading N-by-N lower triangular part of
C             this array contains the lower triangle of the lower left
C             symmetric block Q', and the N-by-N upper triangular part
C             of the submatrix in the columns 2 to N+1 of this array
C             contains the upper triangle of the upper right symmetric
C             block G' of the square-reduced Hamiltonian matrix H'
C             in (2).
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= MAX(1,N).
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,2*N)
C             If COMPU = 'N', then this array is not referenced.
C             If COMPU = 'I' or 'F', then the input contents of this
C             array are not specified.  On output, the leading
C             N-by-(2*N) part of this array contains the first N rows
C             of the orthogonal symplectic matrix U in (3).
C             If COMPU = 'V' or 'A', then, on input, the leading
C             N-by-(2*N) part of this array must contain the first N
C             rows of an orthogonal symplectic matrix S. On output, the
C             leading N-by-(2*N) part of this array contains the first N
C             rows of the product S*U where U is the orthogonal
C             symplectic matrix from (3).
C             The storage scheme implied by (3) is used for orthogonal
C             symplectic matrices, i.e., only the first N rows are
C             stored, as they contain all relevant information.
C
C     LDU     INTEGER
C             The leading dimension of the array U.
C             LDU >= MAX(1,N), if COMPU <> 'N';
C             LDU >= 1,        if COMPU =  'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, then the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The Hamiltonian matrix H is transformed into a square-reduced
C     Hamiltonian matrix H' using the implicit version of Van Loan's
C     method as proposed in [1,2,3].
C
C     REFERENCES
C
C     [1] Van Loan, C. F.
C         A Symplectic Method for Approximating All the Eigenvalues of
C         a Hamiltonian Matrix.
C         Linear Algebra and its Applications, 61, pp. 233-251, 1984.
C
C     [2] Byers, R.
C         Hamiltonian and Symplectic Algorithms for the Algebraic
C         Riccati Equation.
C         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983.
C
C     [3] Benner, P., Byers, R., and Barth, E.
C         Fortran 77 Subroutines for Computing the Eigenvalues of
C         Hamiltonian Matrices. I: The Square-Reduced Method.
C         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000.
C
C     NUMERICAL ASPECTS
C
C     This algorithm requires approximately 20*N**3 flops for
C     transforming H into square-reduced form. If the transformations
C     are required, this adds another 8*N**3 flops. The method is
C     strongly backward stable in the sense that if H' and U are the
C     computed square-reduced Hamiltonian and computed orthogonal
C     symplectic similarity transformation, then there is an orthogonal
C     symplectic matrix T and a Hamiltonian matrix M such that
C
C                  H T  =  T M
C
C        || T - U ||   <=  c1 * eps
C
C        || H' - M ||  <=  c2 * eps * || H ||
C
C     where c1, c2 are modest constants depending on the dimension N and
C     eps is the machine precision.
C
C     Eigenvalues computed by explicitly forming the upper Hessenberg
C     matrix  A'' = A'A' + G'Q', with A', G', and Q' as in (2), and
C     applying the Hessenberg QR iteration to A'' are exactly
C     eigenvalues of a perturbed Hamiltonian matrix H + E,  where
C
C        || E ||  <=  c3 * sqrt(eps) * || H ||,
C
C     and c3 is a modest constant depending on the dimension N and eps
C     is the machine precision.  Moreover, if the norm of H and an
C     eigenvalue lambda are of roughly the same magnitude, the computed
C     eigenvalue is essentially as accurate as the computed eigenvalue
C     from traditional methods.  See [1] or [2].
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany,
C     R. Byers, University of Kansas, Lawrence, USA, and
C     E. Barth, Kalamazoo College, Kalamazoo, USA,
C     Aug. 1998, routine DHASRD.
C     V. Sima, Research Institute for Informatics, Bucharest, Romania,
C     Oct. 1998, SLICOT Library version.
C
C     REVISIONS
C
C     May 2001, A. Varga, German Aeropsce Center, DLR Oberpfaffenhofen.
C     May 2009, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     Orthogonal transformation, (square-reduced) Hamiltonian matrix,
C     symplectic similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
C
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDQG, LDU, N
      CHARACTER         COMPU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), QG(LDQG,*), U(LDU,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  COSINE, SINE, TAU, TEMP, X, Y
      INTEGER           J
      LOGICAL           ACCUM, FORGET, FORM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1), T(2,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      LOGICAL           LSAME
      EXTERNAL          DDOT, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DLARFG, DLARFX, DLARTG,
     $                  DROT, DSYMV, DSYR2, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
      INFO   = 0
      ACCUM  = LSAME( COMPU, 'A' ) .OR. LSAME( COMPU, 'V' )
      FORM   = LSAME( COMPU, 'F' ) .OR. LSAME( COMPU, 'I' )
      FORGET = LSAME( COMPU, 'N' )
C
      IF ( .NOT.ACCUM .AND. .NOT.FORM .AND. .NOT.FORGET ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQG.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDU.LT.1 .OR. ( .NOT.FORGET .AND. LDU.LT.MAX( 1, N ) ) )
     $      THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04ZD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Transform to square-reduced form.
C
      DO 10 J = 1, N - 1
C                         T
C        DWORK <- (Q*A - A *Q)(J+1:N,J).
C
         CALL DCOPY( J-1, QG(J,1), LDQG, DWORK(N+1), 1 )
         CALL DCOPY( N-J+1, QG(J,J), 1, DWORK(N+J), 1 )
         CALL DGEMV( 'Transpose', N, N-J, -ONE, A(1,J+1), LDA,
     $               DWORK(N+1), 1, ZERO, DWORK(J+1), 1 )
         CALL DGEMV( 'NoTranspose', N-J, J, ONE, QG(J+1,1), LDQG,
     $               A(1,J), 1, ONE, DWORK(J+1), 1 )
         CALL DSYMV( 'Lower', N-J, ONE, QG(J+1,J+1), LDQG, A(J+1,J), 1,
     $               ONE, DWORK(J+1), 1 )
C
C        Symplectic reflection to zero (H*H)((N+J+2):2N,J).
C
         CALL DLARFG( N-J, DWORK(J+1), DWORK(J+2), 1, TAU )
         Y = DWORK(J+1)
         DWORK(J+1) = ONE
C
         CALL DLARFX( 'Left', N-J, N, DWORK(J+1), TAU, A(J+1,1), LDA,
     $                DWORK(N+1) )
         CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, A(1,J+1), LDA,
     $                DWORK(N+1) )
C
         CALL DLARFX( 'Left', N-J, J, DWORK(J+1), TAU, QG(J+1,1), LDQG,
     $                DWORK(N+1) )
         CALL DSYMV( 'Lower', N-J, TAU, QG(J+1,J+1), LDQG, DWORK(J+1),
     $               1, ZERO, DWORK(N+J+1), 1 )
         CALL DAXPY( N-J, -TAU*DDOT( N-J, DWORK(N+J+1), 1, DWORK(J+1),
     $               1 )/TWO, DWORK(J+1), 1, DWORK(N+J+1), 1 )
         CALL DSYR2( 'Lower', N-J, -ONE, DWORK(J+1), 1, DWORK(N+J+1), 1,
     $               QG(J+1,J+1), LDQG )
C
         CALL DLARFX( 'Right', J, N-J, DWORK(J+1), TAU, QG(1,J+2), LDQG,
     $                DWORK(N+1) )
         CALL DSYMV( 'Upper', N-J, TAU, QG(J+1,J+2), LDQG, DWORK(J+1),
     $               1, ZERO, DWORK(N+J+1), 1 )
         CALL DAXPY( N-J, -TAU*DDOT( N-J, DWORK(N+J+1), 1, DWORK(J+1),
     $               1 )/TWO, DWORK(J+1), 1, DWORK(N+J+1), 1 )
         CALL DSYR2( 'Upper', N-J, -ONE, DWORK(J+1), 1, DWORK(N+J+1), 1,
     $               QG(J+1,J+2), LDQG )
C
         IF ( FORM ) THEN
C
C           Save reflection.
C
            CALL DCOPY( N-J, DWORK(J+1), 1, U(J+1,J), 1 )
            U(J+1,J) = TAU
C
         ELSE IF ( ACCUM ) THEN
C
C           Accumulate reflection.
C
            CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, U(1,J+1),
     $                   LDU, DWORK(N+1) )
            CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, U(1,N+J+1),
     $                   LDU, DWORK(N+1) )
         END IF
C
C        (X,Y) := ((J+1,J),(N+J+1,J)) component of H*H.
C
         X = DDOT( J, QG(1,J+2), 1, QG(J,1), LDQG ) +
     $       DDOT( N-J, QG(J+1,J+2), LDQG, QG(J+1,J), 1 ) +
     $       DDOT( N, A(J+1,1), LDA, A(1,J), 1 )
C
C        Symplectic rotation to zero (H*H)(N+J+1,J).
C
         CALL DLARTG( X, Y, COSINE, SINE, TEMP )
C
         CALL DROT( J, A(J+1,1), LDA, QG(J+1,1), LDQG, COSINE, SINE )
         CALL DROT( J, A(1,J+1), 1, QG(1,J+2), 1, COSINE, SINE )
         IF( J.LT.N-1 ) THEN
            CALL DROT( N-J-1, A(J+1,J+2), LDA, QG(J+2,J+1), 1,
     $                 COSINE, SINE )
            CALL DROT( N-J-1, A(J+2,J+1), 1, QG(J+1,J+3), LDQG,
     $                 COSINE, SINE )
         END IF
C
         T(1,1) = A(J+1,J+1)
         T(1,2) = QG(J+1,J+2)
         T(2,1) = QG(J+1,J+1)
         T(2,2) = -T(1,1)
         CALL DROT( 2, T(1,1), 1, T(1,2), 1, COSINE, SINE )
         CALL DROT( 2, T(1,1), 2, T(2,1), 2, COSINE, SINE )
         A(J+1,J+1)  = T(1,1)
         QG(J+1,J+2) = T(1,2)
         QG(J+1,J+1) = T(2,1)
C
         IF ( FORM ) THEN
C
C           Save rotation.
C
            U(J,J)   = COSINE
            U(J,N+J) = SINE
C
         ELSE IF ( ACCUM ) THEN
C
C           Accumulate rotation.
C
            CALL DROT( N, U(1,J+1), 1, U(1,N+J+1), 1, COSINE, SINE )
         END IF
C
C        DWORK := (A*A  + G*Q)(J+1:N,J).
C
         CALL DGEMV( 'NoTranspose', N-J, N, ONE, A(J+1,1), LDA, A(1,J),
     $               1, ZERO, DWORK(J+1), 1 )
         CALL DGEMV( 'Transpose', J, N-J, ONE, QG(1,J+2), LDQG, QG(J,1),
     $               LDQG, ONE, DWORK(J+1), 1 )
         CALL DSYMV( 'Upper', N-J, ONE, QG(J+1,J+2), LDQG, QG(J+1,J), 1,
     $               ONE, DWORK(J+1), 1 )
C
C        Symplectic reflection to zero (H*H)(J+2:N,J).
C
         CALL DLARFG( N-J, DWORK(J+1), DWORK(J+2), 1, TAU )
         DWORK(J+1) = ONE
C
         CALL DLARFX( 'Left', N-J, N, DWORK(J+1), TAU, A(J+1,1), LDA,
     $                DWORK(N+1) )
         CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, A(1,J+1), LDA,
     $                DWORK(N+1) )
C
         CALL DLARFX( 'Left', N-J, J, DWORK(J+1), TAU, QG(J+1,1), LDQG,
     $                DWORK(N+1) )
         CALL DSYMV( 'Lower', N-J, TAU, QG(J+1,J+1), LDQG, DWORK(J+1),
     $               1, ZERO, DWORK(N+J+1), 1 )
         CALL DAXPY( N-J, -TAU*DDOT( N-J, DWORK(N+J+1), 1, DWORK(J+1),
     $               1 )/TWO, DWORK(J+1), 1, DWORK(N+J+1), 1 )
         CALL DSYR2( 'Lower', N-J, -ONE, DWORK(J+1), 1, DWORK(N+J+1), 1,
     $               QG(J+1,J+1), LDQG )
C
         CALL DLARFX( 'Right', J, N-J, DWORK(J+1), TAU, QG(1,J+2), LDQG,
     $                DWORK(N+1) )
         CALL DSYMV( 'Upper', N-J, TAU, QG(J+1,J+2), LDQG, DWORK(J+1),
     $               1, ZERO, DWORK(N+J+1), 1 )
         CALL DAXPY( N-J, -TAU*DDOT( N-J, DWORK(N+J+1), 1, DWORK(J+1),
     $               1 )/TWO, DWORK(J+1), 1, DWORK(N+J+1), 1 )
         CALL DSYR2( 'Upper', N-J, -ONE, DWORK(J+1), 1, DWORK(N+J+1), 1,
     $               QG(J+1,J+2), LDQG )
C
         IF ( FORM ) THEN
C
C           Save reflection.
C
            CALL DCOPY( N-J, DWORK(J+1), 1, U(J+1,N+J), 1 )
            U(J+1,N+J) = TAU
C
         ELSE IF ( ACCUM ) THEN
C
C           Accumulate reflection.
C
            CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, U(1,J+1),
     $                   LDU, DWORK(N+1) )
            CALL DLARFX( 'Right', N, N-J, DWORK(J+1), TAU, U(1,N+J+1),
     $                   LDU, DWORK(N+1) )
         END IF
C
   10 CONTINUE
C
      IF ( FORM ) THEN
         DUMMY(1) = ZERO
C
C        Form S by accumulating transformations.
C
         DO 20 J = N - 1, 1, -1
C
C           Initialize (J+1)st column of S.
C
            CALL DCOPY( N, DUMMY, 0, U(1,J+1), 1 )
            U(J+1,J+1) = ONE
            CALL DCOPY( N, DUMMY, 0, U(1,N+J+1), 1 )
C
C           Second reflection.
C
            TAU = U(J+1,N+J)
            U(J+1,N+J) = ONE
            CALL DLARFX( 'Left', N-J, N-J, U(J+1,N+J), TAU,
     $                   U(J+1,J+1), LDU, DWORK(N+1) )
            CALL DLARFX( 'Left', N-J, N-J, U(J+1,N+J), TAU,
     $                   U(J+1,N+J+1), LDU, DWORK(N+1) )
C
C           Rotation.
C
            CALL DROT( N-J, U(J+1,J+1), LDU, U(J+1,N+J+1), LDU,
     $                 U(J,J), U(J,N+J) )
C
C           First reflection.
C
            TAU = U(J+1,J)
            U(J+1,J) = ONE
            CALL DLARFX( 'Left', N-J, N-J, U(J+1,J), TAU, U(J+1,J+1),
     $                   LDU, DWORK(N+1) )
            CALL DLARFX( 'Left', N-J, N-J, U(J+1,J), TAU,
     $                   U(J+1,N+J+1), LDU, DWORK(N+1) )
   20    CONTINUE
C
C        The first column is the first column of identity.
C
         CALL DCOPY( N, DUMMY, 0, U, 1 )
         U(1,1) = ONE
         CALL DCOPY( N, DUMMY, 0, U(1,N+1), 1 )
      END IF
C
      RETURN
C     *** Last line of MB04ZD ***
      END
