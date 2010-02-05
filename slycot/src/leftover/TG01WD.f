      SUBROUTINE TG01WD( N, M, P, A, LDA, E, LDE, B, LDB, C, LDC,
     $                   Q, LDQ, Z, LDZ, ALPHAR, ALPHAI, BETA, DWORK,
     $                   LDWORK, INFO )
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
C     To reduce the pair (A,E) to a real generalized Schur form
C     by using an orthogonal equivalence transformation
C     (A,E) <-- (Q'*A*Z,Q'*E*Z) and to apply the transformation
C     to the matrices B and C: B <-- Q'*B and C <-- C*Z.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e., the order of the matrices A and E.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs, or of columns of B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs, or of rows of C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Q' * A * Z in an upper quasi-triangular form.
C             The elements below the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original descriptor matrix E.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Q' * E * Z in an upper triangular form.
C             The elements below the diagonal are set to zero.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed input matrix Q' * B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed output matrix C * Z.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
C             The leading N-by-N part of this array contains the left
C             orthogonal transformation matrix used to reduce (A,E) to
C             the real generalized Schur form.
C             The columns of Q are the left generalized Schur vectors
C             of the pair (A,E).
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= max(1,N).
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             The leading N-by-N part of this array contains the right
C             orthogonal transformation matrix used to reduce (A,E) to
C             the real generalized Schur form.
C             The columns of Z are the right generalized Schur vectors
C             of the pair (A,E).
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= max(1,N).
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C             On exit, if INFO = 0, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j),
C             j=1,...,N, will be the generalized eigenvalues.
C             ALPHAR(j) + ALPHAI(j)*i, and BETA(j), j=1,...,N, are the
C             diagonals of the complex Schur form that would result if
C             the 2-by-2 diagonal blocks of the real Schur form of
C             (A,E) were further reduced to triangular form using
C             2-by-2 complex unitary transformations.
C             If ALPHAI(j) is zero, then the j-th eigenvalue is real;
C             if positive, then the j-th and (j+1)-st eigenvalues are a
C             complex conjugate pair, with ALPHAI(j+1) negative.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.  LDWORK >= 8*N+16.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, the QZ algorithm failed to compute
C                   the generalized real Schur form; elements i+1:N of
C                   ALPHAR, ALPHAI, and BETA should be correct.
C
C     METHOD
C
C     The pair (A,E) is reduced to a real generalized Schur form using
C     an orthogonal equivalence transformation (A,E) <-- (Q'*A*Z,Q'*E*Z)
C     and the transformation is applied to the matrices B and C:
C     B <-- Q'*B and C <-- C*Z.
C
C     NUMERICAL ASPECTS
C                                     3
C     The algorithm requires about 25N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     KEYWORDS
C
C     Orthogonal transformation, generalized real Schur form, similarity
C     transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDE, LDQ, LDWORK, LDZ,
     $                  M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ALPHAI(*), ALPHAR(*), B(LDB,*),
     $                  BETA(*),  C(LDC,*),  DWORK(*),  E(LDE,*),
     $                  Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           BLAS3, BLOCK
      INTEGER           BL, CHUNK, I, J, MAXWRK, SDIM
C     .. Local Arrays ..
      LOGICAL           BWORK(1)
C     .. External Functions ..
      LOGICAL           LSAME, DELCTG
      EXTERNAL          LSAME, DELCTG
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DGGES, DLACPY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
      INFO = 0
C
C     Check the scalar input parameters.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDWORK.LT.8*N+16 ) THEN
         INFO = -20
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TG01WD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Reduce (A,E) to real generalized Schur form using an orthogonal
C     equivalence transformation (A,E) <-- (Q'*A*Z,Q'*E*Z), accumulate
C     the transformations in Q and Z, and compute the generalized
C     eigenvalues of the pair (A,E) in (ALPHAR, ALPHAI, BETA).
C
C     Workspace:  need   8*N+16;
C                 prefer larger.
C
      CALL DGGES( 'Vectors', 'Vectors', 'Not ordered', DELCTG, N,
     $            A, LDA, E, LDE, SDIM, ALPHAR, ALPHAI, BETA, Q, LDQ,
     $            Z, LDZ, DWORK, LDWORK, BWORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      MAXWRK = INT( DWORK(1) )
C
C     Apply the transformation: B <-- Q'*B. Use BLAS 3, if enough space.
C
      CHUNK = LDWORK / N
      BLOCK = M.GT.1
      BLAS3 = CHUNK.GE.M .AND. BLOCK
C
      IF( BLAS3 ) THEN
C
C        Enough workspace for a fast BLAS 3 algorithm.
C
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
         CALL DGEMM( 'Transpose', 'No transpose', N, M, N, ONE, Q, LDQ,
     $               DWORK, N, ZERO, B, LDB )
C
      ELSE IF ( BLOCK ) THEN
C
C        Use as many columns of B as possible.
C
         DO 10 J = 1, M, CHUNK
            BL = MIN( M-J+1, CHUNK )
            CALL DLACPY( 'Full', N, BL, B(1,J), LDB, DWORK, N )
            CALL DGEMM( 'Transpose', 'NoTranspose', N, BL, N, ONE, Q,
     $                  LDQ, DWORK, N, ZERO, B(1,J), LDB )
   10    CONTINUE
C
      ELSE
C
C        Use a BLAS 2 algorithm. Here, M <= 1.
C
         IF ( M.GT.0 ) THEN
            CALL DCOPY( N, B, 1, DWORK, 1 )
            CALL DGEMV( 'Transpose', N, N, ONE, Q, LDQ, DWORK, 1, ZERO,
     $                  B, 1 )
         END IF
      END IF
      MAXWRK = MAX( MAXWRK, N*M )
C
C     Apply the transformation: C <-- C*Z.  Use BLAS 3, if enough space.
C
      BLOCK = P.GT.1
      BLAS3 = CHUNK.GE.P .AND. BLOCK
C
      IF ( BLAS3 ) THEN
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK, P )
         CALL DGEMM( 'No transpose', 'No transpose', P, N, N, ONE,
     $               DWORK, P, Z, LDZ, ZERO, C, LDC )
C
      ELSE IF ( BLOCK ) THEN
C
C        Use as many rows of C as possible.
C
         DO 20 I = 1, P, CHUNK
            BL = MIN( P-I+1, CHUNK )
            CALL DLACPY( 'Full', BL, N, C(I,1), LDC, DWORK, BL )
            CALL DGEMM( 'NoTranspose', 'NoTranspose', BL, N, N, ONE,
     $                  DWORK, BL, Z, LDZ, ZERO, C(I,1), LDC )
   20    CONTINUE
C
      ELSE
C
C        Use a BLAS 2 algorithm. Here, P <= 1.
C
         IF ( P.GT.0 ) THEN
            CALL DCOPY( N, C, LDC, DWORK, 1 )
            CALL DGEMV( 'Transpose', N, N, ONE, Z, LDZ, DWORK, 1, ZERO,
     $                  C, LDC )
         END IF
C
      END IF
      MAXWRK = MAX( MAXWRK, P*N )
C
      DWORK(1) = DBLE( MAXWRK )
C
      RETURN
C *** Last line of TG01WD ***
      END
