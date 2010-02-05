      SUBROUTINE MB03RD( JOBX, SORT, N, PMAX, A, LDA, X, LDX, NBLCKS,
     $                   BLSIZE, WR, WI, TOL, DWORK, INFO )
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
C     To reduce a matrix A in real Schur form to a block-diagonal form
C     using well-conditioned non-orthogonal similarity transformations.
C     The condition numbers of the transformations used for reduction
C     are roughly bounded by PMAX*PMAX, where PMAX is a given value.
C     The transformations are optionally postmultiplied in a given
C     matrix X. The real Schur form is optionally ordered, so that
C     clustered eigenvalues are grouped in the same block.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBX    CHARACTER*1
C             Specifies whether or not the transformations are
C             accumulated, as follows:
C             = 'N':  The transformations are not accumulated;
C             = 'U':  The transformations are accumulated in X (the
C                     given matrix X is updated).
C
C     SORT    CHARACTER*1
C             Specifies whether or not the diagonal blocks of the real
C             Schur form are reordered, as follows:
C             = 'N':  The diagonal blocks are not reordered;
C             = 'S':  The diagonal blocks are reordered before each
C                     step of reduction, so that clustered eigenvalues
C                     appear in the same block;
C             = 'C':  The diagonal blocks are not reordered, but the
C                     "closest-neighbour" strategy is used instead of
C                     the standard "closest to the mean" strategy
C                     (see METHOD);
C             = 'B':  The diagonal blocks are reordered before each
C                     step of reduction, and the "closest-neighbour"
C                     strategy is used (see METHOD).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and X.  N >= 0.
C
C     PMAX    (input) DOUBLE PRECISION
C             An upper bound for the infinity norm of elementary
C             submatrices of the individual transformations used for
C             reduction (see METHOD).  PMAX >= 1.0D0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A to be block-diagonalized, in real
C             Schur form.
C             On exit, the leading N-by-N part of this array contains
C             the computed block-diagonal matrix, in real Schur
C             canonical form. The non-diagonal blocks are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, if JOBX = 'U', the leading N-by-N part of this
C             array must contain a given matrix X.
C             On exit, if JOBX = 'U', the leading N-by-N part of this
C             array contains the product of the given matrix X and the
C             transformation matrix that reduced A to block-diagonal
C             form. The transformation matrix is itself a product of
C             non-orthogonal similarity transformations having elements
C             with magnitude less than or equal to PMAX.
C             If JOBX = 'N', this array is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of array X.
C             LDX >= 1,        if JOBX = 'N';
C             LDX >= MAX(1,N), if JOBX = 'U'.
C
C     NBLCKS  (output) INTEGER
C             The number of diagonal blocks of the matrix A.
C
C     BLSIZE  (output) INTEGER array, dimension (N)
C             The first NBLCKS elements of this array contain the orders
C             of the resulting diagonal blocks of the matrix A.
C
C     WR,     (output) DOUBLE PRECISION arrays, dimension (N)
C     WI      These arrays contain the real and imaginary parts,
C             respectively, of the eigenvalues of the matrix A.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in the ordering of the diagonal
C             blocks of the real Schur form matrix.
C             If the user sets TOL > 0, then the given value of TOL is
C             used as an absolute tolerance: a block i and a temporarily
C             fixed block 1 (the first block of the current trailing
C             submatrix to be reduced) are considered to belong to the
C             same cluster if their eigenvalues satisfy
C
C               | lambda_1 - lambda_i | <= TOL.
C
C             If the user sets TOL < 0, then the given value of TOL is
C             used as a relative tolerance: a block i and a temporarily
C             fixed block 1 are considered to belong to the same cluster
C             if their eigenvalues satisfy, for j = 1, ..., N,
C
C               | lambda_1 - lambda_i | <= | TOL | * max | lambda_j |.
C
C             If the user sets TOL = 0, then an implicitly computed,
C             default tolerance, defined by TOL = SQRT( SQRT( EPS ) )
C             is used instead, as a relative tolerance, where EPS is
C             the machine precision (see LAPACK Library routine DLAMCH).
C             If SORT = 'N' or 'C', this parameter is not referenced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
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
C     Consider first that SORT = 'N'. Let
C
C            ( A    A   )
C            (  11   12 )
C        A = (          ),
C            ( 0    A   )
C            (       22 )
C
C     be the given matrix in real Schur form, where initially A   is the
C                                                              11
C     first diagonal block of dimension 1-by-1 or 2-by-2. An attempt is
C     made to compute a transformation matrix X of the form
C
C            ( I   P )
C        X = (       )                                               (1)
C            ( 0   I )
C
C     (partitioned as A), so that
C
C                 ( A     0  )
C         -1      (  11      )
C        X  A X = (          ),
C                 ( 0    A   )
C                 (       22 )
C
C     and the elements of P do not exceed the value PMAX in magnitude.
C     An adaptation of the standard method for solving Sylvester
C     equations [1], which controls the magnitude of the individual
C     elements of the computed solution [2], is used to obtain matrix P.
C     When this attempt failed, an 1-by-1 (or 2-by-2) diagonal block of
C     A  , whose eigenvalue(s) is (are) the closest to the mean of those
C      22
C     of A   is selected, and moved by orthogonal similarity
C         11
C     transformations in the leading position of A  ; the moved diagonal
C                                                 22
C     block is then added to the block A  , increasing its order by 1
C                                       11
C     (or 2). Another attempt is made to compute a suitable
C     transformation matrix X with the new definitions of the blocks A
C                                                                     11
C     and A  . After a successful transformation matrix X has been
C          22
C     obtained, it postmultiplies the current transformation matrix
C     (if JOBX = 'U'), and the whole procedure is repeated for the
C     matrix A  .
C             22
C
C     When SORT = 'S', the diagonal blocks of the real Schur form are
C     reordered before each step of the reduction, so that each cluster
C     of eigenvalues, defined as specified in the definition of TOL,
C     appears in adjacent blocks. The blocks for each cluster are merged
C     together, and the procedure described above is applied to the
C     larger blocks. Using the option SORT = 'S' will usually provide
C     better efficiency than the standard option (SORT = 'N'), proposed
C     in [2], because there could be no or few unsuccessful attempts
C     to compute individual transformation matrices X of the form (1).
C     However, the resulting dimensions of the blocks are usually
C     larger; this could make subsequent calculations less efficient.
C
C     When SORT = 'C' or 'B', the procedure is similar to that for
C     SORT = 'N' or 'S', respectively, but the block of A   whose
C                                                        22
C     eigenvalue(s) is (are) the closest to those of A   (not to their
C                                                     11
C     mean) is selected and moved to the leading position of A  . This
C                                                             22
C     is called the "closest-neighbour" strategy.
C
C     REFERENCES
C
C     [1] Bartels, R.H. and Stewart, G.W.  T
C         Solution of the matrix equation A X + XB = C.
C         Comm. A.C.M., 15, pp. 820-826, 1972.
C
C     [2] Bavely, C. and Stewart, G.W.
C         An Algorithm for Computing Reducing Subspaces by Block
C         Diagonalization.
C         SIAM J. Numer. Anal., 16, pp. 359-367, 1979.
C
C     [3] Demmel, J.
C         The Condition Number of Equivalence Transformations that
C         Block Diagonalize Matrix Pencils.
C         SIAM J. Numer. Anal., 20, pp. 599-610, 1983.
C
C     NUMERICAL ASPECTS
C                                       3                     4
C     The algorithm usually requires 0(N ) operations, but 0(N ) are
C     possible in the worst case, when all diagonal blocks in the real
C     Schur form of A are 1-by-1, and the matrix cannot be diagonalized
C     by well-conditioned transformations.
C
C     FURTHER COMMENTS
C
C     The individual non-orthogonal transformation matrices used in the
C     reduction of A to a block-diagonal form have condition numbers
C     of the order PMAX*PMAX. This does not guarantee that their product
C     is well-conditioned enough. The routine can be easily modified to
C     provide estimates for the condition numbers of the clusters of
C     eigenvalues.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998.
C     Partly based on the RASP routine BDIAG by A. Varga, German
C     Aerospace Center, DLR Oberpfaffenhofen.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003.
C
C     KEYWORDS
C
C     Diagonalization, orthogonal transformation, real Schur form,
C     Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBX, SORT
      INTEGER           INFO, LDA, LDX, N, NBLCKS
      DOUBLE PRECISION  PMAX, TOL
C     .. Array Arguments ..
      INTEGER           BLSIZE(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), WI(*), WR(*), X(LDX,*)
C     .. Local Scalars ..
      LOGICAL           LJOBX, LSORN, LSORS, LSORT
      CHARACTER         JOBV
      INTEGER           DA11, DA22, I, IERR, J, K, L, L11, L22, L22M1
      DOUBLE PRECISION  C, CAV, D, EDIF, EMAX, RAV, SAFEMN, SC, THRESH
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLAPY2, DNRM2
      EXTERNAL          DLAMCH, DLAPY2, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLABAD, DLASET, DSCAL, MA02AD, MB03QX,
     $                  MB03RX, MB03RY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO  = 0
      LJOBX = LSAME( JOBX, 'U' )
      LSORN = LSAME( SORT, 'N' )
      LSORS = LSAME( SORT, 'S' )
      LSORT = LSAME( SORT, 'B' ) .OR. LSORS
      IF( .NOT.LJOBX .AND. .NOT.LSAME( JOBX, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSORN .AND. .NOT.LSORT .AND.
     $         .NOT.LSAME( SORT, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( PMAX.LT.ONE ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( ( LDX.LT.1 ) .OR. ( LJOBX .AND. LDX.LT.N ) ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03RD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      NBLCKS = 0
      IF( N.EQ.0 )
     $   RETURN
C
C     Set the "safe" minimum positive number with representable
C     reciprocal, and set JOBV parameter for MB03RX routine.
C
      SAFEMN = DLAMCH( 'Safe minimum' )
      SC = ONE / SAFEMN
      CALL DLABAD( SAFEMN, SC )
      SAFEMN = SAFEMN / DLAMCH( 'Precision' )
      JOBV = JOBX
      IF ( LJOBX )
     $   JOBV = 'V'
C
C     Compute the eigenvalues of A and set the tolerance for reordering
C     the eigenvalues in clusters, if needed.
C
      CALL MB03QX( N, A, LDA, WR, WI, INFO )
C
      IF ( LSORT ) THEN
         THRESH = ABS( TOL )
         IF ( THRESH.EQ.ZERO ) THEN
C
C           Use the default tolerance in ordering the blocks.
C
            THRESH = SQRT( SQRT( DLAMCH( 'Epsilon' ) ) )
         END IF
C
         IF ( TOL.LE.ZERO ) THEN
C
C           Use a relative tolerance. Find max | lambda_j |, j = 1 : N.
C
            EMAX = ZERO
            L = 1
C           WHILE ( L.LE.N ) DO
   10       IF ( L.LE.N ) THEN
               IF ( WI(L).EQ.ZERO ) THEN
                  EMAX = MAX( EMAX, ABS( WR(L) ) )
                  L = L + 1
               ELSE
                  EMAX = MAX( EMAX, DLAPY2( WR(L), WI(L) ) )
                  L = L + 2
               END IF
               GO TO 10
            END IF
C           END WHILE 10
            THRESH = THRESH * EMAX
         END IF
      END IF
C
C     Define the following submatrices of A:
C     A11, the DA11-by-DA11 block in position (L11,L11);
C     A22, the DA22-by-DA22 block in position (L22,L22);
C     A12, the DA11-by-DA22 block in position (L11,L22);
C     A21, the DA22-by-DA11 block in position (L22,L11) (null initially
C                                                        and finally).
C     The following loop uses L11 as loop variable and try to separate a
C     block in position (L11,L11), with possibly clustered eigenvalues,
C     separated by the other eigenvalues (in the block A22).
C
      L11 = 1
C     WHILE ( L11.LE.N ) DO
   20 IF ( L11.LE.N ) THEN
         NBLCKS = NBLCKS + 1
         IF ( WI(L11).EQ.ZERO ) THEN
            DA11 = 1
         ELSE
            DA11 = 2
         END IF
C
         IF ( LSORT ) THEN
C
C           The following loop, using K as loop variable, finds the
C           blocks whose eigenvalues are close to those of A11 and
C           moves these blocks (if any) to the leading position of A22.
C
            L22 = L11 + DA11
            K = L22
C           WHILE ( K.LE.N ) DO
   30       IF ( K.LE.N ) THEN
               EDIF = DLAPY2( WR(L11) - WR(K), WI(L11) - WI(K) )
               IF ( EDIF.LE.THRESH ) THEN
C
C                 An 1x1 or a 2x2 block of A22 has been found so that
C
C                    abs( lambda_1 - lambda_k ) <= THRESH
C
C                 where lambda_1 and lambda_k denote an eigenvalue
C                 of A11 and of that block in A22, respectively.
C                 Try to move that block to the leading position of A22.
C
                  CALL MB03RX( JOBV, N, L22, K, A, LDA, X, LDX, WR, WI,
     $                         DWORK )
C
C                 Extend A11 with the leading block of A22.
C
                  IF ( WI(L22).EQ.ZERO ) THEN
                     DA11 = DA11 + 1
                  ELSE
                     DA11 = DA11 + 2
                  END IF
                  L22 = L11 + DA11
               END IF
               IF ( WI(K).EQ.ZERO ) THEN
                  K = K + 1
               ELSE
                  K = K + 2
               END IF
               GO TO 30
            END IF
C           END WHILE 30
         END IF
C
C        The following loop uses L22 as loop variable and forms a
C        separable DA11-by-DA11 block A11 in position (L11,L11).
C
         L22   = L11 + DA11
         L22M1 = L22 - 1
C        WHILE ( L22.LE.N ) DO
   40    IF ( L22.LE.N ) THEN
            DA22  = N - L22M1
C
C           Try to separate the block A11 of order DA11 by using a
C           well-conditioned similarity transformation.
C
C           First save A12' in the block A21.
C
            CALL MA02AD( 'Full', DA11, DA22, A(L11,L22), LDA,
     $                   A(L22,L11), LDA )
C
C           Solve  -A11*P + P*A22 = A12.
C
            CALL MB03RY( DA11, DA22, PMAX, A(L11,L11), LDA, A(L22,L22),
     $                   LDA, A(L11,L22), LDA, IERR )
C
            IF ( IERR.EQ.1 ) THEN
C
C              The annihilation of A12 failed. Restore A12 and A21.
C
               CALL MA02AD( 'Full', DA22, DA11, A(L22,L11), LDA,
     $                      A(L11,L22), LDA )
               CALL DLASET( 'Full', DA22, DA11, ZERO, ZERO, A(L22,L11),
     $                      LDA )
C
               IF ( LSORN .OR. LSORS ) THEN
C
C                 Extend A11 with an 1x1 or 2x2 block of A22 having the
C                 nearest eigenvalues to the mean of eigenvalues of A11
C                 and resume the loop.
C                 First compute the mean of eigenvalues of A11.
C
                  RAV = ZERO
                  CAV = ZERO
C
                  DO 50 I = L11, L22M1
                     RAV = RAV + WR(I)
                     CAV = CAV + ABS( WI(I) )
   50             CONTINUE
C
                  RAV = RAV/DA11
                  CAV = CAV/DA11
C
C                 Loop to find the eigenvalue of A22 nearest to the
C                 above computed mean.
C
                  D = DLAPY2( RAV-WR(L22), CAV-WI(L22) )
                  K = L22
                  IF ( WI(L22).EQ.ZERO ) THEN
                     L = L22 + 1
                  ELSE
                     L = L22 + 2
                  END IF
C                 WHILE ( L.LE.N ) DO
   60             IF ( L.LE.N ) THEN
                     C = DLAPY2( RAV-WR(L), CAV-WI(L) )
                     IF ( C.LT.D ) THEN
                        D = C
                        K = L
                     END IF
                     IF ( WI(L).EQ.ZERO ) THEN
                        L = L + 1
                     ELSE
                        L = L + 2
                     END IF
                     GO TO 60
                  END IF
C                 END WHILE 60
C
               ELSE
C
C                 Extend A11 with an 1x1 or 2x2 block of A22 having the
C                 nearest eigenvalues to the cluster of eigenvalues of
C                 A11 and resume the loop.
C
C                 Loop to find the eigenvalue of A22 of minimum distance
C                 to the cluster.
C
                  D = SC
                  L = L22
                  K = L22
C                 WHILE ( L.LE.N ) DO
   70             IF ( L.LE.N ) THEN
                     I = L11
C                    WHILE ( I.LE.L22M1 ) DO
   80                IF ( I.LE.L22M1 ) THEN
                        C = DLAPY2( WR(I)-WR(L), WI(I)-WI(L) )
                        IF ( C.LT.D ) THEN
                           D = C
                           K = L
                        END IF
                        IF ( WI(I).EQ.ZERO ) THEN
                           I = I + 1
                        ELSE
                           I = I + 2
                        END IF
                        GO TO 80
                     END IF
C                    END WHILE 80
                     IF ( WI(L).EQ.ZERO ) THEN
                        L = L + 1
                     ELSE
                        L = L + 2
                     END IF
                     GO TO 70
                  END IF
C                 END WHILE 70
               END IF
C
C              Try to move block found to the leading position of A22.
C
               CALL MB03RX( JOBV, N, L22, K, A, LDA, X, LDX, WR, WI,
     $                      DWORK )
C
C              Extend A11 with the leading block of A22.
C
               IF ( WI(L22).EQ.ZERO ) THEN
                  DA11 = DA11 + 1
               ELSE
                  DA11 = DA11 + 2
               END IF
               L22   = L11 + DA11
               L22M1 = L22 - 1
               GO TO 40
            END IF
         END IF
C        END WHILE 40
C
         IF ( LJOBX ) THEN
C
C           Accumulate the transformation in X.
C           Only columns L22, ..., N are modified.
C
            IF ( L22.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, DA22,
     $                     DA11, ONE, X(1,L11), LDX, A(L11,L22), LDA,
     $                     ONE, X(1,L22), LDX )
C
C           Scale to unity the (non-zero) columns of X which will be
C           no more modified and transform A11 accordingly.
C
            DO 90 J = L11, L22M1
               SC = DNRM2( N, X(1,J), 1 )
               IF ( SC.GT.SAFEMN ) THEN
                  CALL DSCAL( DA11, SC, A(J,L11), LDA )
                  SC = ONE/SC
                  CALL DSCAL( N, SC, X(1,J), 1 )
                  CALL DSCAL( DA11, SC, A(L11,J), 1 )
               END IF
   90       CONTINUE
C
         END IF
         IF ( L22.LE.N ) THEN
C
C           Set A12 and A21 to zero.
C
            CALL DLASET( 'Full', DA11, DA22, ZERO, ZERO, A(L11,L22),
     $                   LDA )
            CALL DLASET( 'Full', DA22, DA11, ZERO, ZERO, A(L22,L11),
     $                   LDA )
         END IF
C
C        Store the orders of the diagonal blocks in BLSIZE.
C
         BLSIZE(NBLCKS) = DA11
         L11 = L22
         GO TO 20
      END IF
C     END WHILE 20
C
      RETURN
C *** Last line of MB03RD ***
      END
