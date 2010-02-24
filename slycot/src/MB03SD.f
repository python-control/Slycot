      SUBROUTINE MB03SD( JOBSCL, N, A, LDA, QG, LDQG, WR, WI, DWORK,
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
C     To compute the eigenvalues of an N-by-N square-reduced Hamiltonian
C     matrix
C
C                ( A'   G'  )
C         H'  =  (        T ).                                       (1)
C                ( Q'  -A'  )
C
C     Here, A' is an N-by-N matrix, and G' and Q' are symmetric N-by-N
C     matrices.  It is assumed without a check that H' is square-
C     reduced, i.e., that
C
C           2    ( A''   G'' )
C         H'  =  (         T )    with A'' upper Hessenberg.         (2)
C                ( 0    A''  )
C
C                            T                2
C     (Equivalently, Q'A'- A' Q' = 0, A'' = A' + G'Q', and for i > j+1,
C      A''(i,j) = 0.)  Ordinarily, H' is the output from SLICOT Library
C     routine MB04ZD. The eigenvalues of H' are computed as the square
C     roots of the eigenvalues of A''.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBSCL  CHARACTER*1
C             Specifies whether or not balancing operations should
C             be performed by the LAPACK subroutine DGEBAL on the
C             Hessenberg matrix A'' in (2), as follows:
C             = 'N':  do not use balancing;
C             = 'S':  do scaling in order to equilibrate the rows
C                     and columns of A''.
C             See LAPACK subroutine DGEBAL and Section METHOD below.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             upper left block A' of the square-reduced Hamiltonian
C             matrix H' in (1), as produced by SLICOT Library routine
C             MB04ZD.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1)
C             The leading N-by-N lower triangular part of this array
C             must contain the lower triangle of the lower left
C             symmetric block Q' of the square-reduced Hamiltonian
C             matrix H' in (1), and the N-by-N upper triangular part of
C             the submatrix in the columns 2 to N+1 of this array must
C             contain the upper triangle of the upper right symmetric
C             block G' of the square-reduced Hamiltonian matrix H'
C             in (1), as produced by SLICOT Library routine MB04ZD.
C             So, if i >= j, then Q'(i,j) is stored in QG(i,j) and
C             G'(i,j) is stored in QG(j,i+1).
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= MAX(1,N).
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             The arrays WR and WI contain the real and imaginary parts,
C             respectively, of the N eigenvalues of H' with non-negative
C             real part.  The remaining N eigenvalues are the negatives
C             of these eigenvalues.
C             Eigenvalues are stored in WR and WI in decreasing order of
C             magnitude of the real parts, i.e., WR(I) >= WR(I+1).
C             (In particular, an eigenvalue closest to the imaginary
C              axis is WR(N)+WI(N)i.)
C             In addition, eigenvalues with zero real part are sorted in
C             decreasing order of magnitude of imaginary parts.  Note
C             that non-real eigenvalues with non-zero real part appear
C             in complex conjugate pairs, but eigenvalues with zero real
C             part do not, in general, appear in complex conjugate
C             pairs.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= MAX(1,N*(N+1)).
C             For good performance, LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, then the i-th argument had an illegal
C                   value;
C             > 0:  if INFO =  i, i <= N, then LAPACK subroutine DHSEQR
C                   failed to converge while computing the i-th
C                   eigenvalue.
C
C     METHOD
C
C     The routine forms the upper Hessenberg matrix A'' in (2) and calls
C     LAPACK subroutines to calculate its eigenvalues.  The eigenvalues
C     of H' are the square roots of the eigenvalues of A''.
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
C     The algorithm requires (32/3)*N**3 + O(N**2) floating point
C     operations.
C     Eigenvalues computed by this subroutine are exact eigenvalues
C     of a perturbed Hamiltonian matrix  H' + E  where
C
C                 || E || <= c sqrt(eps) || H' ||,
C
C     c is a modest constant depending on the dimension N and eps is the
C     machine precision. Moreover, if the norm of H' and an eigenvalue
C     are of roughly the same magnitude, the computed eigenvalue is
C     essentially as accurate as the computed eigenvalue obtained by
C     traditional methods. See [1] or [2].
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C     Aug. 1998, routine DHAEVS.
C     V. Sima, Research Institute for Informatics, Bucharest, Romania,
C     Oct. 1998, SLICOT Library version.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2002,
C     May 2009.
C
C     KEYWORDS
C
C     Eigenvalues, (square-reduced) Hamiltonian matrix, symplectic
C     similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDQG, LDWORK, N
      CHARACTER         JOBSCL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), QG(LDQG,*), WI(*), WR(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  SWAP, X, Y
      INTEGER           BL, CHUNK, I, IGNORE, IHI, ILO, J, JW, JWORK, M,
     $                  N2
      LOGICAL           BLAS3, BLOCK, SCALE, SORTED
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     ..
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEBAL, DGEMM, DHSEQR, DLACPY, DLASET,
     $                  DSYMM, DSYMV, MA01AD, MA02ED, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
      INFO  = 0
      N2    = N*N
      SCALE = LSAME( JOBSCL, 'S' )
      IF ( .NOT. ( SCALE .OR. LSAME( JOBSCL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQG.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDWORK.LT.MAX( 1, N2 + N ) ) THEN
         INFO = -10
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03SD', -INFO )
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
      CHUNK = ( LDWORK - N2 ) / N
      BLOCK = MIN( CHUNK, N ).GT.1
      BLAS3 = CHUNK.GE.N
C
      IF ( BLAS3 ) THEN
         JWORK = N2 + 1
      ELSE
         JWORK = 1
      END IF
C                             2
C     Form the matrix A'' = A'  + G'Q'.
C
      CALL DLACPY( 'Lower', N, N, QG, LDQG, DWORK(JWORK), N )
      CALL MA02ED( 'Lower', N, DWORK(JWORK), N )
C
      IF ( BLAS3 ) THEN
C
C        Use BLAS 3 calculation.
C
         CALL DSYMM( 'Left', 'Upper', N, N, ONE, QG(1, 2), LDQG,
     $               DWORK(JWORK), N, ZERO, DWORK, N )
C
      ELSE IF ( BLOCK ) THEN
         JW = N2 + 1
C
C        Use BLAS 3 for as many columns of Q' as possible.
C
         DO 10 J = 1, N, CHUNK
            BL = MIN( N-J+1, CHUNK )
            CALL DSYMM(  'Left', 'Upper', N, BL, ONE, QG(1, 2), LDQG,
     $                   DWORK(1+N*(J-1)), N, ZERO, DWORK(JW), N )
            CALL DLACPY( 'Full', N, BL, DWORK(JW), N, DWORK(1+N*(J-1)),
     $                   N )
   10    CONTINUE
C
      ELSE
C
C        Use BLAS 2 calculation.
C
         DO 20 J = 1, N
            CALL DSYMV( 'Upper', N, ONE, QG(1, 2), LDQG,
     $                  DWORK(1+N*(J-1)), 1, ZERO, WR, 1 )
            CALL DCOPY( N, WR, 1, DWORK(1+N*(J-1)), 1 )
   20    CONTINUE
C
      END IF
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, N, ONE, A, LDA, A,
     $            LDA, ONE, DWORK, N )
      IF ( SCALE .AND. N.GT.2 )
     $   CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, DWORK(3), N )
C                               2
C     Find the eigenvalues of A' + G'Q'.
C
      CALL DGEBAL( JOBSCL, N, DWORK, N, ILO, IHI, DWORK(1+N2), IGNORE )
      CALL DHSEQR( 'Eigenvalues', 'NoSchurVectors', N, ILO, IHI, DWORK,
     $             N, WR, WI, DUMMY, 1, DWORK(1+N2), N, INFO )
      IF ( INFO.EQ.0 ) THEN
C
C        Eigenvalues of H' are the square roots of those computed above.
C
         DO 30 I = 1, N
            X = WR(I)
            Y = WI(I)
            CALL MA01AD( X, Y, WR(I), WI(I) )
   30    CONTINUE
C
C        Sort eigenvalues into decreasing order by real part and, for
C        eigenvalues with zero real part only, decreasing order of
C        imaginary part.  (This simple bubble sort preserves the
C        relative order of eigenvalues with equal but nonzero real part.
C        This ensures that complex conjugate pairs remain
C        together.)
C
         SORTED = .FALSE.
C
         DO 50 M = N, 1, -1
            IF ( SORTED ) GO TO 60
            SORTED = .TRUE.
C
            DO 40 I = 1, M - 1
               IF ( ( ( WR(I).LT.WR(I+1) ) .OR.
     $              ( ( WR(I).EQ.ZERO ) .AND. ( WR(I+1).EQ.ZERO ) .AND.
     $                ( WI(I).LT.WI(I+1) ) ) ) ) THEN
                  SWAP    = WR(I)
                  WR(I)   = WR(I+1)
                  WR(I+1) = SWAP
                  SWAP    = WI(I)
                  WI(I)   = WI(I+1)
                  WI(I+1) = SWAP
C
                  SORTED = .FALSE.
C
               END IF
   40       CONTINUE
C
   50    CONTINUE
C
   60    CONTINUE
C
      END IF
C
      DWORK(1) = 2*N2
      RETURN
C     *** Last line of MB03SD ***
      END
