      SUBROUTINE MB04DY( JOBSCL, N, A, LDA, QG, LDQG, D, DWORK, INFO )
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
C     To perform a symplectic scaling on the Hamiltonian matrix
C
C              ( A    G  )
C          H = (       T ),                                          (1)
C              ( Q   -A  )
C
C     i.e., perform either the symplectic scaling transformation
C
C                                   -1
C                 ( A'   G'  )   ( D   0 ) ( A   G  ) ( D  0   )
C          H' <-- (        T ) = (       ) (      T ) (     -1 ),    (2)
C                 ( Q'  -A'  )   ( 0   D ) ( Q  -A  ) ( 0  D   )
C
C     where D is a diagonal scaling matrix, or the symplectic norm
C     scaling transformation
C
C                  ( A''   G''  )    1  (   A   G/tau )
C          H'' <-- (          T ) = --- (           T ),             (3)
C                  ( Q''  -A''  )   tau ( tau Q   -A  )
C
C     where tau is a real scalar.  Note that if tau is not equal to 1,
C     then (3) is NOT a similarity transformation.  The eigenvalues
C     of H are then tau times the eigenvalues of H''.
C
C     For symplectic scaling (2), D is chosen to give the rows and
C     columns of A' approximately equal 1-norms and to give Q' and G'
C     approximately equal norms.  (See METHOD below for details.) For
C     norm scaling, tau = MAX(1, ||A||, ||G||, ||Q||) where ||.||
C     denotes the 1-norm (column sum norm).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBSCL  CHARACTER*1
C             Indicates which scaling strategy is used, as follows:
C             = 'S'       :  do the symplectic scaling (2);
C             = '1' or 'O':  do the 1-norm scaling (3);
C             = 'N'       :  do nothing; set INFO and return.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, G, and Q.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, if JOBSCL <> 'N', the leading N-by-N part of
C             this array must contain the upper left block A of the
C             Hamiltonian matrix H in (1).
C             On output, if JOBSCL <> 'N', the leading N-by-N part of
C             this array contains the leading N-by-N part of the scaled
C             Hamiltonian matrix H' in (2) or H'' in (3), depending on
C             the setting of JOBSCL.
C             If JOBSCL = 'N', this array is not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,N), if JOBSCL <> 'N';
C             LDA >= 1,        if JOBSCL =  'N'.
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C             (LDQG,N+1)
C             On input, if JOBSCL <> 'N', the leading N-by-N lower
C             triangular part of this array must contain the lower
C             triangle of the lower left symmetric block Q of the
C             Hamiltonian matrix H in (1), and the N-by-N upper
C             triangular part of the submatrix in the columns 2 to N+1
C             of this array must contain the upper triangle of the upper
C             right symmetric block G of H in (1).
C             So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j)
C             and G(i,j) = G(j,i) is stored in QG(j,i+1).
C             On output, if JOBSCL <> 'N', the leading N-by-N lower
C             triangular part of this array contains the lower triangle
C             of the lower left symmetric block Q' or Q'', and the
C             N-by-N upper triangular part of the submatrix in the
C             columns 2 to N+1 of this array contains the upper triangle
C             of the upper right symmetric block G' or G'' of the scaled
C             Hamiltonian matrix H' in (2) or H'' in (3), depending on
C             the setting of JOBSCL.
C             If JOBSCL = 'N', this array is not referenced.
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.
C             LDQG >= MAX(1,N), if JOBSCL <> 'N';
C             LDQG >= 1,        if JOBSCL =  'N'.
C
C     D       (output) DOUBLE PRECISION array, dimension (nd)
C             If JOBSCL = 'S', then nd = N and D contains the diagonal
C             elements of the diagonal scaling matrix in (2).
C             If JOBSCL = '1' or 'O', then nd = 1 and D(1) is set to tau
C             from (3). In this case, no other elements of D are
C             referenced.
C             If JOBSCL = 'N', this array is not referenced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C             If JOBSCL = 'N', this array is not referenced.
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
C     1. Symplectic scaling (JOBSCL = 'S'):
C
C     First, LAPACK subroutine DGEBAL is used to equilibrate the 1-norms
C     of the rows and columns of A using a diagonal scaling matrix D_A.
C     Then, H is similarily transformed by the symplectic diagonal
C     matrix D1 = diag(D_A,D_A**(-1)).  Next, the off-diagonal blocks of
C     the resulting Hamiltonian matrix are equilibrated in the 1-norm
C     using the symplectic diagonal matrix D2 of the form
C
C                 ( I/rho    0   )
C            D2 = (              )
C                 (   0    rho*I )
C
C     where rho is a real scalar. Thus, in (2), D = D1*D2.
C
C     2. Norm scaling (JOBSCL = '1' or 'O'):
C
C     The norm of the matrices A and G of (1) is reduced by setting
C     A := A/tau  and  G := G/(tau**2) where tau is the power of the
C     base of the arithmetic closest to MAX(1, ||A||, ||G||, ||Q||) and
C     ||.|| denotes the 1-norm.
C
C     REFERENCES
C
C     [1] Benner, P., Byers, R., and Barth, E.
C         Fortran 77 Subroutines for Computing the Eigenvalues of
C         Hamiltonian Matrices. I: The Square-Reduced Method.
C         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000.
C
C     NUMERICAL ASPECTS
C
C     For symplectic scaling, the complexity of the used algorithms is
C     hard to estimate and depends upon how well the rows and columns of
C     A in (1) are equilibrated.  In one sweep, each row/column of A is
C     scaled once, i.e., the cost of one sweep is N**2 multiplications.
C     Usually, 3-6 sweeps are enough to equilibrate the norms of the
C     rows and columns of a matrix.  Roundoff errors are possible as
C     LAPACK routine DGEBAL does NOT use powers of the machine base for
C     scaling. The second stage (equilibrating ||G|| and ||Q||) requires
C     N**2 multiplications.
C     For norm scaling, 3*N**2 + O(N) multiplications are required and
C     NO rounding errors occur as all multiplications are performed with
C     powers of the machine base.
C
C     CONTRIBUTOR
C
C     P. Benner, Universitaet Bremen, Germany, and
C     R. Byers, University of Kansas, Lawrence, USA.
C     Aug. 1998, routine DHABL.
C     V. Sima, Research Institute for Informatics, Bucharest, Romania,
C     Oct. 1998, SLICOT Library version.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2009.
C
C     KEYWORDS
C
C     Balancing, Hamiltonian matrix, norms, symplectic similarity
C     transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDQG, N
      CHARACTER         JOBSCL
C    ..
C    .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), D(*), DWORK(*), QG(LDQG,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  ANRM, BASE, EPS, GNRM, OFL, QNRM,
     $                  RHO, SFMAX, SFMIN, TAU, UFL, Y
      INTEGER           I, IERR, IHI, ILO, J
      LOGICAL           NONE, NORM, SYMP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANSY
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DLANGE, DLANSY, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL          DGEBAL, DLABAD, DLASCL, DRSCL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
      INFO = 0
      SYMP = LSAME( JOBSCL, 'S' )
      NORM = LSAME( JOBSCL, '1' ) .OR. LSAME( JOBSCL, 'O' )
      NONE = LSAME( JOBSCL, 'N' )
C
      IF( .NOT.SYMP .AND. .NOT.NORM .AND. .NOT.NONE ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF(  LDA.LT.1 .OR. ( .NOT.NONE .AND.  LDA.LT.N ) ) THEN
         INFO = -4
      ELSE IF( LDQG.LT.1 .OR. ( .NOT.NONE .AND. LDQG.LT.N ) ) THEN
         INFO = -6
      END IF
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04DY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. NONE )
     $   RETURN
C
C     Set some machine dependant constants.
C
      BASE = DLAMCH( 'Base' )
      EPS  = DLAMCH( 'Precision' )
      UFL  = DLAMCH( 'Safe minimum' )
      OFL  = ONE/UFL
      CALL DLABAD( UFL, OFL )
      SFMAX = ( EPS/BASE )/UFL
      SFMIN = ONE/SFMAX
C
      IF ( NORM ) THEN
C
C        Compute norms.
C
         ANRM = DLANGE( '1-norm', N, N, A, LDA, DWORK )
         GNRM = DLANSY( '1-norm', 'Upper', N, QG(1,2), LDQG, DWORK )
         QNRM = DLANSY( '1-norm', 'Lower', N, QG, LDQG, DWORK )
         Y    = MAX( ONE, ANRM, GNRM, QNRM )
         TAU  = ONE
C
C        WHILE ( TAU < Y ) DO
   10    CONTINUE
         IF ( ( TAU.LT.Y ) .AND. ( TAU.LT.SQRT( SFMAX ) ) ) THEN
            TAU = TAU*BASE
            GO TO 10
         END IF
C        END WHILE 10
         IF ( TAU.GT.ONE ) THEN
            IF ( ABS( TAU/BASE - Y ).LT.ABS( TAU - Y ) )
     $         TAU = TAU/BASE
            CALL DLASCL( 'General', 0, 0, TAU, ONE, N, N, A, LDA, IERR )
            CALL DLASCL( 'Upper', 0, 0, TAU, ONE, N, N, QG(1,2), LDQG,
     $                   IERR )
            CALL DLASCL( 'Upper', 0, 0, TAU, ONE, N, N, QG(1,2), LDQG,
     $                   IERR )
         END IF
C
         D(1) = TAU
C
      ELSE
         CALL DGEBAL( 'Scale', N, A, LDA, ILO, IHI, D, IERR )
C
         DO 30 J = 1, N
C
            DO 20 I = J, N
               QG(I,J) = QG(I,J)*D(J)*D(I)
   20       CONTINUE
C
   30    CONTINUE
C
         DO 50 J = 2, N + 1
C
            DO 40 I = 1, J - 1
               QG(I,J) = QG(I,J)/D(J-1)/D(I)
   40       CONTINUE
C
   50    CONTINUE
C
         GNRM = DLANSY( '1-norm', 'Upper', N, QG(1,2), LDQG, DWORK )
         QNRM = DLANSY( '1-norm', 'Lower', N, QG, LDQG, DWORK )
         IF ( GNRM.EQ.ZERO ) THEN
            IF ( QNRM.EQ.ZERO ) THEN
               RHO = ONE
            ELSE
               RHO = SFMAX
            END IF
         ELSE IF ( QNRM.EQ.ZERO ) THEN
            RHO = SFMIN
         ELSE
            RHO = SQRT( QNRM )/SQRT( GNRM )
         END IF
C
         CALL DLASCL( 'Lower', 0, 0, RHO, ONE, N, N, QG, LDQG, IERR )
         CALL DLASCL( 'Upper', 0, 0, ONE, RHO, N, N, QG(1,2), LDQG,
     $                IERR )
         CALL DRSCL( N, SQRT( RHO ), D, 1 )
      END IF
C
      RETURN
C     *** Last line of MB04DY ***
      END
