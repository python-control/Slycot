      SUBROUTINE MB04VD( MODE, JOBQ, JOBZ, M, N, RANKE, A, LDA, E, LDE,
     $                   Q, LDQ, Z, LDZ, ISTAIR, NBLCKS, NBLCKI, IMUK,
     $                   INUK, IMUK0, MNEI, TOL, IWORK, INFO )
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
C     To compute orthogonal transformations Q and Z such that the
C     transformed pencil Q'(sE-A)Z is in upper block triangular form,
C     where E is an M-by-N matrix in column echelon form (see SLICOT
C     Library routine MB04UD) and A is an M-by-N matrix.
C
C     If MODE = 'B', then the matrices A and E are transformed into the
C     following generalized Schur form by unitary transformations Q1
C     and Z1 :
C
C                      | sE(eps,inf)-A(eps,inf) |      X     |
C        Q1'(sE-A)Z1 = |------------------------|------------|.   (1)
C                      |            O           | sE(r)-A(r) |
C
C     The pencil sE(eps,inf)-A(eps,inf) is in staircase form, and it
C     contains all Kronecker column indices and infinite elementary
C     divisors of the pencil sE-A. The pencil sE(r)-A(r) contains all
C     Kronecker row indices and elementary divisors of sE-A.
C     Note: X is a pencil.
C
C     If MODE = 'T', then the submatrices having full row and column
C     rank in the pencil sE(eps,inf)-A(eps,inf) in (1) are
C     triangularized by applying unitary transformations Q2 and Z2 to
C     Q1'*(sE-A)*Z1.
C
C     If MODE = 'S', then the pencil sE(eps,inf)-A(eps,inf) in (1) is
C     separated into sE(eps)-A(eps) and sE(inf)-A(inf) by applying
C     unitary transformations Q3 and Z3 to Q2'*Q1'*(sE-A)*Z1*Z2.
C
C     This gives
C
C                | sE(eps)-A(eps) |        X       |      X     |
C                |----------------|----------------|------------|
C                |        O       | sE(inf)-A(inf) |      X     |
C     Q'(sE-A)Z =|=================================|============| (2)
C                |                                 |            |
C                |                O                | sE(r)-A(r) |
C
C     where Q = Q1*Q2*Q3 and Z = Z1*Z2*Z3.
C     Note: the pencil sE(r)-A(r) is not reduced further.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     MODE    CHARACTER*1
C             Specifies the desired structure of the transformed
C             pencil Q'(sE-A)Z to be computed as follows:
C             = 'B':  Basic reduction given by (1);
C             = 'T':  Further reduction of (1) to triangular form;
C             = 'S':  Further separation of sE(eps,inf)-A(eps,inf)
C                     in (1) into the two pencils in (2).
C
C     JOBQ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Q the orthogonal row transformations, as follows:
C             = 'N':  Do not form Q;
C             = 'I':  Q is initialized to the unit matrix and the
C                     orthogonal transformation matrix Q is returned;
C             = 'U':  The given matrix Q is updated by the orthogonal
C                     row transformations used in the reduction.
C
C     JOBZ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the orthogonal column transformations, as
C             follows:
C             = 'N':  Do not form Z;
C             = 'I':  Z is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z is returned;
C             = 'U':  The given matrix Z is updated by the orthogonal
C                     transformations used in the reduction.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in the matrices A, E and the order of
C             the matrix Q.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns in the matrices A, E and the order
C             of the matrix Z.  N >= 0.
C
C     RANKE   (input) INTEGER
C             The rank of the matrix E in column echelon form.
C             RANKE >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix to be row compressed.
C             On exit, the leading M-by-N part of this array contains
C             the matrix that has been row compressed while keeping
C             matrix E in column echelon form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix in column echelon form to be
C             transformed equivalent to matrix A.
C             On exit, the leading M-by-N part of this array contains
C             the matrix that has been transformed equivalent to matrix
C             A.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,M).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*)
C             On entry, if JOBQ = 'U', then the leading M-by-M part of
C             this array must contain a given matrix Q (e.g. from a
C             previous call to another SLICOT routine), and on exit, the
C             leading M-by-M part of this array contains the product of
C             the input matrix Q and the row transformation matrix used
C             to transform the rows of matrices A and E.
C             On exit, if JOBQ = 'I', then the leading M-by-M part of
C             this array contains the matrix of accumulated orthogonal
C             row transformations performed.
C             If JOBQ = 'N', the array Q is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDQ = 1 and
C             declare this array to be Q(1,1) in the calling program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q. If JOBQ = 'U' or
C             JOBQ = 'I', LDQ >= MAX(1,M); if JOBQ = 'N', LDQ >= 1.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*)
C             On entry, if JOBZ = 'U', then the leading N-by-N part of
C             this array must contain a given matrix Z (e.g. from a
C             previous call to another SLICOT routine), and on exit, the
C             leading N-by-N part of this array contains the product of
C             the input matrix Z and the column transformation matrix
C             used to transform the columns of matrices A and E.
C             On exit, if JOBZ = 'I', then the leading N-by-N part of
C             this array contains the matrix of accumulated orthogonal
C             column transformations performed.
C             If JOBZ = 'N', the array Z is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDZ = 1 and
C             declare this array to be Z(1,1) in the calling program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If JOBZ = 'U' or
C             JOBZ = 'I', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1.
C
C     ISTAIR  (input/output) INTEGER array, dimension (M)
C             On entry, this array must contain information on the
C             column echelon form of the unitary transformed matrix E.
C             Specifically, ISTAIR(i) must be set to +j if the first
C             non-zero element E(i,j) is a corner point and -j
C             otherwise, for i = 1,2,...,M.
C             On exit, this array contains no useful information.
C
C     NBLCKS  (output) INTEGER
C             The number of submatrices having full row rank greater
C             than or equal to 0 detected in matrix A in the pencil
C             sE(x)-A(x),
C                where  x = eps,inf  if MODE = 'B' or 'T',
C                or     x = eps      if MODE = 'S'.
C
C     NBLCKI  (output) INTEGER
C             If MODE = 'S', the number of diagonal submatrices in the
C             pencil sE(inf)-A(inf). If MODE = 'B' or 'T' then
C             NBLCKI = 0.
C
C     IMUK    (output) INTEGER array, dimension (MAX(N,M+1))
C             The leading NBLCKS elements of this array contain the
C             column dimensions mu(1),...,mu(NBLCKS) of the submatrices
C             having full column rank in the pencil sE(x)-A(x),
C                where  x = eps,inf  if MODE = 'B' or 'T',
C                or     x = eps      if MODE = 'S'.
C
C     INUK    (output) INTEGER array, dimension (MAX(N,M+1))
C             The leading NBLCKS elements of this array contain the
C             row dimensions nu(1),...,nu(NBLCKS) of the submatrices
C             having full row rank in the pencil sE(x)-A(x),
C                where  x = eps,inf  if MODE = 'B' or 'T',
C                or     x = eps      if MODE = 'S'.
C
C     IMUK0   (output) INTEGER array, dimension (limuk0),
C             where limuk0 = N if MODE = 'S' and 1, otherwise.
C             If MODE = 'S', then the leading NBLCKI elements of this
C             array contain the dimensions mu0(1),...,mu0(NBLCKI)
C             of the square diagonal submatrices in the pencil
C             sE(inf)-A(inf).
C             Otherwise, IMUK0 is not referenced and can be supplied
C             as a dummy array.
C
C     MNEI    (output) INTEGER array, dimension (3)
C             If MODE = 'B' or 'T' then
C             MNEI(1) contains the row dimension of
C                     sE(eps,inf)-A(eps,inf);
C             MNEI(2) contains the column dimension of
C                     sE(eps,inf)-A(eps,inf);
C             MNEI(3) = 0.
C             If MODE = 'S', then
C             MNEI(1) contains the row    dimension of sE(eps)-A(eps);
C             MNEI(2) contains the column dimension of sE(eps)-A(eps);
C             MNEI(3) contains the order of the regular pencil
C                     sE(inf)-A(inf).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance below which matrix elements are considered
C             to be zero. If the user sets TOL to be less than (or
C             equal to) zero then the tolerance is taken as
C             EPS * MAX( ABS(A(I,J)), ABS(E(I,J)) ), where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH),
C             I = 1,2,...,M and J = 1,2,...,N.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C             > 0:  if incorrect rank decisions were revealed during the
C                   triangularization phase. This failure is not likely
C                   to occur. The possible values are:
C             = 1:  if incorrect dimensions of a full column rank
C                   submatrix;
C             = 2:  if incorrect dimensions of a full row rank
C                   submatrix.
C
C     METHOD
C
C     Let sE - A be an arbitrary pencil. Prior to calling the routine,
C     this pencil must be transformed into a pencil with E in column
C     echelon form. This may be accomplished by calling the SLICOT
C     Library routine MB04UD. Depending on the value of MODE,
C     submatrices of A and E are then reduced to one of the forms
C     described above. Further details can be found in [1].
C
C     REFERENCES
C
C     [1] Beelen, Th. and Van Dooren, P.
C         An improved algorithm for the computation of Kronecker's
C         canonical form of a singular pencil.
C         Linear Algebra and Applications, 105, pp. 9-65, 1988.
C
C     NUMERICAL ASPECTS
C
C     It is shown in [1] that the algorithm is numerically backward
C     stable. The operations count is proportional to (MAX(M,N))**3.
C
C     FURTHER COMMENTS
C
C     The difference mu(k)-nu(k), for k = 1,2,...,NBLCKS, is the number
C     of elementary Kronecker blocks of size k x (k+1).
C
C     If MODE = 'B' or 'T' on entry, then the difference nu(k)-mu(k+1),
C     for k = 1,2,...,NBLCKS, is the number of infinite elementary
C     divisors of degree k (with mu(NBLCKS+1) = 0).
C
C     If MODE = 'S' on entry, then the difference mu0(k)-mu0(k+1),
C     for k = 1,2,...,NBLCKI, is the number of infinite elementary
C     divisors of degree k (with mu0(NBLCKI+1) = 0).
C     In the pencil sE(r)-A(r), the pencils sE(f)-A(f) and
C     sE(eta)-A(eta) can be separated by pertransposing the pencil
C     sE(r)-A(r) and calling the routine with MODE set to 'B'. The
C     result has got to be pertransposed again. (For more details see
C     [1]).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998.
C     Based on Release 3.0 routine MB04TD modified by A. Varga,
C     German Aerospace Research Establishment, Oberpfaffenhofen,
C     Germany, Nov. 1997, as follows:
C     1) NBLCKI is added;
C     2) the significance of IMUK0 and MNEI is changed;
C     3) INUK0 is removed.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, orthogonal transformation,
C     staircase form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBQ, JOBZ, MODE
      INTEGER           INFO, LDA, LDE, LDQ, LDZ, M, N, NBLCKI, NBLCKS,
     $                  RANKE
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IMUK(*), IMUK0(*), INUK(*), ISTAIR(*), IWORK(*),
     $                  MNEI(*)
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           FIRST, FIRSTI, LJOBQI, LJOBZI, LMODEB, LMODES,
     $                  LMODET, UPDATQ, UPDATZ
      INTEGER           I, IFICA, IFIRA, ISMUK, ISNUK, JK, K, NCA, NRA,
     $                  RANKA
      DOUBLE PRECISION  TOLER
C     .. Local Arrays ..
      DOUBLE PRECISION  DWORK(1)
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          DLASET, MB04TT, MB04TY, MB04VX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      INFO = 0
      LMODEB = LSAME( MODE, 'B' )
      LMODET = LSAME( MODE, 'T' )
      LMODES = LSAME( MODE, 'S' )
      LJOBQI = LSAME( JOBQ, 'I' )
      UPDATQ = LJOBQI.OR.LSAME( JOBQ, 'U' )
      LJOBZI = LSAME( JOBZ, 'I' )
      UPDATZ = LJOBZI.OR.LSAME( JOBZ, 'U' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LMODEB .AND. .NOT.LMODET .AND. .NOT.LMODES ) THEN
         INFO = -1
      ELSE IF( .NOT.UPDATQ .AND. .NOT.LSAME( JOBQ, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.UPDATZ .AND. .NOT.LSAME( JOBZ, 'N' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( RANKE.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDE.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( .NOT.UPDATQ .AND. LDQ.LT.1 .OR.
     $              UPDATQ .AND. LDQ.LT.MAX( 1, M ) ) THEN
         INFO = -12
      ELSE IF( .NOT.UPDATZ .AND. LDZ.LT.1 .OR.
     $              UPDATZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB04VD', -INFO )
         RETURN
      END IF
C
C     Initialize Q and Z to the identity matrices, if needed.
C
      IF ( LJOBQI )
     $   CALL DLASET( 'Full', M, M, ZERO, ONE, Q, LDQ )
      IF ( LJOBZI )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Quick return if possible.
C
      NBLCKS = 0
      NBLCKI = 0
C
      IF ( N.EQ.0 ) THEN
         MNEI(1) = 0
         MNEI(2) = 0
         MNEI(3) = 0
         RETURN
      END IF
C
      IF ( M.EQ.0 ) THEN
         NBLCKS = N
         DO 10 I = 1, N
            IMUK(I) = 1
            INUK(I) = 0
   10    CONTINUE
         MNEI(1) = 0
         MNEI(2) = N
         MNEI(3) = 0
         RETURN
      END IF
C
      TOLER = TOL
      IF ( TOLER.LE.ZERO )
     $   TOLER = DLAMCH( 'Epsilon' )*
     $                    MAX( DLANGE( 'M', M, N, A, LDA, DWORK ),
     $                         DLANGE( 'M', M, N, E, LDE, DWORK ) )
C
C     A(k) is the submatrix in A that will be row compressed.
C
C     ISMUK = sum(i=1,..,k) MU(i), ISNUK = sum(i=1,...,k) NU(i),
C     IFIRA, IFICA: first row and first column index of A(k) in A.
C     NRA, NCA: number of rows and columns in A(k).
C
      IFIRA = 1
      IFICA = 1
      NRA = M
      NCA = N - RANKE
      ISNUK = 0
      ISMUK = 0
      K = 0
C
C     Initialization of the arrays INUK and IMUK.
C
      DO 20 I = 1, M + 1
         INUK(I) = -1
   20 CONTINUE
C
C     Note: it is necessary that array INUK has DIMENSION M+1 since it
C           is possible that M = 1 and NBLCKS = 2.
C           Example sE-A = (0 0 s -1).
C
      DO 40 I = 1, N
         IMUK(I) = -1
   40 CONTINUE
C
C     Compress the rows of A while keeping E in column echelon form.
C
C     REPEAT
C
   60 K = K + 1
         CALL MB04TT( UPDATQ, UPDATZ, M, N, IFIRA, IFICA, NCA, A, LDA,
     $                E, LDE, Q, LDQ, Z, LDZ, ISTAIR, RANKA, TOLER,
     $                IWORK )
         IMUK(K) = NCA
         ISMUK = ISMUK + NCA
C
         INUK(K) = RANKA
         ISNUK  = ISNUK  + RANKA
         NBLCKS = NBLCKS + 1
C
C        If the rank of A(k) is nra then A has full row rank;
C        JK = the first column index (in A) after the right most column
C        of matrix A(k+1). (In case A(k+1) is empty, then JK = N+1.)
C
         IFIRA = 1 + ISNUK
         IFICA = 1 + ISMUK
         IF ( IFIRA.GT.M ) THEN
            JK = N + 1
         ELSE
            JK = ABS( ISTAIR(IFIRA) )
         END IF
         NRA = M - ISNUK
         NCA = JK - 1 - ISMUK
C
C        If NCA > 0 then there can be done some more row compression
C        of matrix A while keeping matrix E in column echelon form.
C
         IF ( NCA.GT.0 ) GO TO 60
C     UNTIL NCA <= 0
C
C     Matrix E(k+1) has full column rank since NCA = 0.
C     Reduce A and E by ignoring all rows and columns corresponding
C     to E(k+1). Ignoring these columns in E changes the ranks of the
C     submatrices E(i), (i=1,...,k-1).
C
      MNEI(1) = ISNUK
      MNEI(2) = ISMUK
      MNEI(3) = 0
C
      IF ( LMODEB )
     $   RETURN
C
C     Triangularization of the submatrices in A and E.
C
      CALL MB04TY( UPDATQ, UPDATZ, M, N, NBLCKS, INUK, IMUK, A, LDA, E,
     $             LDE, Q, LDQ, Z, LDZ, INFO )
C
      IF ( INFO.GT.0 .OR. LMODET )
     $   RETURN
C
C     Save the row dimensions of the diagonal submatrices in pencil
C     sE(eps,inf)-A(eps,inf).
C
      DO 80 I = 1, NBLCKS
         IMUK0(I) = INUK(I)
   80 CONTINUE
C
C     Reduction to square submatrices E(k)'s in E.
C
      CALL MB04VX( UPDATQ, UPDATZ, M, N, NBLCKS, INUK, IMUK, A, LDA, E,
     $             LDE, Q, LDQ, Z, LDZ, MNEI )
C
C     Determine the dimensions of the inf diagonal submatrices and
C     update block numbers if necessary.
C
      FIRST  = .TRUE.
      FIRSTI = .TRUE.
      NBLCKI = NBLCKS
      K = NBLCKS
C
      DO 100 I = K, 1, -1
         IMUK0(I) = IMUK0(I) - INUK(I)
         IF ( FIRSTI .AND. IMUK0(I).EQ.0 ) THEN
            NBLCKI = NBLCKI - 1
         ELSE
            FIRSTI = .FALSE.
         END IF
         IF ( FIRST .AND. IMUK(I).EQ.0 ) THEN
            NBLCKS = NBLCKS - 1
         ELSE
            FIRST = .FALSE.
         END IF
  100 CONTINUE
C
      RETURN
C *** Last line of MB04VD ***
      END
