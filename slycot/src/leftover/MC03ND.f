      SUBROUTINE MC03ND( MP, NP, DP, P, LDP1, LDP2, DK, GAM, NULLSP,
     $                   LDNULL, KER, LDKER1, LDKER2, TOL, IWORK, DWORK,
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
C     To compute the coefficients of a minimal polynomial basis
C                                                 DK
C         K(s) = K(0) + K(1) * s + ... + K(DK) * s
C
C     for the right nullspace of the MP-by-NP polynomial matrix of
C     degree DP, given by
C                                                 DP
C         P(s) = P(0) + P(1) * s + ... + P(DP) * s  ,
C
C     which corresponds to solving the polynomial matrix equation
C     P(s) * K(s) = 0.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     MP      (input) INTEGER
C             The number of rows of the polynomial matrix P(s).
C             MP >= 0.
C
C     NP      (input) INTEGER
C             The number of columns of the polynomial matrix P(s).
C             NP >= 0.
C
C     DP      (input) INTEGER
C             The degree of the polynomial matrix P(s).  DP >= 1.
C
C     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1)
C             The leading MP-by-NP-by-(DP+1) part of this array must
C             contain the coefficients of the polynomial matrix P(s).
C             Specifically, P(i,j,k) must contain the (i,j)-th element
C             of P(k-1), which is the cofficient of s**(k-1) of P(s),
C             where i = 1,2,...,MP, j = 1,2,...,NP and k = 1,2,...,DP+1.
C
C     LDP1    INTEGER
C             The leading dimension of array P.  LDP1 >= MAX(1,MP).
C
C     LDP2    INTEGER
C             The second dimension of array P.   LDP2 >= MAX(1,NP).
C
C     DK      (output) INTEGER
C             The degree of the minimal polynomial basis K(s) for the
C             right nullspace of P(s) unless DK = -1, in which case
C             there is no right nullspace.
C
C     GAM     (output) INTEGER array, dimension (DP*MP+1)
C             The leading (DK+1) elements of this array contain
C             information about the ordering of the right nullspace
C             vectors stored in array NULLSP.
C
C     NULLSP  (output) DOUBLE PRECISION array, dimension
C             (LDNULL,(DP*MP+1)*NP)
C             The leading NP-by-SUM(i*GAM(i)) part of this array
C             contains the right nullspace vectors of P(s) in condensed
C             form (as defined in METHOD), where i = 1,2,...,DK+1.
C
C     LDNULL  INTEGER
C             The leading dimension of array NULLSP.
C             LDNULL >= MAX(1,NP).
C
C     KER     (output) DOUBLE PRECISION array, dimension
C             (LDKER1,LDKER2,DP*MP+1)
C             The leading NP-by-nk-by-(DK+1) part of this array contains
C             the coefficients of the minimal polynomial basis K(s),
C             where nk = SUM(GAM(i)) and i = 1,2,...,DK+1. Specifically,
C             KER(i,j,m) contains the (i,j)-th element of K(m-1), which
C             is the coefficient of s**(m-1) of K(s), where i = 1,2,...,
C             NP, j = 1,2,...,nk and m = 1,2,...,DK+1.
C
C     LDKER1  INTEGER
C             The leading dimension of array KER.  LDKER1 >= MAX(1,NP).
C
C     LDKER2  INTEGER
C             The second dimension of array KER.   LDKER2 >= MAX(1,NP).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance below which matrix elements are considered
C             to be zero. If the user sets TOL to be less than
C             10 * EPS * MAX( ||A|| , ||E|| ), then the tolerance is
C                                  F       F
C             taken as 10 * EPS * MAX( ||A|| , ||E|| ), where EPS is the
C                                           F       F
C             machine precision (see LAPACK Library Routine DLAMCH) and
C             A and E are matrices (as defined in METHOD).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (m+2*MAX(n,m+1)+n),
C             where m = DP*MP and n = (DP-1)*MP + NP.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= m*n*n + 2*m*n + 2*n*n.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C             > 0:  if incorrect rank decisions were taken during the
C                   computations. This failure is not likely to occur.
C                   The possible values are:
C                     k, 1 <= k <= DK+1, the k-th diagonal submatrix had
C                           not a full row rank;
C                     DK+2, if incorrect dimensions of a full column
C                           rank submatrix;
C                     DK+3, if incorrect dimensions of a full row rank
C                           submatrix.
C
C     METHOD
C
C     The computation of the right nullspace of the MP-by-NP polynomial
C     matrix P(s) of degree DP given by
C                                                  DP-1            DP
C        P(s) = P(0) + P(1) * s + ... + P(DP-1) * s     + P(DP) * s
C
C     is performed via the pencil s*E - A, associated with P(s), where
C
C            | I              |           | 0         -P(DP) |
C            |   .            |           | I .          .   |
C        A = |     .          |  and  E = |   . .        .   |.      (1)
C            |       .        |           |     . 0      .   |
C            |         I      |           |       I 0 -P(2)  |
C            |           P(0) |           |         I -P(1)  |
C
C     The pencil s*E - A is transformed by unitary matrices Q and Z such
C     that
C
C                     | sE(eps)-A(eps) |        X       |      X     |
C                     |----------------|----------------|------------|
C                     |        0       | sE(inf)-A(inf) |      X     |
C        Q'(s*E-A)Z = |=================================|============|.
C                     |                                 |            |
C                     |                0                | sE(r)-A(r) |
C
C     Since s*E(inf)-A(inf) and s*E(r)-A(r) have full column rank, the
C     minimal polynomial basis for the right nullspace of Q'(s*E-A)Z
C     (and consequently the basis for the right nullspace of s*E - A) is
C     completely determined by s*E(eps)-A(eps).
C
C     Let Veps(s) be a minimal polynomial basis for the right nullspace
C     of s*E(eps)-A(eps). Then
C
C                   | Veps(s) |
C        V(s) = Z * |---------|
C                   |    0    |
C
C     is a minimal polynomial basis for the right nullspace of s*E - A.
C     From the structure of s*E - A it can be shown that if V(s) is
C     partitioned as
C
C               | Vo(s) | (DP-1)*MP
C        V(s) = |------ |
C               | Ve(s) | NP
C
C     then the columns of Ve(s) form a minimal polynomial basis for the
C     right nullspace of P(s).
C
C     The vectors of Ve(s) are computed and stored in array NULLSP in
C     the following condensed form:
C
C        ||      ||      |      ||      |      |      ||      |     |
C        || U1,0 || U2,0 | U2,1 || U3,0 | U3,1 | U3,2 || U4,0 | ... |,
C        ||      ||      |      ||      |      |      ||      |     |
C
C     where Ui,j is an NP-by-GAM(i) matrix which contains the i-th block
C     of columns of K(j), the j-th coefficient of the polynomial matrix
C     representation for the right nullspace
C                                                  DK
C        K(s) = K(0) + K(1) * s + . . . + K(DK) * s  .
C
C     The coefficients K(0), K(1), ..., K(DK) are NP-by-nk matrices
C     given by
C
C        K(0)  = | U1,0 | U2,0 | U3,0 | . . .          | U(DK+1,0) |
C
C        K(1)  = |  0   | U2,1 | U3,1 | . . .          | U(DK+1,1) |
C
C        K(2)  = |  0   |  0   | U3,2 | . . .          | U(DK+1,2) |
C
C          .     .     .     .     .     .     .     .     .     .
C
C        K(DK) = |  0   |  0   |  0   | . . .    |  0  | U(DK+1,DK)|.
C
C     Note that the degree of K(s) satisfies the inequality DK <=
C     DP * MIN(MP,NP) and that the dimension of K(s) satisfies the
C     inequality (NP-MP) <= nk <= NP.
C
C     REFERENCES
C
C     [1] Beelen, Th.G.J.
C         New Algorithms for Computing the Kronecker structure of a
C         Pencil with Applications to Systems and Control Theory.
C         Ph.D.Thesis, Eindhoven University of Technology, 1987.
C
C     [2] Van Den Hurk, G.J.H.H.
C         New Algorithms for Solving Polynomial Matrix Problems.
C         Master's Thesis, Eindhoven University of Technology, 1987.
C
C     NUMERICAL ASPECTS
C
C     The algorithm used by the routine involves the construction of a
C     special block echelon form with pivots considered to be non-zero
C     when they are larger than TOL. These pivots are then inverted in
C     order to construct the columns of the kernel of the polynomial
C     matrix. If TOL is chosen to be too small then these inversions may
C     be sensitive whereas increasing TOL will make the inversions more
C     robust but will affect the block echelon form (and hence the
C     column degrees of the polynomial kernel). Furthermore, if the
C     elements of the computed polynomial kernel are large relative to
C     the polynomial matrix, then the user should consider trying
C     several values of TOL.
C
C     FURTHER COMMENTS
C
C     It also possible to compute a minimal polynomial basis for the
C     right nullspace of a pencil, since a pencil is a polynomial matrix
C     of degree 1. Thus for the pencil (s*E - A), the required input is
C     P(1)  = E and P(0) = -A.
C
C     The routine can also be used to compute a minimal polynomial
C     basis for the left nullspace of a polynomial matrix by simply
C     transposing P(s).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC03BD by A.J. Geurts and MC03BZ by
C     Th.G.J. Beelen, A.J. Geurts, and G.J.H.H. van den Hurk.
C
C     REVISIONS
C
C     Jan. 1998.
C
C     KEYWORDS
C
C     Echelon form, elementary polynomial operations, input output
C     description, polynomial matrix, polynomial operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TEN
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DK, DP, INFO, LDKER1, LDKER2, LDNULL, LDP1,
     $                  LDP2, LDWORK, MP, NP
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           GAM(*), IWORK(*)
      DOUBLE PRECISION  DWORK(*), KER(LDKER1,LDKER2,*),
     $                  NULLSP(LDNULL,*), P(LDP1,LDP2,*)
C     .. Local Scalars ..
      INTEGER           GAMJ, H, I, IDIFF, IFIR, J, JWORKA, JWORKE,
     $                  JWORKQ, JWORKV, JWORKZ, K, M, MUK, N, NBLCKS,
     $                  NBLCKI, NCA, NCV, NRA, NUK, RANKE, SGAMK, TAIL,
     $                  VC1, VR2
      DOUBLE PRECISION  TOLER
C     .. Local Arrays ..
      INTEGER           MNEI(3)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE, DLAPY2
      EXTERNAL          DLAMCH, DLANGE, DLAPY2
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLACPY, DLASET, MB04UD, MB04VD, MC03NX,
     $                  MC03NY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, SQRT
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      M = DP*MP
      H = M - MP
      N = H + NP
      INFO = 0
      IF( MP.LT.0 ) THEN
         INFO = -1
      ELSE IF( NP.LT.0 ) THEN
         INFO = -2
      ELSE IF( DP.LE.0 ) THEN
         INFO = -3
      ELSE IF( LDP1.LT.MAX( 1, MP ) ) THEN
         INFO = -5
      ELSE IF( LDP2.LT.MAX( 1, NP ) ) THEN
         INFO = -6
      ELSE IF( LDNULL.LT.MAX( 1, NP ) ) THEN
         INFO = -10
      ELSE IF( LDKER1.LT.MAX( 1, NP ) ) THEN
         INFO = -12
      ELSE IF( LDKER2.LT.MAX( 1, NP ) ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.( N*( M*N + 2*( M + N ) ) ) ) THEN
         INFO = -17
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC03ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MP.EQ.0 .OR. NP.EQ.0 ) THEN
         DK = -1
         RETURN
      END IF
C
      JWORKA = 1
      JWORKE = JWORKA + M*N
      JWORKZ = JWORKE + M*N
      JWORKV = JWORKZ + N*N
      JWORKQ = JWORKA
C
C     Construct the matrices A and E in the pencil s*E-A in (1).
C     Workspace:  2*M*N.
C
      CALL MC03NX( MP, NP, DP, P, LDP1, LDP2, DWORK(JWORKA), M,
     $             DWORK(JWORKE), M )
C
C     Computation of the tolerance.
C
      TOLER = MAX( DLANGE( 'F', M,  NP, DWORK(JWORKE+H*M), M, DWORK ),
     $             DLANGE( 'F', MP, NP, P, LDP1, DWORK ) )
      TOLER = TEN*DLAMCH( 'Epsilon' )
     $           *DLAPY2( TOLER, SQRT( DBLE( H ) ) )
      IF ( TOLER.LE.TOL ) TOLER = TOL
C
C     Reduction of E to column echelon form E0 = Q' x E x Z and
C     transformation of A, A0 = Q' x A x Z.
C     Workspace:  2*M*N + N*N + max(M,N).
C
      CALL MB04UD( 'No Q', 'Identity Z', M, N, DWORK(JWORKA), M,
     $             DWORK(JWORKE), M, DWORK(JWORKQ), M, DWORK(JWORKZ), N,
     $             RANKE, IWORK, TOLER, DWORK(JWORKV), INFO )
C
C     The contents of ISTAIR is transferred from MB04UD to MB04VD by
C     IWORK(i), i=1,...,M.
C     In the sequel the arrays IMUK and INUK are part of IWORK, namely:
C     IWORK(i), i = M+1,...,M+max(N,M+1), contains IMUK,
C     IWORK(i), i = M+max(N,M+1)+1,...,M+2*max(N,M+1), contains INUK.
C     IWORK(i), i = M+2*max(N,M+1)+1,...,M+2*max(N,M+1)+N, contains
C               IMUK0 (not needed), and is also used as workspace.
C
      MUK  = M + 1
      NUK  = MUK + MAX( N, M+1 )
      TAIL = NUK + MAX( N, M+1 )
C
      CALL MB04VD( 'Separation', 'No Q', 'Update Z', M, N, RANKE,
     $             DWORK(JWORKA), M, DWORK(JWORKE), M, DWORK(JWORKQ), M,
     $             DWORK(JWORKZ), N, IWORK, NBLCKS, NBLCKI, IWORK(MUK),
     $             IWORK(NUK), IWORK(TAIL), MNEI, TOLER, IWORK(TAIL),
     $             INFO )
      IF ( INFO.GT.0 ) THEN
C
C        Incorrect rank decisions.
C
         INFO = INFO + NBLCKS
         RETURN
      END IF
C
C     If NBLCKS < 1, or the column dimension of s*E(eps) - A(eps) is
C     zero, then there is no right nullspace.
C
      IF ( NBLCKS.LT.1 .OR. MNEI(2).EQ.0 ) THEN
         DK = -1
         RETURN
      END IF
C
C     Start of the computation of the minimal basis.
C
      DK = NBLCKS - 1
      NRA = MNEI(1)
      NCA = MNEI(2)
C
C     Determine a minimal basis VEPS(s) for the right nullspace of the
C     pencil s*E(eps)-A(eps) associated with the polynomial matrix P(s).
C     Workspace:  2*M*N + N*N + N*N*(M+1).
C
      CALL MC03NY( NBLCKS, NRA, NCA, DWORK(JWORKA), M, DWORK(JWORKE), M,
     $             IWORK(MUK), IWORK(NUK), DWORK(JWORKV), N, INFO )
C
      IF ( INFO.GT.0 )
     $   RETURN
C
      NCV = IWORK(MUK) - IWORK(NUK)
      GAM(1) = NCV
      IWORK(1) = 0
      IWORK(TAIL) = IWORK(MUK)
C
      DO 20 I = 2, NBLCKS
         IDIFF = IWORK(MUK+I-1) - IWORK(NUK+I-1)
         GAM(I) = IDIFF
         IWORK(I) = NCV
         NCV = NCV + I*IDIFF
         IWORK(TAIL+I-1) = IWORK(TAIL+I-2) + IWORK(MUK+I-1)
   20 CONTINUE
C
C     Determine a basis for the right nullspace of the polynomial
C     matrix P(s). This basis is stored in array NULLSP in condensed
C     form.
C
      CALL DLASET( 'Full', NP, NCV, ZERO, ZERO, NULLSP, LDNULL )
C
C                                                |VEPS(s)|
C     The last NP rows of the product matrix Z x |-------| contain the
C                                                |   0   |
C     polynomial basis for the right nullspace of the polynomial matrix
C     P(s) in condensed form. The multiplication is restricted to the
C     nonzero submatrices Vij,k of VEPS, the result is stored in the
C     array NULLSP.
C
      VC1 = 1
C
      DO 60 I = 1, NBLCKS
         VR2 = IWORK(TAIL+I-1)
C
         DO 40 J = 1, I
C
C           Multiplication of Z(H+1:N,1:VR2) with V.i,j-1 stored in
C           VEPS(1:VR2,VC1:VC1+GAM(I)-1).
C
            CALL DGEMM( 'No transpose', 'No transpose', NP, GAM(I), VR2,
     $                  ONE, DWORK(JWORKZ+H), N,
     $                  DWORK(JWORKV+(VC1-1)*N), N, ZERO, NULLSP(1,VC1),
     $                  LDNULL )
            VC1 = VC1 + GAM(I)
            VR2 = VR2 - IWORK(MUK+I-J)
   40    CONTINUE
C
   60 CONTINUE
C
C     Transfer of the columns of NULLSP to KER in order to obtain the
C     polynomial matrix representation of K(s), the right nullspace
C     of P(s).
C
      SGAMK = 1
C
      DO 100 K = 1, NBLCKS
         CALL DLASET( 'Full', NP, SGAMK-1, ZERO, ZERO, KER(1,1,K),
     $                LDKER1 )
         IFIR = SGAMK
C
C        Copy the appropriate columns of NULLSP into KER(k).
C        SGAMK = 1 + SUM(i=1,..,k-1) GAM(i), is the first nontrivial
C        column of KER(k), the first SGAMK - 1 columns of KER(k) are
C        zero. IFIR denotes the position of the first column in KER(k)
C        in the set of columns copied for a value of J.
C        VC1 is the first column of NULLSP to be copied.
C
         DO 80 J = K, NBLCKS
            GAMJ = GAM(J)
            VC1 = IWORK(J) + (K-1)*GAMJ + 1
            CALL DLACPY( 'Full', NP, GAMJ, NULLSP(1,VC1), LDNULL,
     $                   KER(1,IFIR,K), LDKER1 )
            IFIR = IFIR + GAMJ
   80    CONTINUE
C
         SGAMK = SGAMK + GAM(K)
  100 CONTINUE
C
      RETURN
C *** Last line of MC03ND ***
      END
