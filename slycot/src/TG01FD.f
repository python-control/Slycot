      SUBROUTINE TG01FD( COMPQ, COMPZ, JOBA, L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, RANKE, RNKA22,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To compute for the descriptor system (A-lambda E,B,C)
C     the orthogonal transformation matrices Q and Z such that the
C     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is
C     in a SVD-like coordinate form with
C
C                  ( A11  A12 )             ( Er  0 )
C         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) ,
C                  ( A21  A22 )             (  0  0 )
C
C     where Er is an upper triangular invertible matrix.
C     Optionally, the A22 matrix can be further reduced to the form
C
C                  ( Ar  X )
C            A22 = (       ) ,
C                  (  0  0 )
C
C     with Ar an upper triangular invertible matrix, and X either a full
C     or a zero matrix.
C     The left and/or right orthogonal transformations performed
C     to reduce E and A22 can be optionally accumulated.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPQ   CHARACTER*1
C             = 'N':  do not compute Q;
C             = 'I':  Q is initialized to the unit matrix, and the
C                     orthogonal matrix Q is returned;
C             = 'U':  Q must contain an orthogonal matrix Q1 on entry,
C                     and the product Q1*Q is returned.
C
C     COMPZ   CHARACTER*1
C             = 'N':  do not compute Z;
C             = 'I':  Z is initialized to the unit matrix, and the
C                     orthogonal matrix Z is returned;
C             = 'U':  Z must contain an orthogonal matrix Z1 on entry,
C                     and the product Z1*Z is returned.
C
C     JOBA    CHARACTER*1
C             = 'N':  do not reduce A22;
C             = 'R':  reduce A22 to a SVD-like upper triangular form.
C             = 'T':  reduce A22 to an upper trapezoidal form.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of rows of matrices A, B, and E.  L >= 0.
C
C     N       (input) INTEGER
C             The number of columns of matrices A, E, and C.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of matrix C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix Q'*A*Z. If JOBA = 'T', this matrix
C             is in the form
C
C                           ( A11  *   *  )
C                  Q'*A*Z = (  *   Ar  X  ) ,
C                           (  *   0   0  )
C
C             where A11 is a RANKE-by-RANKE matrix and Ar is a
C             RNKA22-by-RNKA22 invertible upper triangular matrix.
C             If JOBA = 'R' then A has the above form with X = 0.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the descriptor matrix E.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix Q'*E*Z.
C
C                      ( Er  0 )
C             Q'*E*Z = (       ) ,
C                      (  0  0 )
C
C             where Er is a RANKE-by-RANKE upper triangular invertible
C             matrix.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,L).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading L-by-M part of this array must
C             contain the input/state matrix B.
C             On exit, the leading L-by-M part of this array contains
C             the transformed matrix Q'*B.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed matrix C*Z.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,L)
C             If COMPQ = 'N':  Q is not referenced.
C             If COMPQ = 'I':  on entry, Q need not be set;
C                              on exit, the leading L-by-L part of this
C                              array contains the orthogonal matrix Q,
C                              where Q' is the product of Householder
C                              transformations which are applied to A,
C                              E, and B on the left.
C             If COMPQ = 'U':  on entry, the leading L-by-L part of this
C                              array must contain an orthogonal matrix
C                              Q1;
C                              on exit, the leading L-by-L part of this
C                              array contains the orthogonal matrix
C                              Q1*Q.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= 1,        if COMPQ = 'N';
C             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If COMPZ = 'N':  Z is not referenced.
C             If COMPZ = 'I':  on entry, Z need not be set;
C                              on exit, the leading N-by-N part of this
C                              array contains the orthogonal matrix Z,
C                              which is the product of Householder
C                              transformations applied to A, E, and C
C                              on the right.
C             If COMPZ = 'U':  on entry, the leading N-by-N part of this
C                              array must contain an orthogonal matrix
C                              Z1;
C                              on exit, the leading N-by-N part of this
C                              array contains the orthogonal matrix
C                              Z1*Z.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.
C             LDZ >= 1,        if COMPZ = 'N';
C             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'.
C
C     RANKE   (output) INTEGER
C             The estimated rank of matrix E, and thus also the order
C             of the invertible upper triangular submatrix Er.
C
C     RNKA22  (output) INTEGER
C             If JOBA = 'R' or 'T', then RNKA22 is the estimated rank of
C             matrix A22, and thus also the order of the invertible
C             upper triangular submatrix Ar.
C             If JOBA = 'N', then RNKA22 is not referenced.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the rank of E
C             and of A22. If the user sets TOL > 0, then the given
C             value of TOL is used as a lower bound for the
C             reciprocal condition numbers of leading submatrices
C             of R or R22 in the QR decompositions E * P = Q * R of E
C             or A22 * P22 = Q22 * R22 of A22.
C             A submatrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = L*N*EPS,  is used instead, where
C             EPS is the machine precision (see LAPACK Library routine
C             DLAMCH). TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, N+P, MIN(L,N)+MAX(3*N-1,M,L) ).
C             For optimal performance, LDWORK should be larger.
C
C             If LDWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             DWORK array, returns this value as the first entry of
C             the DWORK array, and no error message related to LDWORK
C             is issued by XERBLA.
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
C     The routine computes a truncated QR factorization with column
C     pivoting of E, in the form
C
C                       ( E11 E12 )
C           E * P = Q * (         )
C                       (  0  E22 )
C
C     and finds the largest RANKE-by-RANKE leading submatrix E11 whose
C     estimated condition number is less than 1/TOL. RANKE defines thus
C     the rank of matrix E. Further E22, being negligible, is set to
C     zero, and an orthogonal matrix Y is determined such that
C
C           ( E11 E12 ) = ( Er  0 ) * Y .
C
C     The overal transformation matrix Z results as Z = P * Y' and the
C     resulting transformed matrices Q'*A*Z and Q'*E*Z have the form
C
C                          ( Er  0 )                      ( A11  A12 )
C         E <- Q'* E * Z = (       ) ,  A <- Q' * A * Z = (          ) ,
C                          (  0  0 )                      ( A21  A22 )
C
C     where Er is an upper triangular invertible matrix.
C     If JOBA = 'R' the same reduction is performed on A22 to obtain it
C     in the form
C
C                  ( Ar  0 )
C            A22 = (       ) ,
C                  (  0  0 )
C
C     with Ar an upper triangular invertible matrix.
C     If JOBA = 'T' then A22 is row compressed using the QR
C     factorization with column pivoting to the form
C
C                  ( Ar  X )
C            A22 = (       )
C                  (  0  0 )
C
C     with Ar an upper triangular invertible matrix.
C
C     The transformations are also applied to the rest of system
C     matrices
C
C          B <- Q' * B, C <- C * Z.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable and requires
C     0( L*L*N )  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the RASP routine RPDSSV.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 1999,
C     May 2003, Jan. 2009.
C
C     KEYWORDS
C
C     Descriptor system, matrix algebra, matrix operations,
C     orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOBA
      INTEGER            INFO, L, LDA, LDB, LDC, LDE, LDQ, LDWORK,
     $                   LDZ, M, N, P, RANKE, RNKA22
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK( * ),  E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      LOGICAL            ILQ, ILZ, LQUERY, REDA, REDTR, WITHB, WITHC
      INTEGER            I, ICOMPQ, ICOMPZ, IR1, IRE1, J, K, KW, LA22,
     $                   LH, LN, LWR, NA22, NB, WRKOPT
      DOUBLE PRECISION   SVLMAX, TOLDEF
C     .. Local Arrays ..
      DOUBLE PRECISION   SVAL(3)
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL           DLASET, DORMQR, DORMRZ, DSWAP, DTZRZF, MB03OY,
     $                   XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode COMPQ.
C
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'U' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
C
C     Decode COMPZ.
C
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'U' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
      REDA  = LSAME( JOBA, 'R' )
      REDTR = LSAME( JOBA, 'T' )
      WITHB = M.GT.0
      WITHC = P.GT.0
      LQUERY = ( LDWORK.EQ.-1 )
C
C     Test the input parameters.
C
      LN = MIN( L, N )
      INFO = 0
      WRKOPT = MAX( 1, N+P, LN + MAX( 3*N-1, M, L ) )
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( JOBA, 'N' ) .AND. .NOT.REDA .AND.
     $         .NOT.REDTR ) THEN
         INFO = -3
      ELSE IF( L.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -9
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.1 .OR. ( WITHB .AND. LDB.LT.L ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( ( ILQ .AND. LDQ.LT.L ) .OR. LDQ.LT.1 ) THEN
         INFO = -17
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -19
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -22
      ELSE
         IF( LQUERY ) THEN
            NB = MIN( 64, ILAENV( 1, 'DORMQR', 'LC', L, N, LN, -1 ) )
            WRKOPT = MAX( WRKOPT, LN + N*NB )
            IF( WITHB ) THEN
               NB = MIN( 64, ILAENV( 1, 'DORMQR', 'LC', L, M, LN, -1 ) )
               WRKOPT = MAX( WRKOPT, LN + M*NB )
            END IF
            IF( ILQ ) THEN
               NB = MIN( 64, ILAENV( 1, 'DORMQR', 'RN', L, L, LN, -1 ) )
               WRKOPT = MAX( WRKOPT, LN + L*NB )
            END IF
            NB = ILAENV( 1, 'DGERQF', ' ', L, N, -1, -1 )
            WRKOPT = MAX( WRKOPT, LN + N*NB )
            NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RC', L, N, N, -1 ) )
            WRKOPT = MAX( WRKOPT, N + MAX( 1, L )*NB )
            IF( WITHC ) THEN
               NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RC', P, N, N, -1 ) )
               WRKOPT = MAX( WRKOPT, N + MAX( 1, P )*NB )
            END IF
            IF( ILZ ) THEN
               NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RC', N, N, N, -1 ) )
               WRKOPT = MAX( WRKOPT, N + MAX( 1, N )*NB )
            END IF
         ELSE IF( LDWORK.LT.WRKOPT ) THEN
            INFO = -25
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01FD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Initialize Q and Z if necessary.
C
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', L, L, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Quick return if possible.
C
      IF( L.EQ.0 .OR. N.EQ.0 ) THEN
         DWORK(1) = ONE
         RANKE = 0
         IF( REDA .OR. REDTR ) RNKA22 = 0
         RETURN
      END IF
C
      TOLDEF = TOL
      IF( TOLDEF.LE.ZERO ) THEN
C
C        Use the default tolerance for rank determination.
C
         TOLDEF = DBLE( L*N )*DLAMCH( 'EPSILON' )
      END IF
C
C     Set the estimate of maximum singular value of E to
C     max(||E||,||A||) to detect negligible A or E matrices.
C
      SVLMAX = MAX( DLANGE( 'F', L, N, E, LDE, DWORK ),
     $              DLANGE( 'F', L, N, A, LDA, DWORK ) )
C
C     Compute the rank-revealing QR decomposition of E,
C
C                        ( E11 E12 )
C           E * P = Qr * (         ) ,
C                        (  0  E22 )
C
C     and determine the rank of E using incremental condition
C     estimation.
C     Workspace: MIN(L,N) + 3*N - 1.
C
      LWR = LDWORK - LN
      KW  = LN + 1
C
      CALL MB03OY( L, N, E, LDE, TOLDEF, SVLMAX, RANKE, SVAL, IWORK,
     $             DWORK, DWORK(KW), INFO )
C
C     Apply transformation on the rest of matrices.
C
      IF( RANKE.GT.0 ) THEN
C
C        A <-- Qr' * A.
C        Workspace: need   MIN(L,N) + N;
C                   prefer MIN(L,N) + N*NB.
C
         CALL DORMQR( 'Left', 'Transpose', L, N, RANKE, E, LDE, DWORK,
     $                A, LDA, DWORK(KW), LWR, INFO )
         WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C        B <-- Qr' * B.
C        Workspace: need   MIN(L,N) + M;
C                   prefer MIN(L,N) + M*NB.
C
         IF( WITHB ) THEN
            CALL DORMQR( 'Left', 'Transpose', L, M, RANKE, E, LDE,
     $                   DWORK, B, LDB, DWORK(KW), LWR, INFO )
            WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
         END IF
C
C        Q <-- Q * Qr.
C        Workspace: need   MIN(L,N) + L;
C                   prefer MIN(L,N) + L*NB.
C
         IF( ILQ ) THEN
            CALL DORMQR( 'Right', 'No Transpose', L, L, RANKE, E, LDE,
     $                   DWORK, Q, LDQ, DWORK(KW), LWR, INFO )
            WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
         END IF
C
C        Set lower triangle of E to zero.
C
         IF( L.GE.2 )
     $      CALL DLASET( 'Lower', L-1, RANKE, ZERO, ZERO, E(2,1), LDE )
C
C        Compute A*P, C*P and Z*P by forward permuting the columns of
C        A, C and Z based on information in IWORK.
C
         DO 10 J = 1, N
            IWORK(J) = -IWORK(J)
   10    CONTINUE
         DO 30 I = 1, N
            IF( IWORK(I).LT.0 ) THEN
               J = I
               IWORK(J) = -IWORK(J)
   20          CONTINUE
               K = IWORK(J)
               IF( IWORK(K).LT.0 ) THEN
                  CALL DSWAP( L, A(1,J), 1, A(1,K), 1 )
                  IF( WITHC )
     $               CALL DSWAP( P, C(1,J), 1, C(1,K), 1 )
                  IF( ILZ )
     $               CALL DSWAP( N, Z(1,J), 1, Z(1,K), 1 )
                  IWORK(K) = -IWORK(K)
                  J = K
                  GO TO 20
               END IF
            END IF
   30    CONTINUE
C
C        Determine an orthogonal matrix Y such that
C
C           ( E11 E12 ) = ( Er  0 ) * Y .
C
C        Compute E <-- E*Y', A <-- A*Y', C <-- C*Y', Z <-- Z*Y'.
C
         IF( RANKE.LT.N ) THEN
C
C           Workspace: need   2*N;
C                      prefer N + N*NB.
C
            KW = RANKE + 1
            CALL DTZRZF( RANKE, N, E, LDE, DWORK, DWORK(KW),
     $                   LDWORK-KW+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C           Workspace: need   N + MAX(L,P,N);
C                      prefer N + MAX(L,P,N)*NB.
C
            LH = N - RANKE
            CALL DORMRZ( 'Right', 'Transpose', L, N, RANKE, LH, E, LDE,
     $                   DWORK, A, LDA, DWORK(KW), LDWORK-KW+1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
            IF( WITHC ) THEN
               CALL DORMRZ( 'Right', 'Transpose', P, N, RANKE, LH, E,
     $                      LDE, DWORK, C, LDC, DWORK(KW), LDWORK-KW+1,
     $                      INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
            END IF
            IF( ILZ ) THEN
               CALL DORMRZ( 'Right', 'Transpose', N, N, RANKE, LH, E,
     $                      LDE, DWORK, Z, LDZ, DWORK(KW), LDWORK-KW+1,
     $                      INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
            END IF
C
C           Set E12 and E22 to zero.
C
            CALL DLASET( 'Full', L, LH, ZERO, ZERO, E(1,KW), LDE )
         END IF
      ELSE
         CALL DLASET( 'Full', L, N, ZERO, ZERO, E, LDE )
      END IF
C
C     Reduce A22 if necessary.
C
      IF( REDA .OR. REDTR ) THEN
         LA22 = L - RANKE
         NA22 = N - RANKE
         IF( MIN( LA22, NA22 ).EQ.0 ) THEN
            RNKA22 = 0
         ELSE
C
C           Compute the rank-revealing QR decomposition of A22,
C
C                              ( R11 R12 )
C              A22 * P2 = Q2 * (         ) ,
C                              (  0  R22 )
C
C           and determine the rank of A22 using incremental
C           condition estimation.
C           Workspace: MIN(L,N) + 3*N - 1.
C
            IR1 = RANKE + 1
            CALL MB03OY( LA22, NA22, A(IR1,IR1), LDA, TOLDEF,
     $                   SVLMAX, RNKA22, SVAL, IWORK, DWORK,
     $                   DWORK(KW), INFO )
C
C           Apply transformation on the rest of matrices.
C
            IF( RNKA22.GT.0 ) THEN
C
C              A <-- diag(I, Q2') * A
C              Workspace: need   MIN(L,N) + N;
C                         prefer MIN(L,N) + N*NB.
C
               CALL DORMQR( 'Left', 'Transpose', LA22, RANKE, RNKA22,
     $                      A(IR1,IR1), LDA, DWORK, A(IR1,1), LDA,
     $                      DWORK(KW), LWR, INFO )
C
C              B <-- diag(I, Q2') * B
C              Workspace: need   MIN(L,N) + M;
C                         prefer MIN(L,N) + M*NB.
C
               IF ( WITHB )
     $            CALL DORMQR( 'Left', 'Transpose', LA22, M, RNKA22,
     $                         A(IR1,IR1), LDA, DWORK, B(IR1,1), LDB,
     $                         DWORK(KW), LWR, INFO )
C
C              Q <-- Q * diag(I, Q2)
C              Workspace: need   MIN(L,N) + L;
C                         prefer MIN(L,N) + L*NB.
C
               IF( ILQ )
     $            CALL DORMQR( 'Right', 'No transpose', L, LA22, RNKA22,
     $                         A(IR1,IR1), LDA, DWORK, Q(1,IR1), LDQ,
     $                         DWORK(KW), LWR, INFO )
C
C              Set lower triangle of A22 to zero.
C
               IF( LA22.GE.2 )
     $            CALL DLASET( 'Lower', LA22-1, RNKA22, ZERO, ZERO,
     $                         A(IR1+1,IR1), LDA )
C
C              Compute A*diag(I,P2), C*diag(I,P2) and Z*diag(I,P2)
C              by forward permuting the columns of A, C and Z based
C              on information in IWORK.
C
               DO 40 J = 1, NA22
                  IWORK(J) = -IWORK(J)
   40          CONTINUE
               DO 60 I = 1, NA22
                  IF( IWORK(I).LT.0 ) THEN
                     J = I
                     IWORK(J) = -IWORK(J)
   50                CONTINUE
                     K = IWORK(J)
                     IF( IWORK(K).LT.0 ) THEN
                        CALL DSWAP( RANKE, A(1,RANKE+J), 1,
     $                              A(1,RANKE+K), 1 )
                        IF( WITHC )
     $                     CALL DSWAP( P, C(1,RANKE+J), 1,
     $                                 C(1,RANKE+K), 1 )
                        IF( ILZ )
     $                     CALL DSWAP( N, Z(1,RANKE+J), 1,
     $                                 Z(1,RANKE+K), 1 )
                        IWORK(K) = -IWORK(K)
                        J = K
                        GO TO 50
                     END IF
                  END IF
   60          CONTINUE
C
               IF( REDA .AND. RNKA22.LT.NA22 ) THEN
C
C                 Determine an orthogonal matrix Y2 such that
C
C                 ( R11 R12 ) = ( Ar  0 ) * Y2 .
C
C                 Compute A <-- A*diag(I, Y2'), C <-- C*diag(I, Y2'),
C                         Z <-- Z*diag(I, Y2').
C                 Workspace: need   2*N.
C                            prefer N + N*NB.
C
                  KW = RANKE + 1
                  CALL DTZRZF( RNKA22, NA22, A(IR1,IR1), LDA, DWORK,
     $                         DWORK(KW), LDWORK-KW+1, INFO )
                  WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C                 Workspace: need   N + MAX(P,N);
C                            prefer N + MAX(P,N)*NB.
C
                  LH = NA22 - RNKA22
                  IF( WITHC ) THEN
                     CALL DORMRZ( 'Right', 'Transpose', P, N, RNKA22,
     $                            LH, A(IR1,IR1), LDA, DWORK, C, LDC,
     $                            DWORK(KW), LDWORK-KW+1, INFO )
                     WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
                  END IF
                  IF( ILZ ) THEN
                     CALL DORMRZ( 'Right', 'Transpose', N, N, RNKA22,
     $                            LH, A(IR1,IR1), LDA, DWORK, Z, LDZ,
     $                            DWORK(KW), LDWORK-KW+1, INFO )
                     WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
                  END IF
                  IRE1 = RANKE + RNKA22 + 1
C
C                 Set R12 and R22 to zero.
C
                  CALL DLASET( 'Full', LA22, LH, ZERO, ZERO,
     $                         A(IR1,IRE1), LDA )
               END IF
            ELSE
               CALL DLASET( 'Full', LA22, NA22, ZERO, ZERO,
     $                      A(IR1,IR1), LDA)
            END IF
         END IF
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TG01FD ***
      END
