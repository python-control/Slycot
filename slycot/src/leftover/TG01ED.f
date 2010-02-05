      SUBROUTINE TG01ED( JOBA, L, N, M, P, A, LDA, E, LDE, B, LDB,
     $                   C, LDC, Q, LDQ, Z, LDZ, RANKE, RNKA22, TOL,
     $                   DWORK, LDWORK, INFO )
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
C     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an
C     SVD (singular value decomposition) coordinate form with
C     the system matrices  Q'*A*Z and Q'*E*Z in the form
C
C                  ( A11  A12 )             ( Er  0 )
C         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) ,
C                  ( A21  A22 )             (  0  0 )
C
C     where Er is an invertible diagonal matrix having on the diagonal
C     the decreasingly ordered nonzero singular values of E.
C     Optionally, the A22 matrix can be further reduced to the
C     SVD form
C
C                  ( Ar  0 )
C            A22 = (       ) ,
C                  (  0  0 )
C
C     where Ar is an invertible diagonal matrix having on the diagonal
C     the decreasingly ordered nonzero singular values of A22.
C     The left and/or right orthogonal transformations performed
C     to reduce E and A22 are accumulated.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBA    CHARACTER*1
C             = 'N':  do not reduce A22;
C             = 'R':  reduce A22 to an SVD form.
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
C             the transformed matrix Q'*A*Z. If JOBA = 'R', this matrix
C             is in the form
C
C                           ( A11  *   *  )
C                  Q'*A*Z = (  *   Ar  0  ) ,
C                           (  *   0   0  )
C
C             where A11 is a RANKE-by-RANKE matrix and Ar is a
C             RNKA22-by-RNKA22 invertible diagonal matrix, with
C             decresingly ordered positive diagonal elements.
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
C             where Er is a RANKE-by-RANKE invertible diagonal matrix
C             having on the diagonal the decreasingly ordered positive
C             singular values of E.
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
C     Q       (output) DOUBLE PRECISION array, dimension (LDQ,L)
C             The leading L-by-L part of this array contains the
C             orthogonal matrix Q, which is the accumulated product of
C             transformations applied to A, E, and B on the left.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,L).
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             The leading N-by-N part of this array contains the
C             orthogonal matrix Z, which is the accumulated product of
C             transformations applied to A, E, and C on the right.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= MAX(1,N).
C
C     RANKE   (output) INTEGER
C             The effective rank of matrix E, and thus also the order
C             of the invertible diagonal submatrix Er.
C             RANKE is computed as the number of singular values of E
C             greater than TOL*SVEMAX, where SVEMAX is the maximum
C             singular value of E.
C
C     RNKA22  (output) INTEGER
C             If JOBA = 'R', then RNKA22 is the effective rank of
C             matrix A22, and thus also the order of the invertible
C             diagonal submatrix Ar. RNKA22 is computed as the number
C             of singular values of A22 greater than TOL*SVAMAX,
C             where SVAMAX is an estimate of the maximum singular value
C             of A.
C             If JOBA = 'N', then RNKA22 is not referenced.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the rank of E
C             and of A22. If TOL > 0, then singular values less than
C             TOL*SVMAX are treated as zero, where SVMAX is the maximum
C             singular value of E or an estimate of it for A and E.
C             If TOL <= 0, the default tolerance TOLDEF = EPS*L*N is
C             used instead, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH). TOL < 1.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,MIN(L,N) +
C                           MAX(3*MIN(L,N)+MAX(L,N), 5*MIN(L,N), M, P)).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  the QR algorithm has failed to converge when computing
C                   singular value decomposition. In this case INFO
C                   specifies how many superdiagonals did not converge.
C                   This failure is not likely to occur.
C
C     METHOD
C
C     The routine computes the singular value decomposition (SVD) of E,
C     in the form
C
C                    ( Er  0 )
C           E  = Q * (       ) * Z'
C                    (  0  0 )
C
C     and finds the largest RANKE-by-RANKE leading diagonal submatrix
C     Er whose condition number is less than 1/TOL. RANKE defines thus
C     the effective rank of matrix E.
C     If JOBA = 'R' the same reduction is performed on A22 in the
C     partitioned matrix
C
C                  ( A11  A12 )
C         Q'*A*Z = (          ) ,
C                  ( A21  A22 )
C
C     to obtain it in the form
C
C                  ( Ar  0 )
C            A22 = (       ) ,
C                  (  0  0 )
C
C     with Ar an invertible diagonal matrix.
C
C     The accumulated transformations are also applied to the rest of
C     matrices
C
C          B <- Q' * B,  C <- C * Z.
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
C     Feb. 2000, Oct. 2001, May 2003.
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
      CHARACTER          JOBA
      INTEGER            INFO, L, LDA, LDB, LDC, LDE, LDQ, LDWORK,
     $                   LDZ, M, N, P, RNKA22, RANKE
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK( * ),  E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      LOGICAL            REDA
      INTEGER            I, IR1, J, KW, LA22, LN, LN2, LWR, NA22, WRKOPT
      DOUBLE PRECISION   EPSM, SVEMAX, SVLMAX, TOLDEF
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DGEQRF, DGELQF, DGESVD,
     $                   DLACPY, DLASET, DORMQR, DORMLQ, DSWAP, MA02AD,
     $                   MB03UD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
      REDA = LSAME( JOBA, 'R' )
C
C     Test the input parameters.
C
      INFO = 0
      WRKOPT = MIN( L, N ) +
     $         MAX( M, P, 3*MIN( L, N ) + MAX( L, N ), 5*MIN( L, N ) )
      IF( .NOT.LSAME( JOBA, 'N' ) .AND. .NOT.REDA ) THEN
         INFO = -1
      ELSE IF( L.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -7
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.L ) ) THEN
         INFO = -11
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDQ.LT.MAX( 1, L ) ) THEN
         INFO = -15
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -20
      ELSE IF( LDWORK.LT.MAX( 1, WRKOPT ) ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01ED', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( L.EQ.0 .OR. N.EQ.0 ) THEN
         IF( L.GT.0 )
     $      CALL DLASET( 'Full', L, L, ZERO, ONE, Q, LDQ )
         IF( N.GT.0 )
     $      CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         DWORK(1) = ONE
         RANKE = 0
         IF( REDA ) RNKA22 = 0
         RETURN
      END IF
C
      LN = MIN( L, N )
      EPSM = DLAMCH( 'EPSILON' )
C
      TOLDEF = TOL
      IF( TOLDEF.LE.ZERO ) THEN
C
C        Use the default tolerance for rank determination.
C
         TOLDEF = EPSM * DBLE( L*N )
      END IF
C
C     Set the estimate of the maximum singular value of E to
C     max(||E||,||A||) to detect negligible A or E matrices.
C
      SVLMAX = MAX( DLANGE( 'F', L, N, E, LDE, DWORK ) ,
     $              DLANGE( 'F', L, N, A, LDA, DWORK ) )
C
C     Compute the SVD of E
C
C                    ( Er  0 )
C           E = Qr * (       ) * Zr'
C                    (  0  0 )
C
C     Workspace: needed  MIN(L,N) + MAX(3*MIN(L,N)+MAX(L,N),5*MIN(L,N));
C                prefer larger.
C
      LWR = LDWORK - LN
      KW = LN + 1
C
      CALL DGESVD( 'A', 'A', L, N, E, LDE, DWORK, Q, LDQ, Z, LDZ,
     $             DWORK(KW), LWR, INFO )
      IF( INFO.GT.0 )
     $   RETURN
      WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C     Determine the rank of E.
C
      RANKE = 0
      IF( DWORK(1).GT.SVLMAX*EPSM ) THEN
         RANKE = 1
         SVEMAX = DWORK(1)
         DO 10 I = 2, LN
            IF( DWORK(I).LT.SVEMAX*TOLDEF ) GO TO 20
            RANKE = RANKE + 1
   10    CONTINUE
C
   20    CONTINUE
      END IF
C
C     Apply transformation on the rest of matrices.
C
      IF( RANKE.GT.0 ) THEN
C
C        A <-- Qr' * A * Zr.
C
         CALL DGEMM( 'Transpose', 'No transpose', L, N, L, ONE,
     $               Q, LDQ, A, LDA, ZERO, E, LDE )
         CALL DGEMM( 'No transpose', 'Transpose', L, N, N, ONE,
     $               E, LDE, Z, LDZ, ZERO, A, LDA )
C
C        B <-- Qr' * B.
C        Workspace: need   L;
C                   prefer L*M.
C
         IF( LWR.GT.L*M .AND. M.GT.0 ) THEN
C
            CALL DGEMM( 'Transpose', 'No transpose', L, M, L, ONE,
     $                  Q, LDQ, B, LDB, ZERO, DWORK(KW), L )
            CALL DLACPY( 'Full', L, M, DWORK(KW), L, B, LDB )
         ELSE
            DO 30 J = 1, M
               CALL DGEMV( 'Transpose', L, L, ONE, Q, LDQ, B(1,J), 1,
     $                     ZERO, DWORK(KW), 1 )
               CALL DCOPY( L, DWORK(KW), 1, B(1,J), 1 )
   30       CONTINUE
         END IF
C
C        C <-- C * Zr.
C        Workspace: need   N;
C                   prefer P*N.
C
         IF( LWR.GT.P*N ) THEN
C
            CALL DGEMM( 'No transpose', 'Transpose', P, N, N, ONE,
     $                  C, LDC, Z, LDZ, ZERO, DWORK(KW), MAX( 1, P ) )
            CALL DLACPY( 'Full', P, N, DWORK(KW), MAX( 1, P ), C, LDC )
         ELSE
            DO 40 I = 1, P
               CALL DGEMV( 'No transpose', N, N, ONE, Z, LDZ,
     $                     C(I,1), LDC, ZERO, DWORK(KW), 1 )
               CALL DCOPY( N, DWORK(KW), 1, C(I,1), LDC )
   40       CONTINUE
         END IF
         WRKOPT = MAX( WRKOPT, L*M, P*N )
      END IF
C
C     Reduce A22 if necessary.
C
      IF( REDA ) THEN
         LA22 = L - RANKE
         NA22 = N - RANKE
         LN2 = MIN( LA22, NA22 )
         IF( LN2.EQ.0 ) THEN
            IR1 = 1
            RNKA22 = 0
         ELSE
C
C           Compute the SVD of A22 using a storage saving approach.
C
            IR1 = RANKE + 1
            IF( LA22.GE.NA22 ) THEN
C
C              Compute the QR decomposition of A22 in the form
C
C              A22 = Q2 * ( R2 ) ,
C                         ( 0  )
C
C              where R2 is upper triangular.
C              Workspace: need   MIN(L,N) + N;
C                         prefer MIN(L,N) + N*NB.
C
               CALL DGEQRF( LA22, NA22, A(IR1,IR1), LDA, DWORK(IR1),
     $                      DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Apply transformation Q2 to A, B, and Q.
C
C              A <--diag(I, Q2') * A
C              Workspace: need   MIN(L,N) + N;
C                         prefer MIN(L,N) + N*NB.
C
               CALL DORMQR( 'Left', 'Transpose', LA22, RANKE, LN2,
     $                      A(IR1,IR1), LDA, DWORK(IR1), A(IR1,1), LDA,
     $                      DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              B <-- diag(I, Q2') * B
C              Workspace: need   MIN(L,N) + M;
C                         prefer MIN(L,N) + M*NB.
C
               IF ( M.GT.0 ) THEN
                  CALL DORMQR( 'Left', 'Transpose', LA22, M, LN2,
     $                         A(IR1,IR1), LDA, DWORK(IR1), B(IR1,1),
     $                         LDB, DWORK(KW), LWR, INFO )
                  WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
               END IF
C
C              Q <-- Q * diag(I, Q2)
C              Workspace: need   MIN(L,N) + L;
C                         prefer MIN(L,N) + L*NB.
C
               CALL DORMQR( 'Right', 'No transpose', L, LA22, LN2,
     $                       A(IR1,IR1), LDA, DWORK(IR1), Q(1,IR1), LDQ,
     $                       DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Compute the SVD of the upper triangular submatrix R2 as
C
C                               ( Ar  0 )
C                    R2 = Q2r * (       ) * Z2r' ,
C                               (  0  0 )
C
C              where Q2r is stored in E and Z2r' is stored in A22.
C              Workspace: need   MAX(1,5*MIN(L,N));
C                         prefer larger.
C
               CALL MB03UD( 'Vectors', 'Vectors', LN2, A(IR1,IR1), LDA,
     $                      E(IR1,IR1), LDE, DWORK(IR1), DWORK(KW), LWR,
     $                      INFO )
               IF( INFO.GT.0 )
     $            RETURN
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Determine the rank of A22.
C
               RNKA22 = 0
               IF( DWORK(IR1).GT.SVLMAX*EPSM ) THEN
                  RNKA22 = 1
                  DO 50 I = IR1+1, LN
                     IF( DWORK(I).LE.SVLMAX*TOLDEF ) GO TO 60
                     RNKA22 = RNKA22 + 1
   50             CONTINUE
C
   60             CONTINUE
               END IF
C
C              Apply transformation on the rest of matrices.
C
               IF( RNKA22.GT.0 ) THEN
C
C                 A <-- diag(I,Q2r') * A * diag(I,Zr2)
C
                  CALL DGEMM( 'Transpose', 'No transpose', LN2, RANKE,
     $                        LN2, ONE, E(IR1,IR1), LDE, A(IR1,1), LDA,
     $                        ZERO, E(IR1,1), LDE )
                  CALL DLACPY( 'Full', LN2, RANKE, E(IR1,1), LDE,
     $                         A(IR1,1), LDA )
                  CALL DGEMM( 'No transpose', 'Transpose', RANKE, LN2,
     $                        LN2, ONE, A(1,IR1), LDA, A(IR1,IR1), LDA,
     $                        ZERO, E(1,IR1), LDE )
                  CALL DLACPY( 'Full', RANKE, LN2, E(1,IR1), LDE,
     $                         A(1,IR1), LDA )
C
C                 B <-- diag(I,Q2r') * B
C
                  IF( LWR.GT.LN2*M .AND. M.GT.0 ) THEN
C
                     CALL DGEMM( 'Transpose', 'No transpose', LN2, M,
     $                           LN2, ONE, E(IR1,IR1), LDE, B(IR1,1),
     $                           LDB, ZERO, DWORK(KW), LN2 )
                     CALL DLACPY( 'Full', LN2, M, DWORK(KW), LN2,
     $                            B(IR1,1), LDB )
                  ELSE
                     DO 70 J = 1, M
                        CALL DGEMV( 'Transpose', LN2, LN2, ONE,
     $                              E(IR1,IR1), LDE, B( IR1,J), 1,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, B(IR1,J), 1 )
   70                CONTINUE
                  END IF
C
C                 C <-- C * diag(I,Zr2)
C
                  IF( LWR.GT.P*LN2 .AND. P.GT.0 ) THEN
C
                     CALL DGEMM( 'No transpose', 'Transpose', P, LN2,
     $                           LN2, ONE, C(1,IR1), LDC, A(IR1,IR1),
     $                           LDA, ZERO, DWORK(KW), P )
                     CALL DLACPY( 'Full', P, LN2, DWORK( KW ), P,
     $                            C(1,IR1), LDC )
                  ELSE
                     DO 80 I = 1, P
                        CALL DGEMV( 'No transpose', LN2, LN2, ONE,
     $                              A(IR1,IR1), LDA, C(I,IR1), LDC,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, C(I,IR1), LDC )
   80                CONTINUE
                  END IF
C
C                 Q <-- Q * diag(I, Qr2)
C
                  IF( LWR.GT.L*LN2 ) THEN
C
                     CALL DGEMM( 'No transpose', 'No transpose', L, LN2,
     $                           LN2, ONE, Q(1,IR1), LDQ, E(IR1,IR1),
     $                           LDE, ZERO, DWORK(KW), L )
                     CALL DLACPY( 'Full', L, LN2, DWORK(KW), L,
     $                            Q(1,IR1), LDQ )
                  ELSE
                     DO 90 I = 1, L
                        CALL DGEMV( 'Transpose', LN2, LN2, ONE,
     $                              E(IR1,IR1), LDE, Q(I,IR1), LDQ,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, Q(I,IR1), LDQ )
   90                CONTINUE
                  END IF
C
C                 Z' <-- diag(I, Zr2') * Z'
C
                  IF( LWR.GT.N*LN2 ) THEN
C
                     CALL DGEMM( 'No transpose', 'No transpose', LN2, N,
     $                           LN2, ONE, A(IR1,IR1), LDA, Z(IR1,1),
     $                           LDZ, ZERO, DWORK(KW), LN2 )
                     CALL DLACPY( 'Full', LN2, N, DWORK(KW), LN2,
     $                            Z(IR1,1), LDZ )
                  ELSE
                     DO 100 J = 1, N
                        CALL DGEMV( 'No transpose', LN2, LN2, ONE,
     $                              A(IR1,IR1), LDA, Z(IR1,J), 1,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, Z(IR1,J), 1 )
  100                CONTINUE
                  END IF
               END IF
            ELSE
C
C              Compute the LQ decomposition of A22 in the form
C
C                  A22 = ( L2 0 )* Z2
C
C              where L2 is lower triangular.
C              Workspace: need   MIN(L,N) + L;
C                         prefer MIN(L,N) + L*NB.
C
               CALL DGELQF( LA22, NA22, A(IR1,IR1), LDA, DWORK(IR1),
     $                      DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Apply transformation Z2 to A, C, and Z.
C
C              A <-- A * diag(I, Z2')
C              Workspace: need   2*MIN(L,N);
C                         prefer MIN(L,N) + MIN(L,N)*NB.
C
               CALL DORMLQ( 'Right', 'Transpose', RANKE, NA22, LN2,
     $                      A(IR1,IR1), LDA, DWORK(IR1), A(1,IR1), LDA,
     $                      DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              C <-- C * diag(I, Z2')
C              Workspace: need   MIN(L,N) + P;
C                         prefer MIN(L,N) + P*NB.
C
               IF ( P.GT.0 ) THEN
                  CALL DORMLQ( 'Right', 'Transpose', P, NA22, LN2,
     $                         A(IR1,IR1), LDA, DWORK(IR1), C(1,IR1),
     $                         LDC, DWORK(KW), LWR, INFO )
                  WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
               END IF
C
C              Z' <-  diag(I, Z2) * Z'
C              Workspace: need   MIN(L,N) + N;
C                         prefer MIN(L,N) + N*NB.
C
               CALL DORMLQ( 'Left', 'No transpose', NA22, N, LN2,
     $                       A(IR1,IR1), LDA, DWORK(IR1), Z(IR1,1), LDZ,
     $                       DWORK(KW), LWR, INFO )
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Compute the SVD of the lower triangular submatrix L2 as
C
C                              ( Ar  0 )
C                  L2' = Z2r * (       ) * Q2r'
C                              (  0  0 )
C
C              where Q2r' is stored in E and Z2r is stored in A22.
C              Workspace: need   MAX(1,5*MIN(L,N));
C                         prefer larger.
C
               CALL MA02AD( 'Lower', LN2, LN2, A(IR1,IR1), LDA,
     $                      E(IR1,IR1), LDE )
               CALL MB03UD( 'Vectors', 'Vectors', LN2, E(IR1,IR1), LDE,
     $                      A(IR1,IR1), LDA, DWORK(IR1), DWORK(KW),
     $                      LWR, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               WRKOPT = MAX( WRKOPT, LN + INT( DWORK(KW) ) )
C
C              Determine the rank of A22.
C
               RNKA22 = 0
               IF( DWORK(IR1).GT.SVLMAX*EPSM ) THEN
                  RNKA22 = 1
                  DO 110 I = IR1+1, LN
                     IF( DWORK(I).LE.SVLMAX*TOLDEF ) GO TO 120
                     RNKA22 = RNKA22 + 1
  110             CONTINUE
C
  120             CONTINUE
               END IF
C
C              Apply transformation on the rest of matrices.
C
               IF( RNKA22.GT.0 ) THEN
C
C                 A <-- diag(I,Q2r') * A * diag(I,Zr2)
C
                  CALL DGEMM( 'No transpose', 'No transpose', LN2,
     $                        RANKE, LN2, ONE, E(IR1,IR1), LDE,
     $                        A(IR1,1), LDA, ZERO, E(IR1,1), LDE )
                  CALL DLACPY( 'Full', LN2, RANKE, E(IR1,1), LDE,
     $                         A(IR1,1), LDA )
                  CALL DGEMM( 'No transpose', 'No transpose', RANKE,
     $                        LN2, LN2, ONE, A(1,IR1), LDA,
     $                        A(IR1,IR1), LDA, ZERO, E(1,IR1), LDE )
                  CALL DLACPY( 'Full', RANKE, LN2, E(1,IR1), LDE,
     $                         A(1,IR1), LDA )
C
C                 B <-- diag(I,Q2r') * B
C
                  IF( LWR.GT.LN2*M .AND. M.GT.0 ) THEN
C
                     CALL DGEMM( 'No transpose', 'No transpose', LN2, M,
     $                           LN2, ONE, E(IR1,IR1), LDE, B(IR1,1),
     $                           LDB, ZERO, DWORK(KW), LN2 )
                     CALL DLACPY( 'Full', LN2, M, DWORK(KW), LN2,
     $                            B(IR1,1), LDB )
                  ELSE
                     DO 130 J = 1, M
                        CALL DGEMV( 'No transpose', LN2, LN2, ONE,
     $                              E(IR1,IR1), LDE, B( IR1,J), 1,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, B(IR1,J), 1 )
  130                CONTINUE
                  END IF
C
C                 C <-- C * diag(I,Zr2)
C
                  IF( LWR.GT.P*LN2 .AND. P.GT.0 ) THEN
C
                     CALL DGEMM( 'No transpose', 'No transpose', P, LN2,
     $                           LN2, ONE, C(1,IR1), LDC, A(IR1,IR1),
     $                           LDA, ZERO, DWORK(KW), P )
                     CALL DLACPY( 'Full', P, LN2, DWORK( KW ), P,
     $                            C(1,IR1), LDC )
                  ELSE
                     DO 140 I = 1, P
                        CALL DGEMV( 'Transpose', LN2, LN2, ONE,
     $                              A(IR1,IR1), LDA, C(I,IR1), LDC,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, C(I,IR1), LDC )
  140                CONTINUE
                  END IF
C
C                 Q <-- Q * diag(I, Qr2)
C
                  IF( LWR.GT.L*LN2 ) THEN
C
                     CALL DGEMM( 'No transpose', 'Transpose', L, LN2,
     $                           LN2, ONE, Q(1,IR1), LDQ, E(IR1,IR1),
     $                           LDE, ZERO, DWORK(KW), L )
                     CALL DLACPY( 'Full', L, LN2, DWORK(KW), L,
     $                            Q(1,IR1), LDQ )
                  ELSE
                     DO 150 I = 1, L
                        CALL DGEMV( 'No transpose', LN2, LN2, ONE,
     $                              E(IR1,IR1), LDE, Q(I,IR1), LDQ,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, Q(I,IR1), LDQ )
  150                CONTINUE
                  END IF
C
C                 Z' <-- diag(I, Zr2') * Z'
C
                  IF( LWR.GT.N*LN2 ) THEN
C
                     CALL DGEMM( 'Transpose', 'No transpose', LN2, N,
     $                           LN2, ONE, A(IR1,IR1), LDA, Z(IR1,1),
     $                           LDZ, ZERO, DWORK(KW), LN2 )
                     CALL DLACPY( 'Full', LN2, N, DWORK(KW), LN2,
     $                            Z(IR1,1), LDZ )
                  ELSE
                     DO 160 J = 1, N
                        CALL DGEMV( 'Transpose', LN2, LN2, ONE,
     $                              A(IR1,IR1), LDA, Z(IR1,J), 1,
     $                              ZERO, DWORK(KW), 1 )
                        CALL DCOPY( LN2, DWORK(KW), 1, Z(IR1,J), 1 )
  160                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     Set E.
C
      CALL DLASET( 'Full', L, N, ZERO, ZERO, E, LDE )
      CALL DCOPY( RANKE, DWORK, 1, E, LDE+1 )
C
      IF( REDA ) THEN
C
C        Set A22.
C
         CALL DLASET( 'Full', LA22, NA22, ZERO, ZERO, A(IR1,IR1), LDA )
         CALL DCOPY( RNKA22, DWORK(IR1), 1, A(IR1,IR1), LDA+1 )
      END IF
C
C     Transpose Z.
C
      DO 170 I = 2, N
         CALL DSWAP( I-1, Z(1,I), 1, Z(I,1), LDZ )
  170 CONTINUE
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TG01ED ***
      END
