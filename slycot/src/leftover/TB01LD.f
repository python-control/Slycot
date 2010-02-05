      SUBROUTINE TB01LD( DICO, STDOM, JOBA, N, M, P, ALPHA, A, LDA, B,
     $                   LDB, C, LDC, NDIM, U, LDU, WR, WI, DWORK,
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
C     To reduce the system state matrix A to an ordered upper real
C     Schur form by using an orthogonal similarity transformation
C     A <-- U'*A*U and to apply the transformation to the matrices
C     B and C: B <-- U'*B and C <-- C*U.
C     The leading block of the resulting A has eigenvalues in a
C     suitably defined domain of interest.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     STDOM   CHARACTER*1
C             Specifies whether the domain of interest is of stability
C             type (left part of complex plane or inside of a circle)
C             or of instability type (right part of complex plane or
C             outside of a circle) as follows:
C             = 'S':  stability type domain;
C             = 'U':  instability type domain.
C
C     JOBA    CHARACTER*1
C             Specifies the shape of the state dynamics matrix on entry
C             as follows:
C             = 'S':  A is in an upper real Schur form;
C             = 'G':  A is a general square dense matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state-space representation,
C             i.e. the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs, or of columns of B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs, or of rows of C.  P >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION.
C             Specifies the boundary of the domain of interest for the
C             eigenvalues of A. For a continuous-time system
C             (DICO = 'C'), ALPHA is the boundary value for the real
C             parts of eigenvalues, while for a discrete-time system
C             (DICO = 'D'), ALPHA >= 0 represents the boundary value
C             for the moduli of eigenvalues.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the unreduced state dynamics matrix A.
C             If JOBA = 'S' then A must be a matrix in real Schur form.
C             On exit, the leading N-by-N part of this array contains
C             the ordered real Schur matrix U' * A * U with the elements
C             below the first subdiagonal set to zero.
C             The leading NDIM-by-NDIM part of A has eigenvalues in the
C             domain of interest and the trailing (N-NDIM)-by-(N-NDIM)
C             part has eigenvalues outside the domain of interest.
C             The domain of interest for lambda(A), the eigenvalues
C             of A, is defined by the parameters ALPHA, DICO and STDOM
C             as follows:
C             For a continuous-time system (DICO = 'C'):
C               Real(lambda(A)) < ALPHA if STDOM = 'S';
C               Real(lambda(A)) > ALPHA if STDOM = 'U';
C             For a discrete-time system (DICO = 'D'):
C               Abs(lambda(A)) < ALPHA if STDOM = 'S';
C               Abs(lambda(A)) > ALPHA if STDOM = 'U'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed input matrix U' * B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed output matrix C * U.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     NDIM    (output) INTEGER
C             The number of eigenvalues of A lying inside the domain of
C             interest for eigenvalues.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,N)
C             The leading N-by-N part of this array contains the
C             orthogonal transformation matrix used to reduce A to the
C             real Schur form and/or to reorder the diagonal blocks of
C             real Schur form of A. The first NDIM columns of U form
C             an orthogonal basis for the invariant subspace of A
C             corresponding to the first NDIM eigenvalues.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= max(1,N).
C
C     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N)
C             WR and WI contain the real and imaginary parts,
C             respectively, of the computed eigenvalues of A. The
C             eigenvalues will be in the same order that they appear on
C             the diagonal of the output real Schur form of A. Complex
C             conjugate pairs of eigenvalues will appear consecutively
C             with the eigenvalue having the positive imaginary part
C             first.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.
C             LDWORK >= MAX(1,N)   if JOBA = 'S';
C             LDWORK >= MAX(1,3*N) if JOBA = 'G'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0: successful exit;
C             < 0: if INFO = -i, the i-th argument had an illegal
C                  value;
C             = 1: the QR algorithm failed to compute all the
C                  eigenvalues of A;
C             = 2: a failure occured during the ordering of the real
C                  Schur form of A.
C
C     METHOD
C
C     Matrix A is reduced to an ordered upper real Schur form using an
C     orthogonal similarity transformation A <-- U'*A*U. This
C     transformation is determined so that the leading block of the
C     resulting A has eigenvalues in a suitably defined domain of
C     interest. Then, the transformation is applied to the matrices B
C     and C: B <-- U'*B and C <-- C*U.
C
C     NUMERICAL ASPECTS
C                                     3
C     The algorithm requires about 14N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine SRSFOD.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001.
C
C     KEYWORDS
C
C     Invariant subspace, orthogonal transformation, real Schur form,
C     similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBA, STDOM
      INTEGER          INFO, LDA, LDB, LDC, LDU, LDWORK, M, N, NDIM, P
      DOUBLE PRECISION ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), U(LDU,*),
     $                 WI(*), WR(*)
C     .. Local Scalars ..
      LOGICAL          DISCR, LJOBG
      INTEGER          I, IERR, LDWP, SDIM
      DOUBLE PRECISION WRKOPT
C     .. Local Arrays ..
      LOGICAL          BWORK( 1 )
C     .. External Functions ..
      LOGICAL          LSAME, SELECT
      EXTERNAL         LSAME, SELECT
C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGEES, DGEMM, DGEMV, DLACPY, DLASET,
     $                 MB03QD, MB03QX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        DBLE, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
      DISCR = LSAME( DICO, 'D' )
      LJOBG = LSAME( JOBA, 'G' )
C
C     Check input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( STDOM, 'S' ) .OR.
     $                 LSAME( STDOM, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( JOBA, 'S' ) .OR. LJOBG ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( DISCR .AND. ALPHA.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF( LDWORK.LT.MAX( 1,   N ) .OR.
     $         LDWORK.LT.MAX( 1, 3*N ) .AND. LJOBG ) THEN
         INFO = -20
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01LD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      NDIM = 0
      IF( N.EQ.0 )
     $   RETURN
C
      IF( LSAME( JOBA, 'G' ) ) THEN
C
C        Reduce A to real Schur form using an orthogonal similarity
C        transformation A <- U'*A*U, accumulate the transformation in U
C        and compute the eigenvalues of A in (WR,WI).
C
C        Workspace:  need   3*N;
C                    prefer larger.
C
         CALL DGEES( 'Vectors', 'Not ordered', SELECT, N, A, LDA, SDIM,
     $               WR, WI, U, LDU, DWORK, LDWORK, BWORK, INFO )
         WRKOPT = DWORK( 1 )
         IF( INFO.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
      ELSE
C
C        Initialize U with an identity matrix.
C
         CALL DLASET( 'Full', N, N, ZERO, ONE, U, LDU )
         WRKOPT = 0
      END IF
C
C     Separate the spectrum of A. The leading NDIM-by-NDIM submatrix of
C     A corresponds to the eigenvalues of interest.
C     Workspace:  need   N.
C
      CALL MB03QD( DICO, STDOM, 'Update', N, 1, N, ALPHA, A, LDA,
     $             U, LDU, NDIM, DWORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
C
C     Compute the eigenvalues.
C
      CALL MB03QX( N, A, LDA, WR, WI, IERR )
C
C     Apply the transformation: B <-- U'*B.
C
      IF( LDWORK.LT.N*M ) THEN
C
C        Not enough working space for using DGEMM.
C
         DO 10 I = 1, M
            CALL DCOPY( N, B(1,I), 1, DWORK, 1 )
            CALL DGEMV( 'Transpose', N, N, ONE, U, LDU, DWORK, 1, ZERO,
     $                  B(1,I), 1 )
   10    CONTINUE
C
      ELSE
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
         CALL DGEMM( 'Transpose', 'No transpose', N, M, N, ONE, U, LDU,
     $               DWORK, N, ZERO, B, LDB )
         WRKOPT = MAX( WRKOPT, DBLE( N*M ) )
      END IF
C
C     Apply the transformation: C <-- C*U.
C
      IF( LDWORK.LT.N*P ) THEN
C
C        Not enough working space for using DGEMM.
C
         DO 20 I = 1, P
            CALL DCOPY( N, C(I,1), LDC, DWORK, 1 )
            CALL DGEMV( 'Transpose', N, N, ONE, U, LDU, DWORK, 1, ZERO,
     $                  C(I,1), LDC )
   20    CONTINUE
C
      ELSE
         LDWP = MAX( 1, P )
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK, LDWP )
         CALL DGEMM( 'No transpose', 'No transpose', P, N, N, ONE,
     $               DWORK, LDWP, U, LDU, ZERO, C, LDC )
         WRKOPT = MAX( WRKOPT, DBLE( N*P ) )
      END IF
C
      DWORK( 1 ) = WRKOPT
C
      RETURN
C *** Last line of TB01LD ***
      END
