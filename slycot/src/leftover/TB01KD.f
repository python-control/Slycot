      SUBROUTINE TB01KD( DICO, STDOM, JOBA, N, M, P, ALPHA, A, LDA, B,
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
C     To compute an additive spectral decomposition of the transfer-
C     function matrix of the system (A,B,C) by reducing the system
C     state-matrix A to a block-diagonal form.
C     The system matrices are transformed as
C     A <-- inv(U)*A*U, B <--inv(U)*B and C <-- C*U.
C     The leading diagonal block of the resulting A has eigenvalues
C     in a suitably defined domain of interest.
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
C             (DICO = 'D'), ALPHA >= 0 represents the boundary value for
C             the moduli of eigenvalues.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the unreduced state dynamics matrix A.
C             If JOBA = 'S' then A must be a matrix in real Schur form.
C             On exit, the leading N-by-N part of this array contains a
C             block diagonal matrix inv(U) * A * U with two diagonal
C             blocks in real Schur form with the elements below the
C             first subdiagonal set to zero.
C             The leading NDIM-by-NDIM block of A has eigenvalues in the
C             domain of interest and the trailing (N-NDIM)-by-(N-NDIM)
C             block has eigenvalues outside the domain of interest.
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
C             the transformed input matrix inv(U) * B.
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
C             transformation matrix used to reduce A to the block-
C             diagonal form. The first NDIM columns of U span the
C             invariant subspace of A corresponding to the eigenvalues
C             of its leading diagonal block. The last N-NDIM columns
C             of U span the reducing subspace of A corresponding to
C             the eigenvalues of the trailing diagonal block of A.
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
C                  Schur form of A;
C             = 3: the separation of the two diagonal blocks failed
C                  because of very close eigenvalues.
C
C     METHOD
C
C     A similarity transformation U is determined that reduces the
C     system state-matrix A to a block-diagonal form (with two diagonal
C     blocks), so that the leading diagonal block of the resulting A has
C     eigenvalues in a specified domain of the complex plane. The
C     determined transformation is applied to the system (A,B,C) as
C       A <-- inv(U)*A*U, B <-- inv(U)*B and C <-- C*U.
C
C     REFERENCES
C
C     [1] Safonov, M.G., Jonckheere, E.A., Verma, M., Limebeer, D.J.N.
C         Synthesis of positive real multivariable feedback systems.
C         Int. J. Control, pp. 817-842, 1987.
C
C     NUMERICAL ASPECTS
C                                     3
C     The algorithm requires about 14N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine SADSDC.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Invariant subspace, real Schur form, similarity transformation,
C     spectral factorization.
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
      INTEGER          NDIM1, NR
      DOUBLE PRECISION SCALE
C     .. External Functions ..
      LOGICAL          LSAME
      EXTERNAL         LSAME
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DLASET, DTRSYL, TB01LD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        MAX
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
         CALL XERBLA( 'TB01KD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      NDIM = 0
      IF( N.EQ.0 )
     $   RETURN
C
C     Reduce A to an ordered real Schur form using an orthogonal
C     similarity transformation A <- U'*A*U and accumulate the
C     transformations in U. The reordering of the real Schur form of A
C     is performed in accordance with the values of the parameters DICO,
C     STDOM and ALPHA. Apply the transformation to B and C: B <- U'*B
C     and C <- C*U. The eigenvalues of A are computed in (WR,WI).
C
C     Workspace:  need   3*N (if JOBA = 'G'), or N (if JOBA = 'S');
C                 prefer larger.
C
      CALL TB01LD( DICO, STDOM, JOBA, N, M, P, ALPHA, A, LDA, B, LDB, C,
     $             LDC, NDIM, U, LDU, WR, WI, DWORK, LDWORK, INFO )
C
      IF ( INFO.NE.0 )
     $    RETURN
C
      IF ( NDIM.GT.0 .AND. NDIM.LT.N ) THEN
C
C        Reduce A to a block-diagonal form by a similarity
C        transformation of the form
C               -1                  ( I -X )
C         A <- T  AT,  where    T = (      )  and X satisfies the
C                                   ( 0  I )
C        Sylvester equation
C
C          A11*X - X*A22 = A12.
C
         NR = N - NDIM
         NDIM1 = NDIM + 1
         CALL DTRSYL( 'N', 'N', -1, NDIM, NR, A, LDA, A(NDIM1,NDIM1),
     $                LDA, A(1,NDIM1), LDA, SCALE, INFO )
         IF ( INFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
C                      -1
C        Compute B <- T  B,  C <- CT,  U <- UT.
C
         SCALE = ONE/SCALE
         CALL DGEMM( 'N', 'N', NDIM, M, NR, SCALE, A(1,NDIM1), LDA,
     $               B(NDIM1,1), LDB, ONE, B, LDB )
         CALL DGEMM( 'N', 'N', P, NR, NDIM, -SCALE, C, LDC, A(1,NDIM1),
     $               LDA, ONE, C(1,NDIM1), LDC )
         CALL DGEMM( 'N', 'N', N, NR, NDIM, -SCALE, U, LDU, A(1,NDIM1),
     $               LDA, ONE, U(1,NDIM1), LDU )
C
C        Set A12 to zero.
C
         CALL DLASET( 'Full', NDIM, NR, ZERO, ZERO, A(1,NDIM1), LDA )
      END IF
C
C     Set to zero the lower triangular part under the first subdiagonal
C     of A.
C
      IF ( N.GT.2 )
     $   CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, A( 3, 1 ), LDA )
      RETURN
C *** Last line of TB01KD ***
      END
