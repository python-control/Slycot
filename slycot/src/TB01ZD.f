      SUBROUTINE TB01ZD( JOBZ, N, P, A, LDA, B, C, LDC, NCONT, Z, LDZ,
     $                   TAU, TOL, DWORK, LDWORK, INFO )
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
C     To find a controllable realization for the linear time-invariant
C     single-input system
C
C             dX/dt = A * X + B * U,
C                Y  = C * X,
C
C     where A is an N-by-N matrix, B is an N element vector, C is an
C     P-by-N matrix, and A and B are reduced by this routine to
C     orthogonal canonical form using (and optionally accumulating)
C     orthogonal similarity transformations, which are also applied
C     to C.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBZ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the orthogonal similarity transformations for
C             reducing the system, as follows:
C             = 'N':  Do not form Z and do not store the orthogonal
C                     transformations;
C             = 'F':  Do not form Z, but store the orthogonal
C                     transformations in the factored form;
C             = 'I':  Z is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z is returned.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e. the order of the matrix A.  N >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs, or of rows of C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading NCONT-by-NCONT upper Hessenberg
C             part of this array contains the canonical form of the
C             state dynamics matrix, given by Z' * A * Z, of a
C             controllable realization for the original system. The
C             elements below the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, the original input/state vector B.
C             On exit, the leading NCONT elements of this array contain
C             canonical form of the input/state vector, given by Z' * B,
C             with all elements but B(1) set to zero.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output/state matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed output/state matrix, given by C * Z, and
C             the leading P-by-NCONT part contains the output/state
C             matrix of the controllable realization.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     NCONT   (output) INTEGER
C             The order of the controllable state-space representation.
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If JOBZ = 'I', then the leading N-by-N part of this array
C             contains the matrix of accumulated orthogonal similarity
C             transformations which reduces the given system to
C             orthogonal canonical form.
C             If JOBZ = 'F', the elements below the diagonal, with the
C             array TAU, represent the orthogonal transformation matrix
C             as a product of elementary reflectors. The transformation
C             matrix can then be obtained by calling the LAPACK Library
C             routine DORGQR.
C             If JOBZ = 'N', the array Z is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDZ = 1 and
C             declare this array to be Z(1,1) in the calling program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If JOBZ = 'I' or
C             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The elements of TAU contain the scalar factors of the
C             elementary reflectors used in the reduction of B and A.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the
C             controllability of (A,B). If the user sets TOL > 0, then
C             the given value of TOL is used as an absolute tolerance;
C             elements with absolute value less than TOL are considered
C             neglijible. If the user sets TOL <= 0, then an implicitly
C             computed, default tolerance, defined by
C             TOLDEF = N*EPS*MAX( NORM(A), NORM(B) ) is used instead,
C             where EPS is the machine precision (see LAPACK Library
C             routine DLAMCH).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK. LDWORK >= MAX(1,N,P).
C             For optimum performance LDWORK should be larger.
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
C     The Householder matrix which reduces all but the first element
C     of vector B to zero is found and this orthogonal similarity
C     transformation is applied to the matrix A. The resulting A is then
C     reduced to upper Hessenberg form by a sequence of Householder
C     transformations. Finally, the order of the controllable state-
C     space representation (NCONT) is determined by finding the position
C     of the first sub-diagonal element of A which is below an
C     appropriate zero threshold, either TOL or TOLDEF (see parameter
C     TOL); if NORM(B) is smaller than this threshold, NCONT is set to
C     zero, and no computations for reducing the system to orthogonal
C     canonical form are performed.
C     All orthogonal transformations determined in this process are also
C     applied to the matrix C, from the right.
C
C     REFERENCES
C
C     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D.
C         Orthogonal Invariants and Canonical Forms for Linear
C         Controllable Systems.
C         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981.
C
C     [2] Hammarling, S.J.
C         Notes on the use of orthogonal similarity transformations in
C         control.
C         NPL Report DITC 8/82, August 1982.
C
C     [3] Paige, C.C
C         Properties of numerical algorithms related to computing
C         controllability.
C         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001,
C     Sept. 2003.
C
C     KEYWORDS
C
C     Controllability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBZ
      INTEGER           INFO, LDA, LDC, LDWORK, LDZ, N, NCONT, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(*), C(LDC,*), DWORK(*), TAU(*),
     $                  Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           LJOBF, LJOBI, LJOBZ
      INTEGER           ITAU, J
      DOUBLE PRECISION  ANORM, B1, BNORM, FANORM, FBNORM, H, THRESH,
     $                  TOLDEF, WRKOPT
C     .. Local Arrays ..
      DOUBLE PRECISION  NBLK(1)
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEHRD, DLACPY, DLARF, DLARFG, DLASET, DORGQR,
     $                  DORMHR, MB01PD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX
C     .. Executable Statements ..
C
      INFO = 0
      LJOBF = LSAME( JOBZ, 'F' )
      LJOBI = LSAME( JOBZ, 'I' )
      LJOBZ = LJOBF.OR.LJOBI
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBZ .AND. .NOT.LSAME( JOBZ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( LDZ.LT.1 .OR. ( LJOBZ .AND. LDZ.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDWORK.LT.MAX( 1, N, P ) ) THEN
         INFO = -15
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01ZD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      NCONT = 0
      DWORK(1) = ONE
      IF ( N.EQ.0 )
     $   RETURN
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      WRKOPT = ONE
C
C     Calculate the absolute norms of A and B (used for scaling).
C
      ANORM = DLANGE( 'Max', N, N, A, LDA, DWORK )
      BNORM = DLANGE( 'Max', N, 1, B, N, DWORK )
C
C     Return if matrix B is zero.
C
      IF( BNORM.EQ.ZERO ) THEN
         IF( LJOBF ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ZERO, Z, LDZ )
            CALL DLASET( 'Full', N, 1, ZERO, ZERO, TAU, N )
         ELSE IF( LJOBI ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         END IF
         RETURN
      END IF
C
C     Scale (if needed) the matrices A and B.
C
      CALL MB01PD( 'S', 'G', N, N, 0, 0, ANORM, 0, NBLK, A, LDA, INFO )
      CALL MB01PD( 'S', 'G', N, 1, 0, 0, BNORM, 0, NBLK, B, N, INFO )
C
C     Calculate the Frobenius norm of A and the 1-norm of B (used for
C     controlability test).
C
      FANORM = DLANGE( 'Frobenius', N, N, A, LDA, DWORK )
      FBNORM = DLANGE( '1-norm', N, 1, B, N, DWORK )
C
      TOLDEF = TOL
      IF ( TOLDEF.LE.ZERO ) THEN
C
C        Use the default tolerance in controllability determination.
C
         THRESH = DBLE(N)*DLAMCH( 'EPSILON' )
         TOLDEF = THRESH*MAX( FANORM, FBNORM )
      END IF
C
      ITAU = 1
      IF ( FBNORM.GT.TOLDEF ) THEN
C
C        B is not negligible compared with A.
C
         IF ( N.GT.1 ) THEN
C
C           Transform B by a Householder matrix Z1: store vector
C           describing this temporarily in B and in the local scalar H.
C
            CALL DLARFG( N, B(1), B(2), 1, H )
C
            B1 = B(1)
            B(1) = ONE
C
C           Form Z1 * A * Z1.
C           Workspace: need N.
C
            CALL DLARF( 'Right', N, N, B, 1, H, A, LDA, DWORK )
            CALL DLARF( 'Left',  N, N, B, 1, H, A, LDA, DWORK )
C
C           Form C * Z1.
C           Workspace: need P.
C
            CALL DLARF( 'Right', P, N, B, 1, H, C, LDC, DWORK )
C
            B(1) = B1
            TAU(1) = H
            ITAU = ITAU + 1
         ELSE
            B1 = B(1)
            TAU(1) = ZERO
         END IF
C
C        Reduce modified A to upper Hessenberg form by an orthogonal
C        similarity transformation with matrix Z2.
C        Workspace: need N;  prefer N*NB.
C
         CALL DGEHRD( N, 1, N, A, LDA, TAU(ITAU), DWORK, LDWORK, INFO )
         WRKOPT = DWORK(1)
C
C        Form C * Z2.
C        Workspace: need P;  prefer P*NB.
C
         CALL DORMHR( 'Right', 'No transpose', P, N, 1, N, A, LDA,
     $                TAU(ITAU), C, LDC, DWORK, LDWORK, INFO )
         WRKOPT = MAX( WRKOPT, DWORK(1) )
C
         IF ( LJOBZ ) THEN
C
C           Save the orthogonal transformations used, so that they could
C           be accumulated by calling DORGQR routine.
C
            IF ( N.GT.1 )
     $         CALL DLACPY( 'Full',  N-1, 1, B(2), N-1, Z(2,1), LDZ )
            IF ( N.GT.2 )
     $         CALL DLACPY( 'Lower', N-2, N-2, A(3,1), LDA, Z(3,2),
     $                      LDZ )
            IF ( LJOBI ) THEN
C
C              Form the orthogonal transformation matrix Z = Z1 * Z2.
C              Workspace: need N;  prefer N*NB.
C
               CALL DORGQR( N, N, N, Z, LDZ, TAU, DWORK, LDWORK, INFO )
               WRKOPT = MAX( WRKOPT, DWORK(1) )
            END IF
         END IF
C
C        Annihilate the lower part of A and B.
C
         IF ( N.GT.2 )
     $      CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
         IF ( N.GT.1 )
     $      CALL DLASET( 'Full',  N-1, 1, ZERO, ZERO, B(2), N-1 )
C
C        Find NCONT by checking sizes of the sub-diagonal elements of
C        transformed A.
C
         IF ( TOL.LE.ZERO )
     $      TOLDEF = THRESH*MAX( FANORM, ABS( B1 ) )
C
         J = 1
C
C        WHILE ( J < N and ABS( A(J+1,J) ) > TOLDEF ) DO
C
   10    CONTINUE
         IF ( J.LT.N ) THEN
            IF ( ABS( A(J+1,J) ).GT.TOLDEF ) THEN
               J = J + 1
               GO TO 10
            END IF
         END IF
C
C        END WHILE 10
C
C        First negligible sub-diagonal element found, if any: set NCONT.
C
         NCONT = J
         IF ( J.LT.N )
     $      A(J+1,J) = ZERO
C
C        Undo scaling of A and B.
C
         CALL MB01PD( 'U', 'H', NCONT, NCONT, 0, 0, ANORM, 0, NBLK, A,
     $                LDA, INFO )
         CALL MB01PD( 'U', 'G', 1, 1, 0, 0, BNORM, 0, NBLK, B, N, INFO )
         IF ( NCONT.LT.N )
     $      CALL MB01PD( 'U', 'G', N, N-NCONT, 0, 0, ANORM, 0, NBLK,
     $                   A(1,NCONT+1), LDA, INFO )
      ELSE
C
C        B is negligible compared with A. No computations for reducing
C        the system to orthogonal canonical form have been performed,
C        except scaling (which is undoed).
C
         CALL MB01PD( 'U', 'G', N, N, 0, 0, ANORM, 0, NBLK, A, LDA,
     $                INFO )
         CALL MB01PD( 'U', 'G', N, 1, 0, 0, BNORM, 0, NBLK, B, N, INFO )
         IF( LJOBF ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ZERO, Z, LDZ )
            CALL DLASET( 'Full', N, 1, ZERO, ZERO, TAU, N )
         ELSE IF( LJOBI ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         END IF
      END IF
C
C     Set optimal workspace dimension.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TB01ZD ***
      END
