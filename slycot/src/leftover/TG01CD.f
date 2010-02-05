      SUBROUTINE TG01CD( COMPQ, L, N, M, A, LDA, E, LDE, B, LDB, Q, LDQ,
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
C     To reduce the descriptor system pair (A-lambda E,B) to the
C     QR-coordinate form by computing an orthogonal transformation
C     matrix Q such that the transformed descriptor system pair
C     (Q'*A-lambda Q'*E, Q'*B) has the descriptor matrix Q'*E
C     in an upper trapezoidal form.
C     The left orthogonal transformations performed to reduce E
C     can be optionally accumulated.
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
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of rows of matrices A, B, and E.  L >= 0.
C
C     N       (input) INTEGER
C             The number of columns of matrices A and E.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of matrix B.  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix Q'*A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the descriptor matrix E.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix Q'*E in upper trapezoidal form,
C             i.e.
C
C                      ( E11 )
C               Q'*E = (     ) ,     if L >= N ,
C                      (  0  )
C             or
C
C               Q'*E = ( E11 E12 ),  if L < N ,
C
C             where E11 is an MIN(L,N)-by-MIN(L,N) upper triangular
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
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, MIN(L,N) + MAX(L,N,M)).
C             For optimum performance
C             LWORK >= MAX(1, MIN(L,N) + MAX(L,N,M)*NB),
C             where NB is the optimal blocksize.
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
C     The routine computes the QR factorization of E to reduce it
C     to the upper trapezoidal form.
C
C     The transformations are also applied to the rest of system
C     matrices
C
C         A <- Q' * A ,  B <- Q' * B.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable and requires
C     0( L*L*N )  floating point operations.
C
C     CONTRIBUTOR
C
C     C. Oara, University "Politehnica" Bucharest.
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the RASP routine RPDSQR.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 1999,
C     May 2003.
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
      CHARACTER          COMPQ
      INTEGER            INFO, L, LDA, LDB, LDE, LDQ, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), DWORK( * ),
     $                   E( LDE, * ), Q( LDQ, * )
C     .. Local Scalars ..
      LOGICAL            ILQ
      INTEGER            ICOMPQ, LN, WRKOPT
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLASET, DORMQR, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
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
C     Test the input parameters.
C
      INFO = 0
      WRKOPT = MAX( 1, MIN( L, N ) + MAX( L, N, M ) )
      IF( ICOMPQ.EQ.0 ) THEN
         INFO = -1
      ELSE IF( L.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -6
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.L ) ) THEN
         INFO = -10
      ELSE IF( ( ILQ .AND. LDQ.LT.L ) .OR. LDQ.LT.1 ) THEN
         INFO = -12
      ELSE IF( LDWORK.LT.WRKOPT ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01CD', -INFO )
         RETURN
      END IF
C
C     Initialize Q if necessary.
C
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', L, L, ZERO, ONE, Q, LDQ )
C
C     Quick return if possible.
C
      IF( L.EQ.0 .OR. N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      LN = MIN( L, N )
C
C     Compute the QR decomposition of E.
C
C     Workspace: need   MIN(L,N) + N;
C                prefer MIN(L,N) + N*NB.
C
      CALL DGEQRF( L, N, E, LDE, DWORK, DWORK( LN+1 ), LDWORK-LN, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
C
C     Apply transformation on the rest of matrices.
C
C     A <-- Q' * A.
C     Workspace: need   MIN(L,N) + N;
C                prefer MIN(L,N) + N*NB.
C
      CALL DORMQR( 'Left', 'Transpose', L, N, LN, E, LDE, DWORK,
     $             A, LDA, DWORK( LN+1 ), LDWORK-LN, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
C
C     B <-- Q' * B.
C     Workspace: need   MIN(L,N) + M;
C                prefer MIN(L,N) + M*NB.
C
      IF ( M.GT.0 ) THEN
         CALL DORMQR( 'Left', 'Transpose', L, M, LN, E, LDE, DWORK,
     $                B, LDB, DWORK( LN+1 ), LDWORK-LN, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
      END IF
C
C     Q <-- Q1 * Q.
C     Workspace: need   MIN(L,N) + L;
C                prefer MIN(L,N) + L*NB.
C
      IF( ILQ ) THEN
         CALL DORMQR( 'Right', 'No Transpose', L, L, LN, E, LDE, DWORK,
     $                Q, LDQ, DWORK( LN+1 ), LDWORK-LN, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
      END IF
C
C     Set lower triangle of E to zero.
C
      IF( L.GE.2 )
     $   CALL DLASET( 'Lower', L-1, LN, ZERO, ZERO, E( 2, 1 ), LDE )
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TG01CD ***
      END
