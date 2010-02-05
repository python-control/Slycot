      SUBROUTINE TG01DD( COMPZ, L, N, P, A, LDA, E, LDE, C, LDC, Z, LDZ,
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
C     To reduce the descriptor system pair (C,A-lambda E) to the
C     RQ-coordinate form by computing an orthogonal transformation
C     matrix Z such that the transformed descriptor system pair
C     (C*Z,A*Z-lambda E*Z) has the descriptor matrix E*Z in an upper
C     trapezoidal form.
C     The right orthogonal transformations performed to reduce E can
C     be optionally accumulated.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPZ   CHARACTER*1
C             = 'N':  do not compute Z;
C             = 'I':  Z is initialized to the unit matrix, and the
C                     orthogonal matrix Z is returned;
C             = 'U':  Z must contain an orthogonal matrix Z1 on entry,
C                     and the product Z1*Z is returned.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of rows of matrices A and E.  L >= 0.
C
C     N       (input) INTEGER
C             The number of columns of matrices A, E, and C.  N >= 0.
C
C     P       (input) INTEGER
C             The number of rows of matrix C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix A*Z.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the descriptor matrix E.
C             On exit, the leading L-by-N part of this array contains
C             the transformed matrix E*Z in upper trapezoidal form,
C             i.e.
C
C                      ( E11 )
C                E*Z = (     ) ,  if L >= N ,
C                      (  R  )
C             or
C
C                E*Z = ( 0  R ),  if L < N ,
C
C             where R is an MIN(L,N)-by-MIN(L,N) upper triangular
C             matrix.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,L).
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
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, MIN(L,N) + MAX(L,N,P)).
C             For optimum performance
C             LWORK >= MAX(1, MIN(L,N) + MAX(L,N,P)*NB),
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
C     The routine computes the RQ factorization of E to reduce it
C     the upper trapezoidal form.
C
C     The transformations are also applied to the rest of system
C     matrices
C
C         A <- A * Z,  C <- C * Z.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable and requires
C     0( L*N*N )  floating point operations.
C
C     CONTRIBUTOR
C
C     C. Oara, University "Politehnica" Bucharest.
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the RASP routine RPDSRQ.
C
C     REVISIONS
C
C     July 1999, V. Sima, Research Institute for Informatics, Bucharest.
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
      CHARACTER          COMPZ
      INTEGER            INFO, L, LDA, LDC, LDE, LDWORK, LDZ, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), DWORK( * ),
     $                   E( LDE, * ), Z( LDZ, * )
C     .. Local Scalars ..
      LOGICAL            ILZ
      INTEGER            ICOMPZ, LN, WRKOPT
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DGERQF, DLASET, DORMRQ, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
C
C     .. Executable Statements ..
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
C
C     Test the input parameters.
C
      INFO = 0
      WRKOPT = MAX( 1, MIN( L, N ) + MAX( L, N, P ) )
      IF( ICOMPZ.EQ.0 ) THEN
         INFO = -1
      ELSE IF( L.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -6
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -12
      ELSE IF( LDWORK.LT.WRKOPT ) THEN
         INFO = -14
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'TG01DD', -INFO )
         RETURN
      END IF
C
C     Initialize Q if necessary.
C
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Quick return if possible.
C
      IF( L.EQ.0 .OR. N.EQ.0 ) THEN
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
      LN = MIN( L, N )
C
C     Compute the RQ decomposition of E, E = R*Z.
C
C     Workspace: need   MIN(L,N) + L;
C                prefer MIN(L,N) + L*NB.
C
      CALL DGERQF( L, N, E, LDE, DWORK, DWORK( LN+1 ), LDWORK-LN, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
C
C     Apply transformation on the rest of matrices.
C
C     A <--  A * Z'.
C     Workspace: need   MIN(L,N) + L;
C                prefer MIN(L,N) + L*NB.
C
      CALL DORMRQ( 'Right', 'Transpose', L, N, LN, E( L-LN+1,1 ), LDE,
     $              DWORK, A, LDA, DWORK( LN+1 ), LDWORK-LN, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
C
C     C <-- C * Z'.
C     Workspace: need   MIN(L,N) + P;
C                prefer MIN(L,N) + P*NB.
C
      CALL DORMRQ( 'Right', 'Transpose', P, N, LN, E( L-LN+1,1 ), LDE,
     $              DWORK, C, LDC, DWORK( LN+1 ), LDWORK-LN, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
C
C     Z <-- Z1 * Z'.
C     Workspace: need   MIN(L,N) + N;
C                prefer MIN(L,N) + N*NB.
C
      IF( ILZ ) THEN
         CALL DORMRQ( 'Right', 'Transpose', N, N, LN, E( L-LN+1,1 ),
     $                LDE, DWORK, Z, LDZ, DWORK( LN+1 ), LDWORK-LN,
     $                INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK( LN+1 ) ) + LN )
      END IF
C
C     Set lower triangle of E to zero.
C
      IF( L.LT.N ) THEN
         CALL DLASET( 'Full', L, N-L, ZERO, ZERO, E, LDE )
         IF( L.GE.2 ) CALL DLASET( 'Lower', L-1, L, ZERO, ZERO,
     $                             E( 2, N-L+1 ), LDE )
      ELSE
         IF( N.GE.2 ) CALL DLASET( 'Lower', N-1, N, ZERO, ZERO,
     $                             E( L-N+2, 1 ), LDE )
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TG01DD ***
      END
