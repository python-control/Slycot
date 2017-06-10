      SUBROUTINE TB01ND( JOBU, UPLO, N, P, A, LDA, C, LDC, U, LDU,
     $                   DWORK, INFO )
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
C     To reduce the pair (A,C) to lower or upper observer Hessenberg
C     form using (and optionally accumulating) unitary state-space
C     transformations.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBU    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix U the unitary state-space transformations for
C             reducing the system, as follows:
C             = 'N':  Do not form U;
C             = 'I':  U is initialized to the unit matrix and the
C                     unitary transformation matrix U is returned;
C             = 'U':  The given matrix U is updated by the unitary
C                     transformations used in the reduction.
C
C     UPLO    CHARACTER*1
C             Indicates whether the user wishes the pair (A,C) to be
C             reduced to upper or lower observer Hessenberg form as
C             follows:
C             = 'U':  Upper observer Hessenberg form;
C             = 'L':  Lower observer Hessenberg form.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e. the order of the
C             matrix A.  N >= 0.
C
C     P       (input) INTEGER
C             The actual output dimension, i.e. the number of rows of
C             the matrix C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state transition matrix A to be transformed.
C             On exit, the leading N-by-N part of this array contains
C             the transformed state transition matrix U' * A * U.
C             The annihilated elements are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C to be transformed.
C             On exit, the leading P-by-N part of this array contains
C             the transformed output matrix C * U.
C             The annihilated elements are set to zero.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*)
C             On entry, if JOBU = 'U', then the leading N-by-N part of
C             this array must contain a given matrix U (e.g. from a
C             previous call to another SLICOT routine), and on exit, the
C             leading N-by-N part of this array contains the product of
C             the input matrix U and the state-space transformation
C             matrix which reduces the given pair to observer Hessenberg
C             form.
C             On exit, if JOBU = 'I', then the leading N-by-N part of
C             this array contains the matrix of accumulated unitary
C             similarity transformations which reduces the given pair
C             to observer Hessenberg form.
C             If JOBU = 'N', the array U is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDU = 1 and
C             declare this array to be U(1,1) in the calling program).
C
C     LDU     INTEGER
C             The leading dimension of array U. If JOBU = 'U' or
C             JOBU = 'I', LDU >= MAX(1,N); if JOBU = 'N', LDU >= 1.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (MAX(N,P-1))
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
C     The routine computes a unitary state-space transformation U, which
C     reduces the pair (A,C) to one of the following observer Hessenberg
C     forms:
C
C                                N
C                       |*  . . . . . .  *|
C                       |.               .|
C                       |.               .|
C                       |.               .| N
C                       |*               .|
C            |U'AU|     |   .            .|
C            |----|  =  |     .          .|
C            |CU  |     |       * . . .  *|
C                       -------------------
C                       |         * . .  *|
C                       |           .    .| P
C                       |             .  .|
C                       |                *|
C
C         if UPLO = 'U', or
C
C                               N
C                      |*                |
C                      |.  .             |
C                      |.    .           | P
C                      |*  . . *         |
C            |CU  |    -------------------
C            |----|  = |*  . . . *       |
C            |U'AU|    |.          .     |
C                      |.            .   |
C                      |.               *|
C                      |.               .| N
C                      |.               .|
C                      |.               .|
C                      |*  . . . . . .  *|
C
C     if UPLO = 'L'.
C
C     If P >= N, then the matrix CU is trapezoidal and U'AU is full.
C
C     REFERENCES
C
C     [1] Van Dooren, P. and Verhaegen, M.H.G.
C         On the use of unitary state-space transformations.
C         In : Contemporary Mathematics on Linear Algebra and its Role
C         in Systems Theory, 47, AMS, Providence, 1985.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O((N + P) x N**2) operations and is
C     backward stable (see [1]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TB01BD by M. Vanbegin, and
C     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium.
C
C     REVISIONS
C
C     February 1997.
C
C     KEYWORDS
C
C     Controllability, observer Hessenberg form, orthogonal
C     transformation, unitary transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDC, LDU, N, P
      CHARACTER         JOBU, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*), DWORK(*), U(LDU,*)
C     .. Local Scalars ..
      LOGICAL           LJOBA, LJOBI, LUPLO
      INTEGER           II, J, N1, NJ, P1, PAR1, PAR2, PAR3, PAR4, PAR5,
     $                  PAR6
      DOUBLE PRECISION  DZ
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLARFG, DLASET, SLCT_DLATZM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LUPLO = LSAME( UPLO, 'U' )
      LJOBI = LSAME( JOBU, 'I' )
      LJOBA = LJOBI.OR.LSAME( JOBU, 'U' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBA .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LUPLO .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( .NOT.LJOBA .AND. LDU.LT.1 .OR.
     $              LJOBA .AND. LDU.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return
C
         CALL XERBLA( 'TB01ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. P.EQ.0 )
     $   RETURN
C
      P1 = P + 1
      N1 = N - 1
C
      IF ( LJOBI ) THEN
C
C        Initialize U to the identity matrix.
C
         CALL DLASET( 'Full', N, N, ZERO, ONE, U, LDU )
      END IF
C
C     Perform transformations involving both C and A.
C
      DO 20 J = 1, MIN( P, N1 )
         NJ = N - J
         IF ( LUPLO ) THEN
            PAR1 = P - J + 1
            PAR2 = NJ + 1
            PAR3 = 1
            PAR4 = P - J
            PAR5 = NJ
         ELSE
            PAR1 = J
            PAR2 = J
            PAR3 = J + 1
            PAR4 = P
            PAR5 = N
         END IF
C
         CALL DLARFG( NJ+1, C(PAR1,PAR2), C(PAR1,PAR3), LDC, DZ )
C
C        Update A.
C
         CALL SLCT_DLATZM( 'Left', NJ+1, N, C(PAR1,PAR3), LDC, DZ,
     $                A(PAR2,1), A(PAR3,1), LDA, DWORK )
         CALL SLCT_DLATZM( 'Right', N, NJ+1, C(PAR1,PAR3), LDC, DZ,
     $                A(1,PAR2), A(1,PAR3), LDA, DWORK )
C
         IF ( LJOBA ) THEN
C
C           Update U.
C
            CALL SLCT_DLATZM( 'Right', N, NJ+1, C(PAR1,PAR3), LDC, DZ,
     $                   U(1,PAR2), U(1,PAR3), LDU, DWORK )
         END IF
C
         IF ( J.NE.P ) THEN
C
C           Update C.
C
            CALL SLCT_DLATZM( 'Right', PAR4-PAR3+1, NJ+1,
     $                   C(PAR1,PAR3), LDC,
     $                   DZ, C(PAR3,PAR2), C(PAR3,PAR3), LDC, DWORK )
         END IF
C
         DO 10 II = PAR3, PAR5
            C(PAR1,II) = ZERO
   10    CONTINUE
C
   20 CONTINUE
C
      DO 40 J = P1, N1
C
C        Perform next transformations only involving A.
C
         NJ = N - J
         IF ( LUPLO ) THEN
            PAR1 = N + P1 - J
            PAR2 = NJ + 1
            PAR3 = 1
            PAR4 = NJ
            PAR5 = 1
            PAR6 = N + P - J
         ELSE
            PAR1 = J - P
            PAR2 = J
            PAR3 = J + 1
            PAR4 = N
            PAR5 = J - P + 1
            PAR6 = N
         END IF
C
         IF ( NJ.GT.0 ) THEN
C
            CALL DLARFG( NJ+1, A(PAR1,PAR2), A(PAR1,PAR3), LDA, DZ )
C
C           Update A.
C
            CALL SLCT_DLATZM( 'Left', NJ+1, N, A(PAR1,PAR3), LDA, DZ,
     $                   A(PAR2,1), A(PAR3,1), LDA, DWORK )
            CALL SLCT_DLATZM( 'Right', PAR6-PAR5+1, NJ+1, A(PAR1,PAR3),
     $                    LDA,
     $                    DZ, A(PAR5,PAR2), A(PAR5,PAR3), LDA, DWORK )
C
            IF ( LJOBA ) THEN
C
C              Update U.
C
               CALL SLCT_DLATZM( 'Right', N, NJ+1, A(PAR1,PAR3),
     $                      LDA, DZ,
     $                      U(1,PAR2), U(1,PAR3), LDU, DWORK )
            END IF
C
            DO 30 II = PAR3, PAR4
               A(PAR1,II) = ZERO
   30       CONTINUE
C
         END IF
C
   40 CONTINUE
C
      RETURN
C *** Last line of TB01ND ***
      END
