      SUBROUTINE MB05OY( JOB, N, LOW, IGH, A, LDA, SCALE, INFO )
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
C     To restore a matrix after it has been transformed by applying
C     balancing transformations (permutations and scalings), as
C     determined by LAPACK Library routine DGEBAL.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the type of backward transformation required,
C             as follows:
C             = 'N', do nothing, return immediately;
C             = 'P', do backward transformation for permutation only;
C             = 'S', do backward transformation for scaling only;
C             = 'B', do backward transformations for both permutation
C                    and scaling.
C             JOB must be the same as the argument JOB supplied
C             to DGEBAL.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     LOW     (input) INTEGER
C     IGH     (input) INTEGER
C             The integers LOW and IGH determined by DGEBAL.
C             1 <= LOW <= IGH <= N, if N > 0; LOW=1 and IGH=0, if N=0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix to be back-transformed.
C             On exit, the leading N-by-N part of this array contains
C             the transformed matrix.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     SCALE   (input) DOUBLE PRECISION array, dimension (N)
C             Details of the permutation and scaling factors, as
C             returned by DGEBAL.
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
C     Let P be a permutation matrix, and D a diagonal matrix of scaling
C     factors, both of order N. The routine computes
C                     -1
C        A <-- P D A D  P'.
C
C     where the permutation and scaling factors are encoded in the
C     array SCALE.
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C                               2
C     The algorithm requires O(N ) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05CY.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations.
C
C     ******************************************************************
C
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOB
      INTEGER           IGH, INFO, LDA, LOW, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), SCALE(*)
C     .. Local Scalars ..
      INTEGER           I, II, J, K
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DSCAL, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 )THEN
         INFO = -2
      ELSE IF( LOW.LT.1 .OR. LOW.GT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( IGH.LT.MIN( LOW, N ) .OR. IGH.GT.N ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB05OY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .OR. LSAME( JOB, 'N' ) )
     $   RETURN
C
      IF ( .NOT.LSAME( JOB, 'P' ) .AND. IGH.NE.LOW ) THEN
C
         DO 20 I = LOW, IGH
            CALL DSCAL( N, SCALE(I), A(I,1), LDA )
   20    CONTINUE
C
         DO 40 J = LOW, IGH
            CALL DSCAL( N, ONE/SCALE(J), A(1,J), 1 )
   40    CONTINUE
C
      END IF
C
      IF( .NOT.LSAME( JOB, 'S' ) ) THEN
C
         DO 60 II = 1, N
            I = II
            IF ( I.LT.LOW .OR. I.GT.IGH ) THEN
               IF ( I.LT.LOW ) I = LOW - II
               K = SCALE(I)
               IF ( K.NE.I ) THEN
                  CALL DSWAP( N, A(I,1), LDA, A(K,1), LDA )
                  CALL DSWAP( N, A(1,I), 1, A(1,K), 1 )
               END IF
            END IF
   60    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MB05OY ***
      END
