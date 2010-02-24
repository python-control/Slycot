      SUBROUTINE MB02TZ( NORM, N, HNORM, H, LDH, IPIV, RCOND, DWORK,
     $                   ZWORK, INFO )
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
C     To estimate the reciprocal of the condition number of a complex
C     upper Hessenberg matrix H, in either the 1-norm or the
C     infinity-norm, using the LU factorization computed by MB02SZ.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     NORM    CHARACTER*1
C             Specifies whether the 1-norm condition number or the
C             infinity-norm condition number is required:
C             = '1' or 'O':  1-norm;
C             = 'I':         Infinity-norm.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix H.  N >= 0.
C
C     HNORM   (input) DOUBLE PRECISION
C             If NORM = '1' or 'O', the 1-norm of the original matrix H.
C             If NORM = 'I', the infinity-norm of the original matrix H.
C
C     H       (input) COMPLEX*16 array, dimension (LDH,N)
C             The factors L and U from the factorization H = P*L*U
C             as computed by MB02SZ.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= max(1,N).
C
C     IPIV    (input) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= i <= N, row i of the matrix
C             was interchanged with row IPIV(i).
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal of the condition number of the matrix H,
C             computed as RCOND = 1/(norm(H) * norm(inv(H))).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     ZWORK   COMPLEX*16 array, dimension (2*N)
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
C     An estimate is obtained for norm(inv(H)), and the reciprocal of
C     the condition number is computed as
C        RCOND = 1 / ( norm(H) * norm(inv(H)) ).
C
C     REFERENCES
C
C     -
C
C     NUMERICAL ASPECTS
C                                2
C     The algorithm requires 0( N ) complex operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TB01FY by A.J. Laub, University of
C     Southern California, United States of America, May 1980.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2005.
C
C     KEYWORDS
C
C     Frequency response, Hessenberg form, matrix algebra.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDH, N
      DOUBLE PRECISION   HNORM, RCOND
C     ..
C     .. Array Arguments ..
      INTEGER            IPIV(*)
      DOUBLE PRECISION   DWORK( * )
      COMPLEX*16         H( LDH, * ), ZWORK( * )
C     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, J, JP, KASE, KASE1
C
      DOUBLE PRECISION   HINVNM, SCALE, SMLNUM
      COMPLEX*16         T, ZDUM
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, IZAMAX, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDRSCL, ZLACON, ZLATRS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
C     ..
C     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( HNORM.LT.ZERO ) THEN
         INFO = -3
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02TZ', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( HNORM.EQ.ZERO ) THEN
         RETURN
      END IF
C
      SMLNUM = DLAMCH( 'Safe minimum' )
C
C     Estimate the norm of inv(H).
C
      HINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL ZLACON( N, ZWORK( N+1 ), ZWORK, HINVNM, KASE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
C
C           Multiply by inv(L).
C
            DO 20 J = 1, N - 1
               JP = IPIV( J )
               T = ZWORK( JP )
               IF( JP.NE.J ) THEN
                  ZWORK( JP ) = ZWORK( J )
                  ZWORK( J ) = T
               END IF
               ZWORK( J+1 ) = ZWORK( J+1 ) - T * H( J+1, J )
   20       CONTINUE
C
C           Multiply by inv(U).
C
            CALL ZLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   H, LDH, ZWORK, SCALE, DWORK, INFO )
         ELSE
C
C           Multiply by inv(U').
C
            CALL ZLATRS( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                   NORMIN, N, H, LDH, ZWORK, SCALE, DWORK, INFO )
C
C           Multiply by inv(L').
C
            DO 30 J = N - 1, 1, -1
               ZWORK( J ) = ZWORK( J ) -
     $                      DCONJG( H( J+1, J ) ) * ZWORK( J+1 )
               JP = IPIV( J )
               IF( JP.NE.J ) THEN
                  T = ZWORK( JP )
                  ZWORK( JP ) = ZWORK( J )
                  ZWORK( J ) = T
               END IF
   30       CONTINUE
         END IF
C
C        Divide X by 1/SCALE if doing so will not cause overflow.
C
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IZAMAX( N, ZWORK, 1 )
            IF( SCALE.LT.CABS1( ZWORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO
     $        ) GO TO 40
            CALL ZDRSCL( N, SCALE, ZWORK, 1 )
         END IF
         GO TO 10
      END IF
C
C     Compute the estimate of the reciprocal condition number.
C
      IF( HINVNM.NE.ZERO )
     $   RCOND = ( ONE / HINVNM ) / HNORM
C
   40 CONTINUE
      RETURN
C *** Last line of MB02TZ ***
      END
