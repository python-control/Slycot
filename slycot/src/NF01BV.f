      SUBROUTINE NF01BV( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, J,
     $                   LDJ, JTJ, LDJTJ, DWORK, LDWORK, INFO )
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
C     To compute the matrix J'*J + c*I, for the Jacobian J as received
C     from SLICOT Library routine NF01BY, for one output variable.
C
C     NOTE: this routine must have the same arguments as SLICOT Library
C     routine NF01BU.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     STOR    CHARACTER*1
C             Specifies the storage scheme for the symmetric
C             matrix J'*J + c*I, as follows:
C             = 'F' :  full storage is used;
C             = 'P' :  packed storage is used.
C
C     UPLO    CHARACTER*1
C             Specifies which part of the matrix J'*J + c*I is stored,
C             as follows:
C             = 'U' :  the upper triagular part is stored;
C             = 'L' :  the lower triagular part is stored.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of columns of the Jacobian matrix J.  N >= 0.
C
C     IPAR    (input) INTEGER array, dimension (LIPAR)
C             The integer parameters describing the structure of the
C             matrix J, as follows:
C             IPAR(1) must contain the number of rows M of the Jacobian
C                     matrix J.  M >= 0.
C             IPAR is provided for compatibility with SLICOT Library
C             routine MD03AD.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 1.
C
C     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR)
C             The real parameters needed for solving the problem.
C             The entry DPAR(1) must contain the real scalar c.
C
C     LDPAR   (input) INTEGER
C             The length of the array DPAR.  LDPAR >= 1.
C
C     J       (input) DOUBLE PRECISION array, dimension (LDJ,N)
C             The leading M-by-N part of this array must contain the
C             Jacobian matrix J.
C
C     LDJ     INTEGER
C             The leading dimension of the array J.  LDJ >= MAX(1,M).
C
C     JTJ     (output) DOUBLE PRECISION array,
C                      dimension (LDJTJ,N),    if STOR = 'F',
C                      dimension (N*(N+1)/2),  if STOR = 'P'.
C             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if
C             STOR = 'P') part of this array contains the upper or
C             lower triangle of the matrix J'*J + c*I, depending on
C             UPLO = 'U', or UPLO = 'L', respectively, stored either as
C             a two-dimensional, or one-dimensional array, depending
C             on STOR.
C
C     LDJTJ   INTEGER
C             The leading dimension of the array JTJ.
C             LDJTJ >= MAX(1,N), if STOR = 'F'.
C             LDJTJ >= 1,        if STOR = 'P'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             Currently, this array is not used.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 0.
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
C     The matrix product is computed columnn-wise, exploiting the
C     symmetry. BLAS 3 routine DSYRK is used if STOR = 'F', and BLAS 2
C     routine DGEMV is used if STOR = 'P'.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001.
C
C     REVISIONS
C
C     V. Sima, March 2002.
C
C     KEYWORDS
C
C     Elementary matrix operations, matrix algebra, matrix operations,
C     Wiener system.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         STOR, UPLO
      INTEGER           INFO, LDJ, LDJTJ, LDPAR, LDWORK, LIPAR, N
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DPAR(*), DWORK(*), J(LDJ,*), JTJ(*)
C     .. Local Scalars ..
      LOGICAL           FULL, UPPER
      INTEGER           I, II, M
      DOUBLE PRECISION  C
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DLASET, DSYRK, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     ..
C     .. Executable Statements ..
C
      INFO  = 0
      FULL  = LSAME( STOR, 'F' )
      UPPER = LSAME( UPLO, 'U' )
C
      IF( .NOT.( FULL .OR. LSAME( STOR, 'P' ) ) ) THEN
         INFO = -1
      ELSEIF ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSEIF ( N.LT.0 ) THEN
         INFO = -3
      ELSEIF ( LIPAR.LT.1 ) THEN
         INFO = -5
      ELSEIF ( LDPAR.LT.1 ) THEN
         INFO = -7
      ELSEIF ( LDJTJ.LT.1 .OR. ( FULL .AND. LDJTJ.LT.N ) ) THEN
         INFO = -11
      ELSEIF ( LDWORK.LT.0 ) THEN
         INFO = -13
      ELSE
         M = IPAR(1)
         IF ( M.LT.0 ) THEN
            INFO = -4
         ELSEIF ( LDJ.LT.MAX( 1, M ) ) THEN
            INFO = -9
         ENDIF
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BV', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      C = DPAR(1)
      IF ( N.EQ.0 ) THEN
         RETURN
      ELSE IF ( M.EQ.0 ) THEN
         IF ( FULL ) THEN
            CALL DLASET( UPLO, N, N, ZERO, C, JTJ, LDJTJ )
         ELSE
            DUM(1) = ZERO
            CALL DCOPY( ( N*( N + 1 ) )/2, DUM, 0, JTJ, 1 )
            IF ( UPPER ) THEN
               II = 0
C
               DO 10 I = 1, N
                  II = II + I
                  JTJ(II) = C
   10          CONTINUE
C
            ELSE
               II = 1
C
               DO 20 I = N, 1, -1
                  JTJ(II) = C
                  II = II + I
   20          CONTINUE
C
            ENDIF
         ENDIF
         RETURN
      ENDIF
C
C     Build a triangle of the matrix J'*J + c*I.
C
      IF ( FULL ) THEN
         CALL DLASET( UPLO, N, N, ZERO, C, JTJ, LDJTJ )
         CALL DSYRK(  UPLO, 'Transpose', N, M, ONE, J, LDJ, ONE, JTJ,
     $                LDJTJ )
      ELSEIF ( UPPER ) THEN
         II = 0
C
         DO 30 I = 1, N
            CALL DGEMV( 'Transpose', M, I, ONE, J, LDJ, J(1,I), 1, ZERO,
     $                  JTJ(II+1), 1 )
            II = II + I
            JTJ(II) = JTJ(II) + C
   30    CONTINUE
C
      ELSE
         II = 1
C
         DO 40 I = N, 1, -1
            CALL DGEMV( 'Transpose', M, I, ONE, J(1,N-I+1), LDJ,
     $                  J(1,N-I+1), 1, ZERO, JTJ(II), 1 )
            JTJ(II) = JTJ(II) + C
            II = II + I
   40    CONTINUE
C
      ENDIF
C
      RETURN
C
C *** Last line of NF01BV ***
      END
